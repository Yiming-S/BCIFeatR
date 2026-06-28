suppressMessages({library(gsignal); library(R.matlab); library(edfReader)})
suppressMessages(devtools::load_all("/Users/yiming/Documents/GitHub/BCIFeatR", quiet = TRUE))
source("/Users/yiming/Documents/GitHub/Data-Analytics-Lab-Prof.Degras/SSMs_BCI/code/bci_clf/eegDataLoaders.R")
OUT <- "/private/tmp/claude-501/-Users-yiming-Agent-PhD-working/84c1c68b-ba94-4560-816a-ea07efe156e5/scratchpad"
ZHOU <- "/Volumes/Yiming SSD Drive/EEG_MotorImagery/MNE-zhou-2016"
BNCI <- "/Volumes/Yiming SSD Drive/EEG_MotorImagery/BNCI2014_004/MNE-bnci-data/database/data-sets/004-2014"

## ---- neutral downstream classifier: shrinkage LDA (reuses LW_covariance) ----
lda_fit <- function(X, y) {
  y <- factor(y); classes <- levels(y)
  mu <- colMeans(X); sdv <- apply(X, 2, sd); sdv[!is.finite(sdv) | sdv <= 0] <- 1
  Z <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  resid <- Z; cmeans <- matrix(0, length(classes), ncol(Z), dimnames = list(classes, NULL))
  for (c in classes) { idx <- y == c; cmeans[c, ] <- colMeans(Z[idx, , drop = FALSE]); resid[idx, ] <- sweep(Z[idx, , drop = FALSE], 2, cmeans[c, ], "-") }
  Sw <- LW_covariance(resid); Swinv <- solve(Sw + diag(1e-6, ncol(Sw)))
  list(classes = classes, mu = mu, sdv = sdv, cmeans = cmeans, Swinv = Swinv,
       logprior = log(as.numeric(table(y)) / length(y)))
}
lda_predict <- function(m, X) {
  Z <- sweep(sweep(X, 2, m$mu, "-"), 2, m$sdv, "/")
  M <- m$cmeans %*% m$Swinv; quad <- rowSums(M * m$cmeans)
  score <- sweep(Z %*% t(M), 2, 0.5 * quad - m$logprior, "-")
  factor(m$classes[max.col(score, ties.method = "first")], levels = m$classes)
}
bal_acc <- function(truth, pred) {
  truth <- factor(truth); pred <- factor(pred, levels = levels(truth))
  mean(vapply(levels(truth), function(c) { i <- truth == c; if (!any(i)) NA_real_ else mean(pred[i] == c) }, numeric(1)), na.rm = TRUE)
}

## ---- methods ----
feat_methods <- function(nch, fs) {
  bands <- list(c(8, 12), c(12, 16), c(16, 24), c(24, 30))
  list(
    logvar     = list(f = "logvar",     p = list()),
    CSP        = list(f = "CSP",        p = list(ncomps = min(6L, nch))),
    FBCSP      = list(f = "FBCSP",      p = list(fs = fs, channels = seq_len(nch), frequency_bands = bands, ncomps = 2L)),
    bandpower  = list(f = "bandpower",  p = list(fs = fs, frequency_bands = bands)),
    Hjorth     = list(f = "Hjorth",     p = list()),
    TS         = list(f = "TS",         p = list(cov_type = "oas")),
    ACM_TS     = list(f = "ACM_TS",     p = list(order = 2L, delay = 1L, shrinkage = "oas", use_filter = FALSE)),
    Riemannian = list(f = "Riemannian", p = list(cov_type = "oas", use_filter = FALSE)),
    MVAR       = list(f = "MVAR",       p = list(order = 3L))
  )
}

evaluate_fold <- function(tr, te) {
  tr$y <- factor(tr$y); te$y <- factor(te$y, levels = levels(tr$y))
  nch <- ncol(tr$x[[1]]); fs <- tr$fs
  methods <- feat_methods(nch, fs)
  acc <- setNames(rep(NA_real_, length(methods) + 2L), c(names(methods), "MDM", "CSP_MSVAR"))
  tm  <- acc
  for (nm in names(methods)) {
    g <- tryCatch({
      t0 <- Sys.time()
      fit <- featEx4Train(tr$x, tr$y, methods[[nm]]$f, methods[[nm]]$p)
      Xte <- featEx4Test(te$x, fit$object, methods[[nm]]$f)
      a <- bal_acc(te$y, lda_predict(lda_fit(fit$features, tr$y), Xte))
      list(a = a, t = as.numeric(Sys.time() - t0, units = "secs"))
    }, error = function(e) list(a = NA_real_, t = NA_real_))
    acc[nm] <- g$a; tm[nm] <- g$t
  }
  acc["MDM"] <- tryCatch(bal_acc(te$y, mdm_predict(mdm_train(tr$x, tr$y, cov_type = "oas", metric = "logeuclid"), te$x)),
                         error = function(e) NA_real_)
  acc["CSP_MSVAR"] <- tryCatch({
    cf <- featEx4Train(tr$x, tr$y, "CSP", list(ncomps = min(4L, ncol(tr$x[[1]]))))
    W <- cf$object$filter
    ptr <- lapply(tr$x, function(M) M %*% W); pte <- lapply(te$x, function(M) M %*% W)
    mm <- msvar_train(ptr, tr$y, M = 2L, p = 1L, seed = 1L,
                      control = list(maxit = 8L, tol = 1e-4, init = "conditional", nseg = 20L))
    bal_acc(te$y, predict(mm, pte))
  }, error = function(e) NA_real_)
  list(acc = acc, tm = tm)
}

cat("=== Loading + evaluating zhou2016 (leave-one-session-out, binary L/R) ===\n")
zhou_rows <- list(); zhou_tm <- list()
for (sj in 1:4) {
  sess <- lapply(1:3, function(k) tryCatch(load_eeg_zhou2016(subject = as.character(sj), sess_index = k, data_dir = ZHOU, nclass = 2L, verbose = FALSE), error = function(e) NULL))
  sess <- Filter(Negate(is.null), sess)
  if (length(sess) < 2L) { cat("  sub", sj, "insufficient sessions\n"); next }
  for (ho in seq_along(sess)) {
    te <- sess[[ho]]; trs <- sess[-ho]
    tr <- list(x = do.call(c, lapply(trs, `[[`, "x")), y = unlist(lapply(trs, `[[`, "y")), fs = trs[[1]]$fs)
    r <- evaluate_fold(tr, te)
    zhou_rows[[length(zhou_rows) + 1L]] <- data.frame(subject = paste0("sub", sj), fold = ho, t(r$acc))
    zhou_tm[[length(zhou_tm) + 1L]] <- r$tm
    cat(sprintf("  sub%d fold%d done (n_tr=%d n_te=%d)\n", sj, ho, length(tr$x), length(te$x)))
  }
}
zhou <- do.call(rbind, zhou_rows)
write.csv(zhou, file.path(OUT, "real_zhou.csv"), row.names = FALSE)

cat("=== Loading + evaluating BNCI2014_004 (train T, test E, binary L/R) ===\n")
bnci_rows <- list(); bnci_tm <- list()
for (sj in sprintf("B%02d", 1:9)) {
  tr <- tryCatch(load_eeg_bnci004(subject = sj, sess_index = 1:3, data_dir = BNCI, nclass = 2L), error = function(e) NULL)
  te <- tryCatch(load_eeg_bnci004(subject = sj, sess_index = 4:5, data_dir = BNCI, nclass = 2L), error = function(e) NULL)
  if (is.null(tr) || is.null(te)) { cat("  ", sj, "load failed\n"); next }
  r <- evaluate_fold(tr, te)
  bnci_rows[[length(bnci_rows) + 1L]] <- data.frame(subject = sj, fold = 1L, t(r$acc))
  bnci_tm[[length(bnci_tm) + 1L]] <- r$tm
  cat(sprintf("  %s done (n_tr=%d n_te=%d)\n", sj, length(tr$x), length(te$x)))
}
bnci <- do.call(rbind, bnci_rows)
write.csv(bnci, file.path(OUT, "real_bnci.csv"), row.names = FALSE)

summ <- function(df) {
  cols <- setdiff(names(df), c("subject", "fold"))
  per_subj <- aggregate(df[cols], list(subject = df$subject), mean, na.rm = TRUE)
  data.frame(method = cols,
             mean = round(sapply(cols, function(c) mean(per_subj[[c]], na.rm = TRUE)), 3),
             sd   = round(sapply(cols, function(c) sd(per_subj[[c]], na.rm = TRUE)), 3))
}
cat("\n===== ZHOU2016 summary (mean +/- sd balanced acc across 4 subjects) =====\n")
zs <- summ(zhou); print(zs[order(-zs$mean), ], row.names = FALSE)
cat("\n===== BNCI2014_004 summary (mean +/- sd across 9 subjects) =====\n")
bs <- summ(bnci); print(bs[order(-bs$mean), ], row.names = FALSE)
cat("\nmean fold time (s) zhou:", round(mean(sapply(zhou_tm, sum), na.rm=TRUE),1), " bnci:", round(mean(sapply(bnci_tm, sum), na.rm=TRUE),1), "\n")
saveRDS(list(zhou=zhou, bnci=bnci, zhou_tm=zhou_tm, bnci_tm=bnci_tm), file.path(OUT, "real_results.rds"))
cat("DONE\n")
