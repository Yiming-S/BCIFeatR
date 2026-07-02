# Within-session benchmark: stratified 5-fold CV *inside each session*.
# Mirrors realdata_bench.R (same methods, same shrinkage-LDA) so the numbers are
# directly comparable to the cross-session results -> quantifies session drift.
suppressMessages({library(gsignal); library(R.matlab); library(edfReader)})
suppressMessages(devtools::load_all("/Users/yiming/Documents/GitHub/BCIFeatR", quiet = TRUE))
source("/Users/yiming/Documents/GitHub/Data-Analytics-Lab-Prof.Degras/SSMs_BCI/code/bci_clf/eegDataLoaders.R")
OUT <- "/private/tmp/claude-501/-Users-yiming-Agent-PhD-working/84c1c68b-ba94-4560-816a-ea07efe156e5/scratchpad"
ZHOU <- "/Volumes/Yiming SSD Drive/EEG_MotorImagery/MNE-zhou-2016"
BNCI <- "/Volumes/Yiming SSD Drive/EEG_MotorImagery/BNCI2014_004/MNE-bnci-data/database/data-sets/004-2014"

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
feat_methods <- function(nch, fs) {
  bands <- list(c(8, 12), c(12, 16), c(16, 24), c(24, 30))
  list(
    logvar = list(f="logvar", p=list()),
    CSP = list(f="CSP", p=list(ncomps=min(6L,nch))),
    FBCSP = list(f="FBCSP", p=list(fs=fs, channels=seq_len(nch), frequency_bands=bands, ncomps=2L)),
    bandpower = list(f="bandpower", p=list(fs=fs, frequency_bands=bands)),
    Hjorth = list(f="Hjorth", p=list()),
    TS = list(f="TS", p=list(cov_type="oas")),
    ACM_TS = list(f="ACM_TS", p=list(order=2L, delay=1L, shrinkage="oas", use_filter=FALSE)),
    Riemannian = list(f="Riemannian", p=list(cov_type="oas", use_filter=FALSE)),
    MVAR = list(f="MVAR", p=list(order=3L))
  )
}
evaluate_fold <- function(tr, te) {
  tr$y <- factor(tr$y); te$y <- factor(te$y, levels = levels(tr$y))
  nch <- ncol(tr$x[[1]]); fs <- tr$fs; methods <- feat_methods(nch, fs)
  acc <- setNames(rep(NA_real_, length(methods) + 1L), c(names(methods), "MDM"))
  for (nm in names(methods)) {
    acc[nm] <- tryCatch({
      fit <- featEx4Train(tr$x, tr$y, methods[[nm]]$f, methods[[nm]]$p)
      bal_acc(te$y, lda_predict(lda_fit(fit$features, tr$y), featEx4Test(te$x, fit$object, methods[[nm]]$f)))
    }, error = function(e) NA_real_)
  }
  acc["MDM"] <- tryCatch(bal_acc(te$y, mdm_predict(mdm_train(tr$x, tr$y, cov_type="oas", metric="logeuclid"), te$x)), error=function(e) NA_real_)
  acc
}
strat_folds <- function(y, k = 5L, seed = 1L) {
  set.seed(seed); y <- factor(y); folds <- integer(length(y))
  for (lv in levels(y)) { idx <- sample(which(y == lv)); folds[idx] <- rep_len(seq_len(k), length(idx)) }
  folds
}
# within-session CV on one loaded session
within_cv <- function(sess, k = 5L) {
  y <- factor(sess$y)
  if (nlevels(y) < 2L || min(table(y)) < k) return(NULL)
  folds <- strat_folds(y, k)
  accs <- lapply(seq_len(k), function(f) {
    tr <- list(x = sess$x[folds != f], y = y[folds != f], fs = sess$fs)
    te <- list(x = sess$x[folds == f], y = y[folds == f], fs = sess$fs)
    evaluate_fold(tr, te)
  })
  colMeans(do.call(rbind, accs), na.rm = TRUE)   # mean over folds
}

cat("=== zhou2016 within-session (5-fold CV inside each session) ===\n")
zrows <- list()
for (sj in 1:4) for (k in 1:3) {
  sess <- tryCatch(load_eeg_zhou2016(subject=as.character(sj), sess_index=k, data_dir=ZHOU, nclass=2L, verbose=FALSE), error=function(e) NULL)
  if (is.null(sess)) next
  a <- within_cv(sess); if (is.null(a)) next
  zrows[[length(zrows)+1L]] <- data.frame(subject=paste0("sub",sj), session=k, t(a))
  cat(sprintf("  sub%d ses%d done (n=%d)\n", sj, k, length(sess$x)))
}
zhou <- do.call(rbind, zrows); write.csv(zhou, file.path(OUT, "within_zhou.csv"), row.names=FALSE)

cat("=== BNCI2014_004 within-session (5-fold CV inside each of 5 sessions) ===\n")
brows <- list()
for (sj in sprintf("B%02d", 1:9)) for (k in 1:5) {
  sess <- tryCatch(load_eeg_bnci004(subject=sj, sess_index=k, data_dir=BNCI, nclass=2L), error=function(e) NULL)
  if (is.null(sess)) next
  a <- within_cv(sess); if (is.null(a)) next
  brows[[length(brows)+1L]] <- data.frame(subject=sj, session=k, t(a))
  cat(sprintf("  %s ses%d done (n=%d)\n", sj, k, length(sess$x)))
}
bnci <- do.call(rbind, brows); write.csv(bnci, file.path(OUT, "within_bnci.csv"), row.names=FALSE)

summ <- function(df) {
  cols <- setdiff(names(df), c("subject","session"))
  per_subj <- aggregate(df[cols], list(subject=df$subject), mean, na.rm=TRUE)
  data.frame(method=cols, mean=round(sapply(cols, function(c) mean(per_subj[[c]],na.rm=TRUE)),3),
             sd=round(sapply(cols, function(c) sd(per_subj[[c]],na.rm=TRUE)),3))
}
cat("\n===== ZHOU2016 within-session summary =====\n"); zs<-summ(zhou); print(zs[order(-zs$mean),], row.names=FALSE)
cat("\n===== BNCI2014_004 within-session summary =====\n"); bs<-summ(bnci); print(bs[order(-bs$mean),], row.names=FALSE)
saveRDS(list(zhou=zhou, bnci=bnci), file.path(OUT, "within_results.rds"))
cat("DONE\n")
