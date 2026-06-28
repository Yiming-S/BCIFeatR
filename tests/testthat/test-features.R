make_session <- function(n_trials = 16L, n_time = 128L, n_ch = 8L, seed = 1L) {
  set.seed(seed)
  y <- factor(rep(1:2, length.out = n_trials))
  x <- lapply(seq_len(n_trials), function(i) {
    cls <- as.integer(y[i])
    t <- seq(0, 1, length.out = n_time)
    base <- matrix(rnorm(n_time * n_ch, sd = 0.8), nrow = n_time, ncol = n_ch)
    base[, 1] <- base[, 1] + sin(2 * pi * (8 + cls) * t) * (0.8 + 0.2 * cls)
    base[, 2] <- base[, 2] + cos(2 * pi * 12 * t) * (0.5 + 0.1 * cls)
    base
  })
  list(x = x, y = y)
}

run_one_feature <- function(feature, params, train_sess, test_sess) {
  fit <- featEx4Train(train_sess$x, train_sess$y, feature, params)
  x_train <- fit$features
  x_test <- featEx4Test(test_sess$x, fit$object, feature)

  expect_true(is.matrix(x_train))
  expect_true(is.matrix(x_test))
  expect_equal(nrow(x_train), length(train_sess$x))
  expect_equal(nrow(x_test), length(test_sess$x))
  expect_equal(ncol(x_train), ncol(x_test))
  expect_false(anyNA(x_train))
  expect_false(anyNA(x_test))
  expect_true(is.list(fit$object))
  expect_true(is.list(fit$object$metadata))
  expect_identical(fit$object$metadata$feature, feature)
  expect_identical(fit$object$metadata$n_features, ncol(x_train))
  list(fit = fit, x_train = x_train, x_test = x_test)
}

train_sess <- make_session(n_trials = 16L, seed = 7L)
test_sess <- make_session(n_trials = 8L, seed = 77L)

feature_cases <- list(
  logvar = list(params = list()),
  logvar_pca = list(params = list(ncomps = 4L)),
  CSP = list(params = list(ncomps = 4L)),
  FBCSP = list(params = list(
    ncomps = 2L, fs = 128L, channels = 1:8,
    frequency_bands = list(c(8, 12), c(12, 20))
  )),
  FBCSSP = list(params = list(
    ncomps = 2L, fs = 128L, channels = 1:8,
    frequency_bands = list(c(8, 12), c(12, 20))
  )),
  TS = list(params = list(cov_type = "oas")),
  ACM_TS = list(params = list(
    order = 2L, delay = 1L, shrinkage = "oas", use_filter = FALSE
  )),
  Riemannian = list(params = list(
    cov_type = "oas", use_filter = FALSE, jitter_sd = 1e-8
  )),
  ATM = list(params = list(
    z = 2.5, bin = 1L, tuneBin = FALSE, minLen = 2L,
    useAbs = TRUE, includeDiag = TRUE, binCandidates = c(1L, 2L, 3L)
  )),
  bandpower = list(params = list(
    fs = 128L, frequency_bands = list(c(8, 12), c(12, 30))
  )),
  bandpower_relative = list(feature = "bandpower", params = list(
    fs = 128L, frequency_bands = list(c(8, 12), c(12, 30)), relative = TRUE
  )),
  Hjorth = list(params = list()),
  Hjorth_log = list(feature = "Hjorth", params = list(log_activity = TRUE)),
  MVAR = list(params = list(order = 2L)),
  MVAR_noQ = list(feature = "MVAR", params = list(order = 2L, include_Q = FALSE)),
  logvar_simam = list(feature = "logvar", params = list(simam = "robust")),
  TS_simam = list(feature = "TS", params = list(cov_type = "oas", simam = TRUE)),
  MSVAR = list(params = list(
    M = 2L, p = 1L, seed = 1L,
    control = list(maxit = 3L, tol = 1e-4, init = "conditional", nseg = 15L)
  ))
)

for (case_name in names(feature_cases)) {
  local({
    case <- feature_cases[[case_name]]
    feat <- if (!is.null(case$feature)) case$feature else case_name
    test_that(sprintf("%s train/test roundtrip", case_name), {
      run_one_feature(feat, case$params, train_sess, test_sess)
    })
  })
}

test_that("session-list input contract", {
  multi_fit <- featEx4Train(
    x = list(train_sess$x, test_sess$x),
    y = list(train_sess$y, test_sess$y),
    feature = "logvar", params = list()
  )
  expect_true(is.list(multi_fit) && length(multi_fit) == 2L)
  expect_true(is.matrix(multi_fit[[1]]$features))
  expect_true(is.matrix(multi_fit[[2]]$features))

  single_fit <- featEx4Train(
    x = list(train_sess$x), y = list(train_sess$y),
    feature = "logvar", params = list()
  )
  expect_true(is.list(single_fit) && length(single_fit) == 1L)
  expect_true(is.matrix(single_fit[[1]]$features))
})

test_that("ATM session-list with y = NULL", {
  atm_multi_fit <- featEx4Train(
    x = list(train_sess$x, test_sess$x), y = NULL,
    feature = "ATM", params = list(bin = 1L, tuneBin = FALSE)
  )
  expect_true(is.list(atm_multi_fit) && length(atm_multi_fit) == 2L)
  expect_true(all(vapply(atm_multi_fit, function(z) is.matrix(z$features), logical(1))))
})

test_that("Riemannian reproducibility", {
  riem_params <- list(cov_type = "oas", use_filter = FALSE, jitter_sd = 1e-8)
  fit1 <- featEx4Train(train_sess$x, train_sess$y, "Riemannian", riem_params)
  fit2 <- featEx4Train(train_sess$x, train_sess$y, "Riemannian", riem_params)
  expect_equal(fit1$features, fit2$features, tolerance = 0)

  test1 <- featEx4Test(test_sess$x, fit1$object, "Riemannian")
  test2 <- featEx4Test(test_sess$x, fit1$object, "Riemannian")
  expect_equal(test1, test2, tolerance = 0)
})

test_that("preprocessing center/scale/whiten", {
  for (pp in c("center", "scale", "whiten")) {
    fit_pp <- featEx4Train(
      train_sess$x, train_sess$y, "TS",
      params = list(cov_type = "oas", preprocess = pp)
    )
    pred_pp <- featEx4Test(test_sess$x, fit_pp$object, "TS")
    expect_identical(fit_pp$object$preprocess$method, pp)
    expect_equal(ncol(fit_pp$features), ncol(pred_pp))
    expect_false(anyNA(pred_pp))
  }
})

test_that("input contract validation", {
  expect_error(featEx4Train(train_sess$x, train_sess$y[-1], "logvar", list()))

  x_bad <- train_sess$x
  x_bad[[2]] <- x_bad[[2]][, 1:7, drop = FALSE]
  expect_error(featEx4Train(x_bad, train_sess$y, "logvar", list()))

  x_na <- train_sess$x
  x_na[[1]][1, 1] <- NA_real_
  expect_error(featEx4Train(x_na, train_sess$y, "logvar", list()))
})

test_that("feature-specific parameter validation", {
  expect_error(featEx4Train(train_sess$x, train_sess$y, "logvar_pca", list()))
  expect_error(featEx4Train(train_sess$x, train_sess$y, "FBCSP",
                            list(ncomps = 2L, fs = 128L)))
  acm_def <- featEx4Train(train_sess$x, train_sess$y, "ACM_TS",
                           params = list(use_filter = FALSE))
  expect_true(is.matrix(acm_def$features))
})

test_that("degenerate input stability", {
  expect_true(all(is.finite(LW_covariance(matrix(1, 128, 8)))))

  deg <- make_session(n_trials = 8L, seed = 99L)
  deg$x <- lapply(deg$x, function(M) { M[, 1] <- 0; M })
  deg_fit <- featEx4Train(deg$x, deg$y, "FBCSP", params = list(
    ncomps = 2L, fs = 128L, channels = 1:8,
    frequency_bands = list(c(8, 12), c(12, 20))
  ))
  deg_test <- featEx4Test(test_sess$x, deg_fit$object, "FBCSP")
  expect_true(all(is.finite(deg_fit$features)))
  expect_true(all(is.finite(deg_test)))
})

test_that("classifier and feature-selection contracts", {
  expect_error(fisher_lda(matrix(rnorm(4 * 9), nrow = 4), rep(1:3, each = 3)))

  fs_data <- list(
    matrix(rnorm(30), nrow = 10, ncol = 3),
    matrix(rnorm(36), nrow = 12, ncol = 3)
  )
  f_sess <- fisher(fs_data, labels = c("A", "B"), k = 0.5)
  f_samp <- fisher(
    fs_data,
    labels = c(rep("A", nrow(fs_data[[1]])), rep("B", nrow(fs_data[[2]]))),
    k = 0.5
  )
  expect_equal(nrow(f_sess$fisher_df), nrow(f_samp$fisher_df))

  expect_error(pca(fs_data, k = 1.5))
  expect_error(fisher(fs_data, labels = c("A", "B"), k = 1.5))
  expect_error(pca(fs_data, method = "random", num_random = 0, k = 0.5))
  expect_error(fisher(fs_data, labels = c("A", "B"), method = "random",
                      num_random = 0, k = 0.5))
})

test_that("edge-case guardrails", {
  expect_error(compute_ACM(train_sess$x, order = 0, delay = 1, shrinkage = "oas"))
  expect_error(compute_ACM(train_sess$x, order = 2, delay = 0, shrinkage = "oas"))

  atm_robust <- atmTrain(
    train_sess$x,
    params = list(tuneBin = TRUE, bin = 1L, binCandidates = c(1L, 999L))
  )
  expect_true(all(is.finite(atm_robust$features)))

  P_filt <- geodesic_filtering(diag(3), diag(3),
                                matrix(0, nrow = choose(4, 2), ncol = 2))
  expect_true(all(is.finite(P_filt)))
})

test_that("serialization consistency", {
  riem_params <- list(cov_type = "oas", use_filter = FALSE, jitter_sd = 1e-8)
  fit1 <- featEx4Train(train_sess$x, train_sess$y, "Riemannian", riem_params)
  test1 <- featEx4Test(test_sess$x, fit1$object, "Riemannian")

  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit1$object, tmp)
  obj_loaded <- readRDS(tmp)
  test_loaded <- featEx4Test(test_sess$x, obj_loaded, "Riemannian")
  expect_equal(test1, test_loaded, tolerance = 0)
})

test_that("bandpower feature dimensions and naming", {
  bands <- list(c(8, 12), c(12, 30), c(13, 28))
  fit <- featEx4Train(train_sess$x, train_sess$y, "bandpower",
                      params = list(fs = 128L, frequency_bands = bands))
  n_ch <- ncol(train_sess$x[[1]])
  expect_equal(ncol(fit$features), n_ch * length(bands))
  expect_true(all(grepl("^bp_b", colnames(fit$features))))

  rel <- featEx4Train(train_sess$x, train_sess$y, "bandpower",
                      params = list(fs = 128L, frequency_bands = bands, relative = TRUE))
  # Relative band power is a within-channel fraction across bands -> rows of the
  # per-channel band blocks sum to 1.
  M <- rel$features[1, ]
  per_ch <- matrix(M, nrow = n_ch, ncol = length(bands))
  expect_equal(unname(rowSums(per_ch)), rep(1, n_ch), tolerance = 1e-8)
})

test_that("Hjorth feature dimensions and finiteness", {
  fit <- featEx4Train(train_sess$x, train_sess$y, "Hjorth", params = list())
  n_ch <- ncol(train_sess$x[[1]])
  expect_equal(ncol(fit$features), 3L * n_ch)
  expect_true(all(is.finite(fit$features)))
  expect_true(all(grepl("^hjorth_(activity|mobility|complexity)_", colnames(fit$features))))
})

test_that("MDM train/predict roundtrip", {
  for (metric in c("logeuclid", "euclid")) {
    m <- mdm_train(train_sess$x, train_sess$y, cov_type = "oas", metric = metric)
    expect_s3_class(m, "mdm")
    pred <- mdm_predict(m, test_sess$x)
    expect_true(is.factor(pred))
    expect_length(pred, length(test_sess$x))
    expect_setequal(levels(pred), levels(train_sess$y))
    D <- attr(pred, "distances")
    expect_equal(dim(D), c(length(test_sess$x), nlevels(train_sess$y)))
    expect_false(anyNA(D))
  }
  # Predicted label is the nearest class mean.
  m <- mdm_train(train_sess$x, train_sess$y, cov_type = "oas", metric = "logeuclid")
  pred <- mdm_predict(m, test_sess$x)
  D <- attr(pred, "distances")
  expect_equal(as.character(pred),
               m$classes[max.col(-D, ties.method = "first")])
})

test_that("vectorized atmCounts matches the reference double loop", {
  set.seed(3)
  B <- matrix(as.integer(runif(40 * 6) > 0.5), nrow = 40, ncol = 6)
  segs <- atmFindAvalanches(B, minLen = 2L)
  got <- atmCounts(B, segs)

  D <- ncol(B)
  cij <- matrix(0L, D, D); ci <- integer(D)
  for (se in segs) {
    s <- se[1]; e <- se[2]
    if (e <= s) next
    for (t in s:(e - 1L)) {
      ii <- which(B[t, ] != 0L)
      if (!length(ii)) next
      jj <- which(B[t + 1L, ] != 0L)
      if (length(jj)) cij[ii, jj] <- cij[ii, jj] + 1L
      ci[ii] <- ci[ii] + 1L
    }
  }
  expect_equal(unname(got$countIJ), matrix(as.numeric(cij), D, D))
  expect_equal(unname(got$countI), as.numeric(ci))
})

test_that("mibif matches the reference greedy max-relevance selection", {
  skip_if_not_installed("infotheo")
  set.seed(11)
  n <- 30L
  mk <- function(shift) {
    M <- matrix(rnorm(n * 5), n, 5)
    M[, 2] <- M[, 2] + shift * 3
    M[, 4] <- M[, 4] + shift * 2
    colnames(M) <- paste0("f", 1:5)
    M
  }
  data_list <- lapply(c(0, 1, 2), mk)
  labs <- c(1L, 2L, 3L)

  # Reference: the original per-iteration greedy max-relevance loop. The fast
  # path computes each feature's MI once instead, which must give the same rank.
  cv <- unlist(mapply(rep, labs, vapply(data_list, nrow, integer(1))))
  cd <- do.call(rbind, data_list)
  k  <- min(ceiling(ncol(cd) * 0.4), ncol(cd))
  cm <- function(f) infotheo::mutinformation(infotheo::discretize(f),
                                             infotheo::discretize(cv))
  sel <- logical(ncol(cd)); ref <- integer(0)
  while (length(ref) < k) {
    mv <- vapply(seq_len(ncol(cd)),
                 function(i) if (sel[i]) -Inf else cm(cd[, i]), numeric(1))
    b <- which.max(mv); sel[b] <- TRUE; ref <- c(ref, b)
  }

  res <- mibif(data_list, labels = labs, k = 0.4)
  expect_equal(res$selected_indices, ref)
  expect_equal(length(res$selected_indices), k)
  expect_equal(colnames(res$selected_features), colnames(cd)[res$selected_indices])
  # method = "combine" is deterministic.
  expect_equal(mibif(data_list, labels = labs, k = 0.4)$selected_indices,
               res$selected_indices)
})

test_that("splitTimeRange segments by label runs and fixed windows", {
  lab <- c(rep(1L, 5), rep(2L, 3), rep(1L, 4))
  tr <- splitTimeRange(lab, by_label = TRUE)
  expect_length(tr, 3L)
  expect_equal(vapply(tr, length, integer(1)), c(5L, 3L, 4L))
  expect_equal(tr[[1]], 1:5)
  expect_equal(tr[[2]], 6:8)
  expect_equal(unlist(tr), seq_along(lab))           # full, ordered coverage

  # Fixed windows within a single run; trailing short window is dropped.
  win <- splitTimeRange(rep(1L, 10), by_label = FALSE, seglen = 4L)
  expect_length(win, 2L)
  expect_true(all(vapply(win, length, integer(1)) == 4L))
  expect_error(splitTimeRange(rep(1L, 10), by_label = FALSE))  # seglen required

  # snake_case alias matches.
  expect_identical(split_time_range(lab), tr)
})

test_that("freqBandSelect recovers a discriminative band", {
  set.seed(1)
  fs <- 128L; L <- 256L; nblk <- 6L; nc <- 4L
  blocks <- lapply(seq_len(nblk), function(b) {
    cls <- if (b %% 2L == 1L) 1L else 2L
    t <- seq(0, (L - 1) / fs, length.out = L)
    M <- matrix(rnorm(L * nc, sd = 0.5), L, nc)
    M[, 1] <- M[, 1] + 2 * sin(2 * pi * (if (cls == 1L) 10 else 22) * t)
    list(x = M, y = rep(cls, L))
  })
  x <- do.call(rbind, lapply(blocks, `[[`, "x"))
  y <- unlist(lapply(blocks, `[[`, "y"))

  fsel <- freqBandSelect(x, y, fs = fs, flo = 4, fhi = 40, winlen = 256L,
                         by_label = TRUE)
  expect_length(fsel$trials, nblk)
  expect_equal(length(fsel$freqrange), 2L)
  expect_true(fsel$freqrange[1] <= fsel$freqrange[2])
  expect_true(all(is.finite(fsel$freqrange)))
  # The planted contrast sits near 22 Hz (vs 10 Hz baseline); the selected band
  # should land in that neighbourhood, well away from 10 Hz.
  expect_lte(abs(mean(fsel$freqrange) - 22), 4)
})

test_that("simAM variants preserve shape and structure", {
  M <- train_sess$x[[1]]
  for (m in c("vanilla", "robust", "corr")) {
    W <- simAM(M, method = m)
    expect_equal(dim(W), dim(M))
    expect_false(anyNA(W))
  }
  lst <- simAM(train_sess$x, method = "vanilla")
  expect_true(is.list(lst) && length(lst) == length(train_sess$x))
  expect_true(is.matrix(lst[[1]]))
})

test_that("simAM robust default is auto-scaled and does not crush the signal", {
  set.seed(7)
  M <- matrix(rnorm(400 * 8), 400, 8)
  wr <- simAM(M, method = "robust")     # default lambda -> auto-scaled spread
  wv <- simAM(M, method = "vanilla")
  # Robust must retain a substantial fraction of the signal (the old lambda=0.1
  # mis-scaling collapsed weights toward 0, retaining ~0.1% of the energy).
  expect_gt(sd(wr) / sd(M), 0.3)
  # And it should scale comparably to vanilla on clean Gaussian data.
  expect_lt(abs(log(sd(wr) / sd(wv))), 0.7)
})

test_that("simAM attention integrates with the dispatcher deterministically", {
  fit <- featEx4Train(train_sess$x, train_sess$y, "logvar",
                      params = list(simam = TRUE))
  expect_true(isTRUE(fit$object$simam$enabled))
  p1 <- featEx4Test(test_sess$x, fit$object, "logvar")
  p2 <- featEx4Test(test_sess$x, fit$object, "logvar")
  expect_equal(p1, p2, tolerance = 0)

  base <- featEx4Train(train_sess$x, train_sess$y, "logvar", params = list())
  # Attention re-weighting changes the resulting features.
  expect_false(isTRUE(all.equal(unname(fit$features), unname(base$features))))
  # And the disabled default leaves features untouched.
  expect_equal(featEx4Train(train_sess$x, train_sess$y, "logvar",
                            params = list(simam = FALSE))$features,
               base$features)
})

test_that("fit_var recovers a known VAR(1) process", {
  set.seed(2)
  N <- 3000L
  A_true <- matrix(c(0.5, 0.2, -0.3, 0.4), 2, 2)
  y <- matrix(0, N, 2)
  for (t in 2:N) y[t, ] <- A_true %*% y[t - 1, ] + rnorm(2, sd = 0.1)
  fv <- fit_var(y, 1L)
  expect_equal(dim(fv$theta$A), c(2L, 2L))
  expect_equal(length(fv$theta$mu), 2L)
  expect_lt(max(abs(fv$theta$A - A_true)), 0.05)
  expect_error(fit_var(y, nrow(y) + 1L))
})

test_that("extract_mvar dimensions and naming", {
  n_ch <- ncol(train_sess$x[[1]])
  fit <- featEx4Train(train_sess$x, train_sess$y, "MVAR", params = list(order = 2L))
  expect_equal(ncol(fit$features), n_ch * n_ch * 2L + n_ch)   # A coefs + logQ
  noQ <- extract_mvar(train_sess$x, order = 2L, include_Q = FALSE)
  expect_equal(ncol(noQ), n_ch * n_ch * 2L)
  expect_true(all(grepl("^var_l", colnames(noQ))))
})

test_that("mcca aligns correlated views", {
  set.seed(4)
  n <- 120L
  s <- rnorm(n)                                   # shared latent factor
  V1 <- cbind(s + rnorm(n, sd = 0.1), matrix(rnorm(n * 3), n, 3))
  V2 <- cbind(s + rnorm(n, sd = 0.1), matrix(rnorm(n * 2), n, 2))
  m <- mcca_train(list(V1, V2), ncomp = 1L, reg = 1e-3)
  expect_s3_class(m, "mcca")

  tr <- mcca_transform(m, list(V1, V2))
  expect_equal(dim(tr$consensus), c(n, 1L))
  expect_equal(dim(tr$concat), c(n, 2L))
  expect_length(tr$scores, 2L)
  # The two views' leading canonical scores should be strongly correlated.
  cc <- abs(stats::cor(tr$scores[[1]][, 1], tr$scores[[2]][, 1]))
  expect_gt(cc, 0.8)

  expect_error(mcca_train(list(V1)))                            # needs >= 2 views
  expect_error(mcca_transform(m, list(V1)))                     # view count mismatch
})

# --- MSVAR engine / classifier / dispatcher feature -------------------------

make_var_classes <- function(seed, n_per = 16L, T = 180L) {
  A1 <- matrix(c(0.5, 0.1, 0, -0.2, 0.4, 0.1, 0, 0.1, 0.3), 3, 3)
  A2 <- matrix(c(-0.3, 0.2, 0.1, 0.1, -0.4, 0, 0, 0.1, 0.5), 3, 3)
  gen <- function(s, A) {
    set.seed(s)
    lapply(seq_len(n_per), function(i) {
      y <- matrix(0, T, 3)
      for (t in 2:T) y[t, ] <- A %*% y[t - 1, ] + rnorm(3, sd = 0.3)
      y
    })
  }
  list(x = c(gen(seed, A1), gen(seed + 100L, A2)),
       y = factor(rep(c("L", "R"), each = n_per)))
}

msvar_ctrl <- list(maxit = 6L, tol = 1e-4, init = "conditional", nseg = 15L)

test_that("fit_msvar returns a well-formed single-trial fit", {
  d <- make_var_classes(1L, n_per = 1L)
  f <- fit_msvar(d$x[[1]], M = 2L, p = 1L, control = msvar_ctrl)
  expect_equal(dim(f$theta$A), c(3L, 3L, 1L, 2L))
  expect_equal(dim(f$Ms), c(nrow(d$x[[1]]), 2L))
  expect_length(f$x, nrow(d$x[[1]]))
  expect_true(is.finite(f$LL))
})

test_that("msvar_train/predict separates distinct VAR dynamics", {
  tr <- make_var_classes(1L)
  te <- make_var_classes(50L)
  m <- msvar_train(tr$x, tr$y, M = 2L, p = 1L, control = msvar_ctrl, seed = 42L)
  expect_s3_class(m, "msvar")
  expect_setequal(m$classes, c("L", "R"))

  pr <- msvar_predict(m, te$x)
  expect_true(is.factor(pr$y_hat))
  expect_equal(dim(pr$loglik), c(length(te$x), 2L))
  expect_false(anyNA(pr$loglik))
  # Distinct dynamics -> should classify well above chance.
  expect_gt(mean(pr$y_hat == te$y), 0.75)
  # Predicted label is the higher-likelihood class.
  expect_equal(as.character(pr$y_hat),
               m$classes[max.col(pr$loglik, ties.method = "first")])
})

test_that("MSVAR dispatcher feature: shape, occupancy, determinism, serialization", {
  tr <- make_var_classes(1L)
  te <- make_var_classes(50L)

  fit <- featEx4Train(tr$x, tr$y, "MSVAR",
                      params = list(M = 2L, p = 1L, seed = 7L, control = msvar_ctrl))
  xtr <- fit$features
  xte <- featEx4Test(te$x, fit$object, "MSVAR")
  expect_equal(ncol(xtr), 2L)               # one log-likelihood per class
  expect_equal(ncol(xte), 2L)
  expect_false(anyNA(xte))
  expect_true(all(grepl("^msvar_LL_", colnames(xtr))))

  # Occupancy option appends nClass * M columns.
  fit_occ <- featEx4Train(tr$x, tr$y, "MSVAR",
                          params = list(M = 2L, p = 1L, seed = 7L,
                                        include_occupancy = TRUE, control = msvar_ctrl))
  expect_equal(ncol(fit_occ$features), 2L + 2L * 2L)

  # Same seed -> identical trained features.
  fit2 <- featEx4Train(tr$x, tr$y, "MSVAR",
                       params = list(M = 2L, p = 1L, seed = 7L, control = msvar_ctrl))
  expect_equal(fit$features, fit2$features, tolerance = 0)

  # Scoring is deterministic and survives a save/load round-trip.
  tmp <- tempfile(fileext = ".rds")
  saveRDS(fit$object, tmp)
  xte_loaded <- featEx4Test(te$x, readRDS(tmp), "MSVAR")
  expect_equal(xte, xte_loaded, tolerance = 0)
})

test_that("MSVAR feature requires at least two classes", {
  tr <- make_var_classes(1L)
  one_class <- factor(rep("L", length(tr$x)))
  expect_error(featEx4Train(tr$x, one_class, "MSVAR",
                            params = list(M = 2L, p = 1L, control = msvar_ctrl)))
})

# ===== Improvement-review regression & coverage tests (v0.3.0 polish) =========

test_that("FBCSP/FBCSSP work with a proper channel subset (regression)", {
  bands <- list(c(8, 12), c(12, 20))
  for (feat in c("FBCSP", "FBCSSP")) {
    fit <- featEx4Train(train_sess$x, train_sess$y, feat,
                        params = list(fs = 128L, channels = c(1L, 2L, 3L, 4L),
                                      frequency_bands = bands, ncomps = 2L))
    xte <- featEx4Test(test_sess$x, fit$object, feat)
    expect_true(is.matrix(fit$features))
    expect_equal(ncol(fit$features), ncol(xte))
    expect_false(anyNA(xte))
    expect_equal(fit$object$channels, c(1L, 2L, 3L, 4L))
  }
})

test_that("Riemannian jitter is batch-invariant", {
  obj <- featEx4Train(train_sess$x, train_sess$y, "Riemannian",
                      list(cov_type = "oas", use_filter = FALSE, jitter_sd = 1e-8))$object
  full  <- featEx4Test(test_sess$x, obj, "Riemannian")
  alone <- featEx4Test(test_sess$x[1], obj, "Riemannian")
  expect_equal(full[1, ], alone[1, ], tolerance = 0)
})

test_that("Riemannian/ACM_TS use_filter=TRUE path runs (riem_var_select + geodesic)", {
  for (feat in c("Riemannian", "ACM_TS")) {
    p <- if (feat == "Riemannian") list(cov_type = "oas", use_filter = TRUE) else
      list(order = 2L, delay = 1L, shrinkage = "oas", use_filter = TRUE)
    fit <- featEx4Train(train_sess$x, train_sess$y, feat, p)
    xte <- featEx4Test(test_sess$x, fit$object, feat)
    expect_true(is.matrix(fit$features))
    expect_equal(ncol(fit$features), ncol(xte))
    expect_false(anyNA(fit$features)); expect_false(anyNA(xte))
    expect_false(is.null(fit$object$varsel))
  }
})

test_that("ATM bin>1 (vectorized binning) round-trips", {
  fit <- featEx4Train(train_sess$x, train_sess$y, "ATM",
                      params = list(z = 2.0, bin = 2L, tuneBin = FALSE))
  xte <- featEx4Test(test_sess$x, fit$object, "ATM")
  expect_equal(ncol(fit$features), ncol(xte))
  expect_false(anyNA(xte))
})

test_that("MDM metric-aware distance is coherent for every metric", {
  for (mt in c("euclid", "logeuclid", "riemann")) {
    m <- mdm_train(train_sess$x, train_sess$y, cov_type = "oas", metric = mt)
    D <- predict(m, test_sess$x, type = "distance")
    expect_equal(dim(D), c(length(test_sess$x), nlevels(train_sess$y)))
    expect_false(anyNA(D))
    expect_true(all(D >= 0))
    # predict() class agrees with the legacy mdm_predict() factor
    expect_equal(as.character(predict(m, test_sess$x)),
                 as.character(mdm_predict(m, test_sess$x)))
  }
})

test_that("MDM cov_type='cov' guards against p>=n trials", {
  short <- lapply(1:6, function(i) matrix(rnorm(3 * 5), 3, 5))  # 3 samples, 5 ch
  y <- factor(rep(1:2, 3))
  expect_error(mdm_train(short, y, cov_type = "cov"))
  expect_s3_class(mdm_train(short, y, cov_type = "oas"), "mdm")  # shrinkage is fine
})

test_that("S3 predict/print methods for mdm/msvar/mcca", {
  m  <- mdm_train(train_sess$x, train_sess$y, metric = "logeuclid")
  ms <- msvar_train(train_sess$x, train_sess$y, M = 2L, p = 1L, seed = 1L,
                    control = list(maxit = 2L, nseg = 10L))
  expect_true(is.factor(predict(m, test_sess$x)))
  expect_true(is.matrix(predict(m, test_sess$x, type = "distance")))
  expect_true(is.factor(predict(ms, test_sess$x)))
  expect_true(is.matrix(predict(ms, test_sess$x, type = "loglik")))
  v1 <- featEx4Train(train_sess$x, train_sess$y, "logvar", list())$features
  v2 <- featEx4Train(train_sess$x, train_sess$y, "Hjorth", list())$features
  mc <- mcca_train(list(v1, v2), ncomp = 2L)
  expect_output(print(m), "MDM")
  expect_output(print(ms), "Markov-switching")
  expect_output(print(mc), "MCCA")
})

test_that("multiclass_EL trains, predicts, and exposes a working closure", {
  fit <- featEx4Train(train_sess$x, train_sess$y, "logvar", list())
  cl <- multiclass_EL(fit$features, train_sess$y, alpha = 0.5, lambda = 0.1, maxit = 80L)
  expect_equal(nrow(cl$coef), ncol(fit$features))
  expect_equal(ncol(cl$coef), nlevels(train_sess$y))
  pred <- cl$predict(featEx4Test(test_sess$x, fit$object, "logvar"))
  expect_true(is.factor(pred))
  expect_length(pred, length(test_sess$x))
  expect_setequal(levels(pred), levels(train_sess$y))
})

test_that("logvar_pca tolerates a constant/dead channel (no prcomp crash)", {
  xc <- lapply(train_sess$x, function(M) { M[, 3] <- 0; M })
  fit <- featEx4Train(xc, train_sess$y, "logvar_pca", list(ncomps = 4L))
  xte <- featEx4Test(test_sess$x, fit$object, "logvar_pca")
  expect_true(is.matrix(fit$features))
  expect_equal(ncol(fit$features), ncol(xte))
  expect_false(anyNA(fit$features))
})

test_that("mcca center/scale toggles run and snake aliases dispatch", {
  v1 <- featEx4Train(train_sess$x, train_sess$y, "logvar", list())$features
  v2 <- featEx4Train(train_sess$x, train_sess$y, "Hjorth", list())$features
  for (ce in c(TRUE, FALSE)) for (sc in c(TRUE, FALSE)) {
    mc <- mcca_train(list(v1, v2), ncomp = 2L, center = ce, scale = sc)
    expect_equal(dim(mcca_transform(mc, list(v1, v2))$consensus), c(nrow(v1), 2L))
  }
  # snake_case aliases delegate to the camelCase entry points
  a <- feat_ex_train(train_sess$x, train_sess$y, "logvar", list())
  expect_equal(a$features,
               featEx4Train(train_sess$x, train_sess$y, "logvar", list())$features)
})

test_that("freqBandSelect by_label=FALSE (fixed-window) path runs", {
  set.seed(2)
  fs <- 128L; L <- 512L; nblk <- 6L; nc <- 4L
  blocks <- lapply(seq_len(nblk), function(b) {
    cls <- if (b %% 2L == 1L) 1L else 2L
    t <- seq(0, (L - 1) / fs, length.out = L)
    M <- matrix(rnorm(L * nc, sd = 0.5), L, nc)
    M[, 1] <- M[, 1] + 2 * sin(2 * pi * (if (cls == 1L) 10 else 22) * t)
    list(x = M, y = rep(cls, L))
  })
  X <- do.call(rbind, lapply(blocks, `[[`, "x"))
  Y <- unlist(lapply(blocks, `[[`, "y"))
  fsel <- freqBandSelect(X, Y, fs = fs, flo = 4, fhi = 40, winlen = 256L,
                         by_label = FALSE)
  expect_equal(length(fsel$freqrange), 2L)
  expect_true(all(is.finite(fsel$freqrange)))
  expect_true(length(fsel$trials) >= 2L)
})
