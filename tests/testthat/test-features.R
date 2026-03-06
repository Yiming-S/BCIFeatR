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
  ))
)

for (feat in names(feature_cases)) {
  test_that(sprintf("%s train/test roundtrip", feat), {
    run_one_feature(feat, feature_cases[[feat]]$params, train_sess, test_sess)
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
