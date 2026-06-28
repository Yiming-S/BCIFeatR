# Module: Lightweight classifier utilities.
# Scope: Provides binary Fisher LDA projection and multiclass one-vs-rest
# elastic-net logistic models used in smoke experiments.

#' Compute Fisher LDA direction for binary classes.
#'
#' @param x Feature matrix with shape `p x n`.
#' @param label Binary label vector of length `n`.
#' @return List with projection `direction` and normalized sample `score`.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(6 * 40), 6, 40)        # features x samples
#' lab <- rep(c("a", "b"), each = 20)
#' fl <- fisher_lda(X, lab)
#' length(fl$direction)
#' @export
fisher_lda <- function(x, label) {
  n <- ncol(x); p <- nrow(x)
  stopifnot(length(label) == n)
  ulabel <- unique(label)
  if (length(ulabel) != 2L) {
    stop("fisher_lda currently supports binary labels only.")
  }
  label1 <- (label == ulabel[1]); label2 <- !label1
  x1bar <- rowMeans(x[, label1, drop = FALSE])
  x2bar <- rowMeans(x[, label2, drop = FALSE])
  x[, label1] <- x[, label1] - x1bar
  x[, label2] <- x[, label2] - x2bar

  # Solve in an SVD basis to avoid explicit inversion of within-class scatter.
  svdx <- svd(x, nv = 0)
  rmv  <- (svdx$d < 1e-14)
  if (any(rmv)) {
    svdx$d <- svdx$d[!rmv]
    svdx$u <- svdx$u[, !rmv, drop = FALSE]
  }
  phi <- svdx$u %*% (crossprod(svdx$u, x1bar - x2bar) / svdx$d)
  phi <- as.vector(phi); phi <- phi / sqrt(sum(phi^2))

  score <- as.vector(crossprod(x, phi))
  score_sd <- stats::sd(score)
  if (!is.finite(score_sd) || score_sd == 0) score_sd <- 1
  score <- score / score_sd
  list(direction = phi, score = score)
}

#' Train multiclass one-vs-rest elastic-net logistic models.
#'
#' @param X Training matrix (`n_samples x n_features`).
#' @param y Class labels.
#' @param alpha Elastic-net mixing: 1=l1, 0=l2.
#' @param lambda Overall regularization strength.
#' @param lr Gradient step size.
#' @param maxit Maximum iterations per binary model.
#' @param tol Early-stop tolerance on objective improvement.
#' @return Model list with coefficients, class levels, and prediction closure.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 5), 40, 5)
#' y <- factor(rep(c("a", "b"), each = 20))
#' cl <- multiclass_EL(X, y, alpha = 0.5, lambda = 0.1)
#' cl$predict(X)
#' @export
multiclass_EL <- function(X, y, alpha = 0.5, lambda = 0.1,
                          lr = 0.01, maxit = 1000, tol = 1e-6) {

  # Internal helper: proximal gradient for one binary elastic-net logistic model.
  train_binary_enet <- function(Xb, yb, alpha, lambda, lr, maxit, tol) {
    n <- nrow(Xb); p <- ncol(Xb)
    w <- rep(0, p); cost_prev <- Inf
    for (it in 1:maxit) {
      z   <- as.vector(Xb %*% w)
      p_i <- 1 / (1 + exp(-z))
      loss_logit <- -mean(yb * log(p_i + 1e-15) + (1 - yb) * log(1 - p_i + 1e-15))
      pen_l1 <- sum(abs(w)); pen_l2 <- sum(w^2)
      cost <- loss_logit + lambda * (alpha * pen_l1 + (1 - alpha) * pen_l2 / 2)
      if (abs(cost_prev - cost) < tol) break
      cost_prev <- cost

      grad <- (t(Xb) %*% (p_i - yb)) / n + lambda * (1 - alpha) * w
      w_tilde <- w - lr * as.vector(grad)

      # Soft-thresholding step implements the l1 proximal operator.
      thr <- lr * lambda * alpha
      w <- sign(w_tilde) * pmax(0, abs(w_tilde) - thr)
    }
    w
  }

  # Internal helper: one-vs-rest wrapper over all class levels.
  train_ovr <- function(Xm, ym, alpha, lambda, lr, maxit, tol) {
    classes <- levels(factor(ym)); K <- length(classes); p <- ncol(Xm)
    B <- matrix(0, nrow = p, ncol = K)
    for (k in seq_along(classes)) {
      y_bin <- ifelse(ym == classes[k], 1, 0)
      B[, k] <- train_binary_enet(Xm, y_bin, alpha, lambda, lr, maxit, tol)
    }
    list(coef = B, classes = classes)
  }

  # Internal helper: convert OvR probabilities to final class labels.
  predict_ovr <- function(model, Xnew) {
    Z <- Xnew %*% model$coef
    P <- 1 / (1 + exp(-Z))
    idx <- max.col(P, ties.method = "first")
    factor(model$classes[idx], levels = model$classes)
  }

  y <- factor(y)
  if (nrow(X) != length(y)) stop("Number of rows in X must match length(y).")

  model_ovr <- train_ovr(X, y, alpha, lambda, lr, maxit, tol)
  list(
    coef = model_ovr$coef,
    classes = model_ovr$classes,
    alpha = alpha,
    lambda = lambda,
    predict = function(Xnew) predict_ovr(model_ovr, Xnew)
  )
}

# Internal: resolve a per-trial covariance estimator for MDM. Extends the
# shrinkage estimators with the plain sample covariance ("cov").
.mdm_cov_fun <- function(cov_type) {
  switch(tolower(cov_type),
         oas = oas_covariance,
         lw  = LW_covariance,
         cov = function(m) stats::cov(as.matrix(m)),
         stop("`cov_type` must be 'oas', 'lw', or 'cov'."))
}

# Internal: symmetric matrix logarithm / inverse square root via eigen.
.sym_logm <- function(S, eps = 1e-12) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% (t(e$vectors) * log(pmax(e$values, eps)))
}
.sym_inv_sqrt <- function(S, eps = 1e-12) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% (t(e$vectors) / sqrt(pmax(e$values, eps)))
}

# Internal: trial-by-class distance matrix under the model's mean metric. The
# distance geometry matches the metric used to build the class means (Frobenius
# for euclid, log-Euclidean for logeuclid, affine-invariant for riemann), so the
# nearest-mean rule minimizes a single divergence. Per-class quantities are
# precomputed once (not per trial).
.mdm_distances <- function(model, x) {
  if (!is.list(x)) x <- list(x)
  cov_fun <- .mdm_cov_fun(model$cov_type)
  eps <- if (is.null(model$epsilon)) 1e-12 else model$epsilon
  metric <- if (is.null(model$metric)) "riemann" else model$metric
  means <- model$means; K <- length(means)
  prep <- switch(metric,
                 euclid    = means,
                 logeuclid = lapply(means, .sym_logm, eps = eps),
                 riemann   = lapply(means, .sym_inv_sqrt, eps = eps))
  N <- length(x)
  D <- matrix(NA_real_, N, K, dimnames = list(NULL, model$classes))
  for (i in seq_len(N)) {
    C <- cov_fun(as.matrix(x[[i]]))
    Clog <- if (metric == "logeuclid") .sym_logm(C, eps) else NULL
    for (k in seq_len(K)) {
      pk <- prep[[k]]
      D[i, k] <- switch(metric,
        euclid    = sqrt(sum((C - pk)^2)),
        logeuclid = sqrt(sum((Clog - pk)^2)),
        riemann   = {
          W  <- pk %*% C %*% pk
          ew <- eigen((W + t(W)) / 2, symmetric = TRUE, only.values = TRUE)$values
          sqrt(sum(log(pmax(ew, eps))^2))
        })
    }
  }
  D
}

#' Train a Riemannian Minimum Distance to Mean (MDM) classifier.
#'
#' Estimates one SPD class mean per label from per-trial covariance matrices.
#' Test trials are then assigned to the nearest class mean under the
#' affine-invariant Riemannian distance. The geometry reuses
#' [riemannian_mean()] and [riemannian_distance()].
#'
#' @param x List of trial matrices (`samples x channels`).
#' @param y Class labels aligned with `x`.
#' @param cov_type Per-trial covariance estimator: `"oas"`, `"lw"`, or `"cov"`.
#' @param metric Class-mean metric forwarded to [riemannian_mean()]:
#'   `"riemann"` (affine-invariant Fréchet mean, default), `"logeuclid"`, or
#'   `"euclid"`.
#' @param epsilon Numerical floor for SPD operations.
#' @return An object of class `"mdm"`: a list with per-class `means`, the
#'   `classes` levels, and the `metric`/`cov_type`/`epsilon` used.
#' @seealso [mdm_predict()]
#' @examples
#' set.seed(1)
#' x <- lapply(1:10, function(i) matrix(rnorm(200 * 3), 200, 3))
#' y <- factor(rep(c("a", "b"), each = 5))
#' m <- mdm_train(x, y, metric = "logeuclid")
#' predict(m, x)
#' @export
mdm_train <- function(x, y,
                      cov_type = c("oas", "lw", "cov"),
                      metric = c("riemann", "logeuclid", "euclid"),
                      epsilon = 1e-6) {
  cov_type <- match.arg(cov_type)
  metric   <- match.arg(metric)
  if (!is.list(x) || length(x) == 0L) stop("`x` must be a non-empty trial list.")
  if (length(y) != length(x)) stop("length(y) must equal the number of trials.")
  if (anyNA(y)) stop("`y` contains NA labels.")
  if (cov_type == "cov") {
    bad <- which(vapply(x, function(m) nrow(as.matrix(m)) <= ncol(as.matrix(m)),
                        logical(1)))
    if (length(bad)) {
      stop(sprintf(paste0("cov_type='cov' needs more samples than channels per ",
                          "trial, but trial(s) %s are rank-deficient. Use 'oas' ",
                          "or 'lw' shrinkage instead."),
                   paste(utils::head(bad, 5), collapse = ", ")))
    }
  }

  cov_fun  <- .mdm_cov_fun(cov_type)
  cov_list <- lapply(x, function(tr) cov_fun(as.matrix(tr)))
  p <- nrow(cov_list[[1]])

  y <- factor(y)
  classes <- levels(y)
  means <- lapply(classes, function(cl) {
    sub <- cov_list[which(y == cl)]
    # Force a 3D array so single-trial classes keep a valid stack dimension.
    arr <- array(unlist(sub), dim = c(p, p, length(sub)))
    riemannian_mean(arr, epsilon = epsilon, metric = metric)
  })
  names(means) <- classes

  structure(
    list(means = means, classes = classes,
         metric = metric, cov_type = cov_type, epsilon = epsilon),
    class = "mdm"
  )
}

#' Predict class labels with a trained MDM classifier.
#'
#' Assigns each trial to the nearest class mean under the model's mean metric
#' (Frobenius for `euclid`, log-Euclidean for `logeuclid`, affine-invariant for
#' `riemann`), so the rule minimizes a single divergence for every metric.
#'
#' @param model An `"mdm"` object from [mdm_train()].
#' @param x Trial list or a single trial matrix.
#' @return A factor of predicted classes (levels match the training classes),
#'   with the trial-by-class distance matrix attached as attribute `"distances"`.
#'   See also the [predict.mdm()] S3 method.
#' @seealso [mdm_train()], [predict.mdm()]
#' @examples
#' set.seed(1)
#' x <- lapply(1:10, function(i) matrix(rnorm(200 * 3), 200, 3))
#' y <- factor(rep(c("a", "b"), each = 5))
#' m <- mdm_train(x, y, metric = "logeuclid")
#' pred <- mdm_predict(m, x)        # factor; attr(pred, "distances") for scores
#' @export
mdm_predict <- function(model, x) {
  if (!inherits(model, "mdm")) {
    stop("`model` must be an 'mdm' object from mdm_train().")
  }
  D <- .mdm_distances(model, x)
  out <- factor(model$classes[max.col(-D, ties.method = "first")],
                levels = model$classes)
  attr(out, "distances") <- D
  out
}

#' Predict method for MDM classifiers.
#'
#' @param object An `"mdm"` object from [mdm_train()].
#' @param x Trial list or a single trial matrix.
#' @param type `"class"` for the predicted-class factor (default) or
#'   `"distance"` for the trial-by-class distance matrix.
#' @param ... Unused.
#' @return A factor of class labels, or the distance matrix when
#'   `type = "distance"`.
#' @seealso [mdm_train()], [mdm_predict()]
#' @export
predict.mdm <- function(object, x, type = c("class", "distance"), ...) {
  type <- match.arg(type)
  D <- .mdm_distances(object, x)
  if (type == "distance") return(D)
  factor(object$classes[max.col(-D, ties.method = "first")],
         levels = object$classes)
}

#' @export
print.mdm <- function(x, ...) {
  cat(sprintf("Riemannian MDM classifier (metric = %s, cov = %s)\n",
              x$metric, x$cov_type))
  cat(sprintf("  %d classes: %s\n", length(x$classes),
              paste(x$classes, collapse = ", ")))
  invisible(x)
}
