# Module: Lightweight classifier utilities.
# Scope: Provides binary Fisher LDA projection and multiclass one-vs-rest
# elastic-net logistic models used in smoke experiments.

#' Compute Fisher LDA direction for binary classes.
#'
#' @param x Feature matrix with shape `p x n`.
#' @param label Binary label vector of length `n`.
#' @return List with projection `direction` and normalized sample `score`.
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
