# Module: Covariance shrinkage estimators.
# Scope: Provides Ledoit-Wolf and OAS covariance estimators used by manifold
# and tangent-space feature branches.

#' Estimate covariance with Ledoit-Wolf shrinkage toward scaled identity.
#'
#' @param x Numeric matrix (`samples x variables`).
#' @return Symmetric positive semi-definite covariance estimate.
#' @export
LW_covariance <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x); p <- ncol(x)
  if (n < 2L || p < 1L) stop("LW_covariance requires at least 2 rows and 1 column.")
  x <- x - matrix(colMeans(x), n, p, TRUE)
  S <- crossprod(x) / n
  m <- mean(diag(S))
  # `d2` is target distance, `bbar2` is variance of sample covariance entries.
  d2    <- sum((S - diag(m, p))^2) / p
  bbar2 <- (sum(rowSums(x^2)^2) - n * sum(S^2)) / (p * n^2)
  if (!is.finite(d2) || d2 <= .Machine$double.eps) {
    out <- diag(m, p)
    return((out + t(out)) / 2)
  }
  b2 <- min(bbar2, d2)
  rho <- b2 / d2
  if (!is.finite(rho)) rho <- 1
  rho <- max(min(rho, 1), 0)
  S <- (1 - rho) * S
  diag(S) <- diag(S) + rho * m
  (S + t(S)) / 2
}

#' Estimate covariance with Oracle Approximating Shrinkage (OAS).
#'
#' @param X Numeric matrix (`samples x variables`).
#' @param eps Small threshold for safe denominator handling.
#' @return Symmetric positive semi-definite covariance estimate.
#' @export
oas_covariance <- function(X, eps = 1e-12) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (n < 2L || p < 1L) stop("oas_covariance requires at least 2 rows and 1 column.")
  Xc <- scale(X, center = TRUE, scale = FALSE)
  emp_cov <- crossprod(Xc) / n
  emp_cov <- (emp_cov + t(emp_cov)) / 2

  trS  <- sum(diag(emp_cov))
  trS2 <- sum(emp_cov * emp_cov)

  # Closed-form OAS shrinkage intensity.
  num <- (1 - 2 / p) * trS2 + trS^2
  den <- (n + 1 - 2 / p) * (trS2 - (trS^2) / p)
  rho <- if (!is.finite(den) || abs(den) < eps) 1 else num / den
  rho <- max(min(rho, 1), 0)

  mu <- trS / p
  F  <- diag(mu, p)
  Sigma <- (1 - rho) * emp_cov + rho * F
  (Sigma + t(Sigma)) / 2
}
