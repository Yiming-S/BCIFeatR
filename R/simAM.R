# Module: SimAM parameter-free attention re-weighting for EEG trials.
# Scope: Per-trial, per-sample attention weights derived from a local energy
# function, used as an optional conditioning step before feature extraction.

#' Re-weight EEG trials with the parameter-free SimAM attention module.
#'
#' Computes per-element attention weights from a local energy function and
#' multiplies the signal by `1 / (1 + energy)`. The transform carries no fitted
#' state, so train and test trials can be re-weighted independently.
#'
#' Reference: Yang, Zhang, Li & Xie (2021), "SimAM: A Simple, Parameter-Free
#' Attention Module for Convolutional Neural Networks", ICML.
#'
#' @param x A numeric trial matrix (`samples x channels`) or a list of such
#'   matrices.
#' @param method Energy variant: `"vanilla"` (mean/variance), `"robust"`
#'   (median/MAD), or `"corr"` (correlation-aware, a.k.a. Co-SimAM).
#' @param lambda Tuning parameter; if `NULL` a per-method default is used. For
#'   `"robust"` it is a scale-invariant relative bandwidth on the data's robust
#'   spread (default `1`; the spread auto-adapts to the data scale, so the
#'   energy is `O(1)` for ordinary samples as in the vanilla variant). For
#'   `"corr"` it is the cross-channel coupling strength (default `0.1`; `0`
#'   reduces to vanilla). Ignored by `"vanilla"`.
#' @param eps Numerical-stability constant.
#' @return The re-weighted input, matching the input structure (matrix or list).
#' @examples
#' M <- matrix(rnorm(256 * 4), 256, 4)
#' W <- simAM(M, method = "vanilla")
#' dim(W)
#' @export
simAM <- function(x,
                  method = c("vanilla", "robust", "corr"),
                  lambda = NULL,
                  eps = 1e-4) {
  method <- match.arg(method)
  if (is.null(lambda)) lambda <- if (method == "corr") 0.1 else 1
  is_single_mat <- is.matrix(x)
  if (is_single_mat) x <- list(x)
  stopifnot(is.list(x), all(vapply(x, is.matrix, logical(1))))

  out <- lapply(x, function(mat) {
    if (method == "vanilla") {
      mu <- mean(mat)
      v  <- mean((mat - mu)^2)
      energy <- (mat - mu)^2 / (4 * v + eps)

    } else if (method == "robust") {
      med <- stats::median(mat)
      # Consistent (sigma-scale) robust spread via mad() (constant 1.4826), so
      # the denominator auto-adapts to the data scale and the energy is O(1) for
      # ordinary samples -- median/MAD then resist outliers. `lambda` is a
      # scale-invariant relative bandwidth (default 1 ~ vanilla scaling).
      s_rob <- stats::mad(mat) * lambda
      energy <- (abs(mat - med))^2 / (4 * s_rob^2 + eps)

    } else {                                           # "corr" (Co-SimAM)
      mu <- mean(mat)
      v  <- mean((mat - mu)^2)
      e_base <- (mat - mu)^2 / (4 * v + eps)
      R <- suppressWarnings(stats::cor(mat))
      R[!is.finite(R)] <- 0
      diag(R) <- 0
      energy <- e_base + lambda * (e_base %*% R)       # lambda = coupling strength
    }

    weight <- 1 / (1 + energy)
    weighted <- mat * weight
    dimnames(weighted) <- dimnames(mat)
    weighted
  })

  if (is_single_mat) out[[1]] else out
}
