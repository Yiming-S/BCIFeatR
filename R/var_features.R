# Module: Multivariate autoregressive (VAR) feature extraction.
# Scope: Least-squares VAR(p) fit per trial; the AR coefficients (and optional
# residual log-variances) form a per-trial feature vector. This is the
# in-scope, single-regime feature realization of the VAR/SVAR family (the full
# Markov-switching MSVAR is provided separately via msvar_train/msvar_predict and
# the "MSVAR" dispatcher feature). The VAR fit itself lives in `fit_var()`
# (see msvar_engine.R).

#' Extract per-trial VAR(p) features.
#'
#' Fits a VAR(p) to each trial and returns the vectorized AR coefficients,
#' optionally appending the per-channel log innovation variance.
#'
#' @param data List of trial matrices (`samples x channels`).
#' @param order AR order `p` (positive integer).
#' @param include_Q When `TRUE`, append the per-channel log of the innovation
#'   covariance diagonal.
#' @return Numeric matrix (`n_trials x n_features`) with
#'   `channels^2 * order` coefficient features (plus `channels` residual features
#'   when `include_Q = TRUE`).
#' @examples
#' x <- lapply(1:6, function(i) matrix(rnorm(256 * 3), 256, 3))
#' mv <- extract_mvar(x, order = 2L)
#' dim(mv)
#' @export
extract_mvar <- function(data, order = 3, include_Q = TRUE) {
  if (!is.list(data) || length(data) == 0L) {
    stop("`data` must be a non-empty trial list.")
  }
  order <- as.integer(order)
  M1 <- as.matrix(data[[1]])
  n  <- ncol(M1)
  ch <- colnames(M1)
  if (is.null(ch)) ch <- paste0("ch", seq_len(n))

  nA <- n * n * order
  p_total <- if (isTRUE(include_Q)) nA + n else nA

  feats <- t(vapply(data, function(tr) {
    M <- as.matrix(tr)
    if (nrow(M) <= order) {
      stop("Each trial must have more samples than the VAR order.")
    }
    fv <- fit_var(M, order)
    a  <- as.vector(fv$theta$A)
    if (isTRUE(include_Q)) c(a, log(pmax(diag(fv$theta$Q), 1e-16))) else a
  }, numeric(p_total)))

  # A is channels x (channels*order); as.vector() is column-major, so column c
  # encodes lag = ((c-1) %/% n)+1 and source channel = ((c-1) %% n)+1.
  cols   <- n * order
  dst    <- rep(seq_len(n), times = cols)
  colidx <- rep(seq_len(cols), each = n)
  lag    <- ((colidx - 1L) %/% n) + 1L
  src    <- ((colidx - 1L) %% n) + 1L
  anames <- paste0("var_l", lag, "_", ch[src], "to", ch[dst])
  colnames(feats) <- if (isTRUE(include_Q)) {
    c(anames, paste0("var_logQ_", ch))
  } else {
    anames
  }
  feats
}
