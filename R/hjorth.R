# Module: Hjorth-parameter feature extraction.
# Scope: Per-channel Hjorth activity, mobility, and complexity descriptors,
# a classic low-cost time-domain feature family for EEG.

#' Extract Hjorth parameters (activity, mobility, complexity) per channel.
#'
#' Activity is the signal variance, mobility the square-root ratio of the
#' first-difference variance to the signal variance, and complexity the ratio of
#' the first-difference mobility to the signal mobility. Variance ratios are
#' floored by `epsilon` so degenerate (flat) channels stay finite.
#'
#' @param data List of trial matrices (`samples x channels`).
#' @param log_activity When `TRUE`, return log activity instead of raw variance
#'   so the activity term shares the scale of the mobility/complexity terms.
#' @param epsilon Numerical floor guarding the variance ratios.
#' @return Numeric matrix (`n_trials x (3 * n_channels)`): the activity block,
#'   then the mobility block, then the complexity block.
#' @export
hjorth_parameters <- function(data, log_activity = FALSE, epsilon = 1e-12) {
  if (!is.list(data) || length(data) == 0L) {
    stop("`data` must be a non-empty trial list.")
  }
  M1 <- as.matrix(data[[1]])
  nch <- ncol(M1)
  ch_names <- colnames(M1)
  if (is.null(ch_names)) ch_names <- paste0("ch", seq_len(nch))

  feats <- t(vapply(data, function(tr) {
    M  <- as.matrix(tr)
    d1 <- diff(M)
    d2 <- diff(d1)
    var0 <- .col_vars(M)
    var1 <- .col_vars(d1)
    var2 <- .col_vars(d2)
    # Short trials (few rows) can yield NA differences; treat them as flat.
    var0[!is.finite(var0)] <- epsilon
    var1[!is.finite(var1)] <- epsilon
    var2[!is.finite(var2)] <- epsilon
    var0s <- pmax(var0, epsilon)
    var1s <- pmax(var1, epsilon)
    mobility   <- sqrt(var1s / var0s)
    complexity <- sqrt(pmax(var2, epsilon) / var1s) / mobility
    activity   <- if (isTRUE(log_activity)) log(pmax(var0, 1e-16)) else var0
    c(activity, mobility, complexity)
  }, numeric(3 * nch)))

  colnames(feats) <- c(paste0("hjorth_activity_", ch_names),
                       paste0("hjorth_mobility_", ch_names),
                       paste0("hjorth_complexity_", ch_names))
  feats
}
