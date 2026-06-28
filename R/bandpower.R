# Module: Band-power feature extraction.
# Scope: Per-channel log band power (or within-channel relative band power)
# computed over a Butterworth filter bank, with one filter design reused
# across all trials on the hot path.

#' Extract per-channel band-power features over a frequency-band bank.
#'
#' Each band is band-pass filtered and the per-channel variance is taken as a
#' band-power estimate (Parseval's relation). The log band power is returned by
#' default, or the within-channel relative band power (fraction of total power
#' across the bank) when `relative = TRUE`.
#'
#' @param data List of trial matrices (`samples x channels`).
#' @param fs Sampling rate in Hz.
#' @param frequency_bands Non-empty list of `c(low, high)` bands in Hz.
#' @param order Butterworth filter order.
#' @param relative When `TRUE`, divide each channel's per-band powers by their
#'   sum across the bank (per trial) and return the fractions in `(0, 1)`;
#'   otherwise return log band power.
#' @param filters Optional pre-built filter bank from `build_bandpass_bank()`;
#'   when supplied, `fs` and `order` are ignored. Use this to avoid redesigning
#'   identical filters once per trial.
#' @return Numeric matrix (`n_trials x (n_bands * n_channels)`) in band-major
#'   order: all channels of band 1, then band 2, and so on.
#' @export
extract_bandpower <- function(data, fs, frequency_bands,
                              order = 3, relative = FALSE, filters = NULL) {
  if (!is.list(data) || length(data) == 0L) {
    stop("`data` must be a non-empty trial list.")
  }
  if (is.null(filters)) filters <- build_bandpass_bank(fs, frequency_bands, order)
  nbands <- length(filters)
  M1 <- as.matrix(data[[1]])
  nch <- ncol(M1)
  ch_names <- colnames(M1)
  if (is.null(ch_names)) ch_names <- paste0("ch", seq_len(nch))

  feats <- t(vapply(data, function(tr) {
    M <- as.matrix(tr)
    # Per-band power = variance of the band-pass filtered signal (nch x nbands).
    pw <- vapply(filters, function(f) .col_vars(gsignal::filtfilt(f, M)),
                 numeric(nch))
    if (isTRUE(relative)) {
      denom <- rowSums(pw)
      denom[!is.finite(denom) | denom <= 0] <- 1
      val <- pw / denom
    } else {
      val <- log(pmax(pw, 1e-16))
    }
    as.vector(val)
  }, numeric(nch * nbands)))

  colnames(feats) <- paste0("bp_b", rep(seq_len(nbands), each = nch),
                            "_", rep(ch_names, times = nbands))
  feats
}
