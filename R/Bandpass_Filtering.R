# Module: Frequency-domain helpers for filter-bank pipelines.
# Scope: Applies Butterworth filter banks and estimates discriminative
# pass bands from PSD-label correlations.

#' Build a Butterworth filter bank once for reuse across trials.
#'
#' @param fs Sampling rate in Hz.
#' @param bands List of pass bands, each as `c(low, high)` in Hz.
#' @param order Filter order.
#' @return List of filter design objects from `gsignal::butter`.
#' @export
build_bandpass_bank <- function(fs, bands, order = 3) {
  lapply(bands,
         function(band) gsignal::butter(order, band / (fs / 2), type = "pass"))
}

#' Apply a Butterworth filter bank to one trial matrix.
#'
#' @param eeg Trial matrix (`samples x channels`).
#' @param fs Sampling rate in Hz.
#' @param order Filter order.
#' @param bands List of pass bands, each as `c(low, high)` in Hz.
#' @param filters Optional pre-built filter bank from `build_bandpass_bank()`;
#'   when supplied, `fs`, `order`, and `bands` are ignored. Use this to avoid
#'   redesigning identical filters once per trial in filter-bank pipelines.
#' @return List of filtered trial matrices, one per band.
#' @export
freqBank <- function(eeg, fs, order = 3,
                     bands = list(c(8, 13), c(13, 18), c(18, 23), c(23, 28)),
                     filters = NULL) {
  if (is.null(filters)) filters <- build_bandpass_bank(fs, bands, order)
  lapply(filters, gsignal::filtfilt, x = eeg)
}

#' Split a labelled time series into per-trial index segments.
#'
#' Segments a sample-indexed label/marker vector into trials, either by maximal
#' runs of constant label (`by_label = TRUE`) or into fixed-length windows
#' within each constant-label run (`by_label = FALSE`).
#'
#' @param labels Integer-coercible vector of per-sample labels/markers.
#' @param by_label When `TRUE`, each maximal run of identical labels becomes one
#'   trial. When `FALSE`, each run is further chopped into windows of `seglen`
#'   samples, keeping only windows of at least `minlen` samples.
#' @param minlen Minimum window length for `by_label = FALSE`; defaults to
#'   `seglen`.
#' @param seglen Window length in samples; required when `by_label = FALSE`.
#' @return A list of integer index vectors, one per trial.
#' @examples
#' labels <- c(rep(1, 5), rep(2, 3), rep(1, 4))
#' splitTimeRange(labels, by_label = TRUE)
#' @export
splitTimeRange <- function(labels, by_label = TRUE,
                           minlen = NULL, seglen = NULL) {
  labels <- as.integer(labels)
  n <- length(labels)
  if (n < 1L) stop("`labels` must be non-empty.")

  if (by_label) {
    # Segment on label-value changes: each maximal constant run is one trial.
    breakpoints <- c(0L, which(diff(labels) != 0L), n)
    return(lapply(seq_len(length(breakpoints) - 1L), function(i) {
      (breakpoints[i] + 1L):breakpoints[i + 1L]
    }))
  }

  # Fixed-length windowing within each constant-label run.
  if (is.null(seglen)) stop("`seglen` is required when by_label = FALSE.")
  seglen <- as.integer(seglen)
  if (!is.finite(seglen) || seglen < 1L) stop("`seglen` must be an integer >= 1.")
  minlen <- if (is.null(minlen)) seglen else min(as.integer(minlen), seglen)

  idx   <- which(diff(labels) != 0L)
  start <- c(1L, idx + 1L)
  end   <- c(idx, n)
  trials <- list()
  for (s in seq_along(start)) {
    starts <- seq(start[s], end[s], by = seglen)
    ends   <- c(starts[-1] - 1L, end[s])
    valid  <- (ends - starts + 1L) >= minlen
    for (j in which(valid)) {
      trials <- c(trials, list(starts[j]:ends[j]))
    }
  }
  trials
}

#' Select an informative frequency range from PSD-label association.
#'
#' @param x Signal matrix (`samples x channels`).
#' @param y Label/time marker vector aligned with rows of `x`.
#' @param fs Sampling rate in Hz.
#' @param flo Lower frequency bound.
#' @param fhi Upper frequency bound.
#' @param winlen Welch window length; inferred when NULL.
#' @param overlap Welch overlap ratio in [0, 1).
#' @param by_label Whether to segment trials by label runs.
#' @param trials Optional precomputed trial index list.
#' @return List with selected `freqrange`, full `freqgrid`, and `trials`.
#' @export
freqBandSelect <- function(x, y, fs, 
                           flo = 1, fhi = fs/2, 
                           winlen = NULL, overlap = 0.5, 
                           by_label = TRUE, trials = NULL) {
  stopifnot(NROW(x) == length(y), is.matrix(x),
            is.numeric(y), is.numeric(fs),
            flo > 0, flo < fhi, fhi <= (fs / 2),
            overlap >= 0 && overlap < 1)
  nc <- ncol(x)
  ns <- nrow(x)
  if (is.null(winlen)) {
    winlen <- ceiling(2 / flo * fs)
  } else {
    stopifnot(winlen <= ns)
  }
  if (by_label) {
    if (is.null(trials)) {
      trials <- splitTimeRange(y, by_label = TRUE)
    }
    ntrials <- length(trials)
    y <- sapply(trials, function(idx) y[idx[1]])
  } else {
    if (is.null(trials)) {
      trials <- splitTimeRange(y, by_label = FALSE, seglen = winlen)
    }
    ntrials <- length(trials)
    y <- sapply(trials, function(idx) y[idx[1]])
  }
  if (!length(trials)) {
    stop("No trials were generated for frequency-band selection.")
  }
  if (length(unique(y)) < 2L) {
    stop("freqBandSelect requires at least two distinct trial labels.")
  }
  freq <- pwelch(x[trials[[1]], , drop = FALSE], fs = fs,
                 window = winlen, overlap = overlap)$freq
  keep <- (freq >= flo & freq <= fhi)
  if (!any(keep)) {
    stop("No frequencies fall inside [flo, fhi].")
  }
  freq <- freq[keep]
  nfreq <- length(freq)
  psd <- array(dim = c(nfreq, nc, ntrials))

  for (i in 1:ntrials) {
    # drop = FALSE / as.matrix keep the single-channel spectrum 2-D.
    spec <- as.matrix(pwelch(x[trials[[i]], , drop = FALSE], fs = fs,
                             window = winlen, overlap = overlap)$spec)
    psd[, , i] <- spec[keep, , drop = FALSE]
  }
  # Score each (frequency, channel) bin by correlation with trial labels.
  dB <- log10(psd)
  score <- apply(dB, 1:2, function(v) suppressWarnings(stats::cor(v, y, use = "complete.obs")))
  score[!is.finite(score)] <- 0
  top_k <- min(3L, nc)
  top_channels <- order(colSums(score^2), decreasing = TRUE)[seq_len(top_k)]
  score <- score[,top_channels, drop = FALSE]
  fmax <- which.max(abs(rowSums(score)))
  scorestar <- score
  flip <- score[fmax,] < 0
  flip[is.na(flip)] <- FALSE
  scorestar[,flip] <- -score[,flip]
  fscore <- rowSums(scorestar)
  fscore[!is.finite(fscore)] <- 0
  fstarmax <- which.max(fscore)
  f0 <- f1 <- fstarmax
  peak <- max(fscore[fstarmax], .Machine$double.eps)
  
  # Expand around peak while score remains above 5% of the peak response.
  while (f0 > 1 && isTRUE(fscore[f0-1] >= 0.05 * peak)) f0 <- f0 - 1L
  while (f1 < nfreq && isTRUE(fscore[f1+1] >= 0.05 * peak)) f1 <- f1 + 1L
  
  list(freqrange = c(freq[f0], freq[f1]), freqgrid = freq, trials = trials)
}
