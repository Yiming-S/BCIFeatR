# Module: Frequency-domain helpers for filter-bank pipelines.
# Scope: Applies Butterworth filter banks and estimates discriminative
# pass bands from PSD-label correlations.

#' Apply a Butterworth filter bank to one trial matrix.
#'
#' @param eeg Trial matrix (`samples x channels`).
#' @param fs Sampling rate in Hz.
#' @param order Filter order.
#' @param bands List of pass bands, each as `c(low, high)` in Hz.
#' @return List of filtered trial matrices, one per band.
#' @export
freqBank <- function(eeg, fs, order = 3,
                     bands = list(c(8, 13), c(13, 18), c(18, 23), c(23, 28))) 
{	
  band_filters <- lapply(
    bands,
    function(band) gsignal::butter(order, band / (fs / 2), type = "pass")
  )
  filtered_eeg <- lapply(band_filters, gsignal::filtfilt, x = eeg)
  
  return(filtered_eeg)
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
  if (!exists("pwelch", mode = "function"))
    stop("Function 'pwelch' is required but not found. Please define or source it.")
  if (!exists("splitTimeRange", mode = "function"))
    stop("Function 'splitTimeRange' is required but not found. Please define or source it.")
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
  freq <- pwelch(x[trials[[1]],], fs = fs, 
                 window = winlen, overlap = overlap)$freq
  keep <- (freq >= flo & freq <= fhi)
  if (!any(keep)) {
    stop("No frequencies fall inside [flo, fhi].")
  }
  freq <- freq[keep]
  nfreq <- length(freq)
  psd <- array(dim = c(nfreq, nc, ntrials))
  
  for (i in 1:ntrials) {
    psd[,,i] <- pwelch(x[trials[[i]],], fs = fs, 
                       window = winlen, overlap = overlap)$spec[keep,]
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
