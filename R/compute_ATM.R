# Module: Avalanche Transition Matrix (ATM) feature extraction.
# Scope: Converts z-scored binary avalanche activity into transition-probability
# features and provides train/test wrappers with consistent output dimensions.

#' Fit channel-wise z-score parameters across training trials.
#'
#' @param trialList List of trial matrices.
#' @return List with channel means (`mu`) and standard deviations (`sd`).
#' @export
atmZscoreFit <- function(trialList) {
  Xall <- do.call(rbind, lapply(trialList, function(m) as.matrix(m)))
  mu <- colMeans(Xall, na.rm = TRUE)
  sdv <- apply(Xall, 2, function(v) stats::sd(v, na.rm = TRUE))
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(mu = mu, sd = sdv)
}

#' Apply pre-fitted z-score normalization to one trial matrix.
#'
#' @param M Trial matrix.
#' @param zfit Output object from `atmZscoreFit()`.
#' @return Z-scored matrix.
#' @export
atmZscoreApply <- function(M, zfit) {
  sweep(sweep(as.matrix(M), 2, zfit$mu, "-"), 2, zfit$sd, "/")
}

#' Time-bin binary activity matrix by OR aggregation.
#'
#' @param B Binary activity matrix (`time x channels`).
#' @param bin Bin width in time steps.
#' @return Binned binary matrix.
#' @export
atmBinTimeBinary <- function(B, bin = 1L) {
  bin <- as.integer(bin)
  if (!is.finite(bin) || length(bin) != 1L || bin < 1L) {
    stop("`bin` must be an integer >= 1.")
  }
  if (bin == 1L) return(B)
  Tn <- nrow(B); D <- ncol(B)
  k  <- floor(Tn / bin)
  if (k < 1L) stop("Too few time points for the requested bin.")
  B2 <- B[seq_len(k * bin), , drop = FALSE]
  dim(B2) <- c(bin, k, D)
  Y <- apply(B2, c(2, 3), function(v) as.integer(any(v != 0L)))
  matrix(Y, nrow = k, ncol = D)
}

#' Detect contiguous avalanche segments.
#'
#' @param B Binary activity matrix.
#' @param minLen Minimum segment length to keep.
#' @return List of `[start, end]` index pairs.
#' @export
atmFindAvalanches <- function(B, minLen = 2L) {
  active <- rowSums(B) > 0L
  if (!any(active)) return(list())
  r <- rle(active)
  ends <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  idx <- which(r$values)
  segs <- Map(function(s, e) c(s, e), starts[idx], ends[idx])
  Filter(function(se) (se[2] - se[1] + 1L) >= minLen, segs)
}

#' Count channel transitions within avalanche segments.
#'
#' @param B Binary activity matrix.
#' @param segs Avalanche segments from `atmFindAvalanches()`.
#' @return List with pair counts `countIJ` and source counts `countI`.
#' @export
atmCounts <- function(B, segs) {
  D <- ncol(B)
  countIJ <- matrix(0L, D, D)
  countI  <- integer(D)
  for (se in segs) {
    s <- se[1]; e <- se[2]
    if (e <= s) next
    for (t in s:(e - 1L)) {
      iIdx <- which(B[t, ] != 0L)
      if (!length(iIdx)) next
      jIdx <- which(B[t + 1L, ] != 0L)
      if (length(jIdx)) countIJ[iIdx, jIdx] <- countIJ[iIdx, jIdx] + 1L
      countI[iIdx] <- countI[iIdx] + 1L
    }
  }
  list(countIJ = countIJ, countI = countI)
}

#' Convert transition counts to row-normalized transition probabilities.
#'
#' @param countIJ Transition count matrix.
#' @param countI Source activation counts.
#' @return Transition probability matrix.
#' @export
atmProbMatrix <- function(countIJ, countI) {
  D <- length(countI)
  P <- matrix(0, D, D)
  nz <- countI > 0L
  P[nz, ] <- sweep(countIJ[nz, , drop = FALSE], 1, countI[nz], "/")
  P
}

#' Estimate branching ratio from avalanche segments.
#'
#' @param B Binary activity matrix.
#' @param segs Avalanche segments.
#' @return Ratio of next-step to current-step active nodes, or NA if undefined.
#' @export
atmBranchingRatio <- function(B, segs) {
  num <- 0L; den <- 0L
  for (se in segs) {
    s <- se[1]; e <- se[2]
    if (e <= s) next
    for (t in s:(e - 1L)) {
      den <- den + sum(B[t, ] != 0L)
      num <- num + sum(B[t + 1L, ] != 0L)
    }
  }
  if (den == 0L) return(NA_real_)
  as.numeric(num) / as.numeric(den)
}

#' Choose the bin size whose branching ratio is closest to criticality.
#'
#' @param Z Z-scored matrix.
#' @param zThresh Threshold used to binarize activity.
#' @param useAbs Whether thresholding is absolute (`abs(Z)`).
#' @param minLen Minimum avalanche length.
#' @param binCandidates Candidate bin widths.
#' @return Selected bin width.
#' @export
atmChooseBin <- function(Z, zThresh = 3, useAbs = TRUE, minLen = 2L, binCandidates = c(1L, 2L, 3L)) {
  B0 <- if (useAbs) (abs(Z) >= zThresh) else (Z >= zThresh)
  binCandidates <- as.integer(binCandidates)
  binCandidates <- unique(binCandidates[is.finite(binCandidates) & binCandidates >= 1L])
  if (!length(binCandidates)) stop("`binCandidates` must contain at least one integer >= 1.")

  bestBin <- binCandidates[1]; bestGap <- Inf
  any_valid <- FALSE
  for (b in binCandidates) {
    Bb <- tryCatch(atmBinTimeBinary(B0, bin = b), error = function(e) NULL)
    if (is.null(Bb)) next
    segs <- atmFindAvalanches(Bb, minLen = minLen)
    if (!length(segs)) next
    sigma <- atmBranchingRatio(Bb, segs)
    if (!is.finite(sigma)) next
    any_valid <- TRUE
    # Target sigma close to 1.0 as a simple criticality heuristic.
    gap <- abs(sigma - 1.0)
    if (gap < bestGap) { bestGap <- gap; bestBin <- b }
  }
  if (!any_valid) {
    warning("No valid bin candidate produced finite branching ratio; using smallest candidate.")
    bestBin <- min(binCandidates)
  }
  bestBin
}

#' Compute ATM feature vector for one trial.
#'
#' @param trial Trial matrix.
#' @param zfit Z-score fit object.
#' @param zThresh Activity threshold.
#' @param bin Time-bin width.
#' @param minLen Minimum avalanche length.
#' @param useAbs Whether to threshold absolute z-scores.
#' @param includeDiag Whether diagonal transition terms are included.
#' @return Named numeric feature vector for one trial.
#' @export
atmTrialFeature <- function(trial, zfit, zThresh = 3, bin = 1L, minLen = 2L,
                            useAbs = TRUE, includeDiag = TRUE) {
  Z  <- atmZscoreApply(trial, zfit)
  B0 <- if (useAbs) (abs(Z) >= zThresh) else (Z >= zThresh)
  Bb <- atmBinTimeBinary(B0, bin = bin)
  segs <- atmFindAvalanches(Bb, minLen = minLen)
  
  if (!length(segs)) {
    # No avalanches: return a zero transition matrix with consistent dimensionality.
    D <- ncol(Bb); P <- matrix(0, D, D)
  } else {
    cnt <- atmCounts(Bb, segs)
    P   <- atmProbMatrix(cnt$countIJ, cnt$countI)
  }
  D <- ncol(P)
  if (includeDiag) {
    vec <- as.vector(t(P))
    cn  <- NULL
  } else {
    up <- upper.tri(P, diag = FALSE)
    vec <- P[up]
    cn  <- NULL
  }
  ch <- colnames(trial)
  if (!is.null(ch)) {
    safe <- function(s) gsub("[^A-Za-z0-9]+", "", s)
    ch <- vapply(ch, safe, FUN.VALUE = "")
    if (includeDiag) {
      cn <- paste0("atm", rep(ch, each = D), "to", rep(ch, times = D))
    } else {
      ij <- which(upper.tri(matrix(0, D, D), diag = FALSE), arr.ind = TRUE)
      cn <- paste0("atm", ch[ij[,1]], "to", ch[ij[,2]])
    }
  } else {
    if (is.null(cn)) {
      if (includeDiag) {
        cn <- paste0("atmC", rep(seq_len(D), each = D), "toC", rep(seq_len(D), times = D))
      } else {
        ij <- which(upper.tri(matrix(0, D, D), diag = FALSE), arr.ind = TRUE)
        cn <- paste0("atmC", ij[,1], "toC", ij[,2])
      }
    }
  }
  names(vec) <- cn
  vec
}

#' Train ATM extractor and produce training features.
#'
#' @param x List of trial matrices.
#' @param params Optional ATM parameter list.
#' @return List with `features` matrix and reusable `object` for testing.
#' @export
atmTrain <- function(x, params = NULL) {
  if (is.null(params)) params <- list()
  if (is.null(params$z))            params$z <- 3
  if (is.null(params$bin))          params$bin <- 1L
  if (is.null(params$tuneBin))      params$tuneBin <- FALSE
  if (is.null(params$minLen))       params$minLen <- 2L
  if (is.null(params$useAbs))       params$useAbs <- TRUE
  if (is.null(params$includeDiag))  params$includeDiag <- TRUE
  if (is.null(params$binCandidates)) params$binCandidates <- c(1L, 2L, 3L)
  params$bin <- as.integer(params$bin)
  if (!is.finite(params$bin) || params$bin < 1L) {
    stop("`params$bin` must be an integer >= 1.")
  }
  
  if (exists("perf_log_event", mode = "function")) perf_log_event("ATM_train_start")
  
  zfit <- atmZscoreFit(x)
  
  binSel <- params$bin
  if (isTRUE(params$tuneBin)) {
    # Tune bin on pooled normalized data to keep one consistent bin across trials.
    Xall <- do.call(rbind, lapply(x, function(M) atmZscoreApply(M, zfit)))
    binSel <- atmChooseBin(
      Z = Xall, zThresh = params$z, useAbs = params$useAbs,
      minLen = params$minLen, binCandidates = params$binCandidates
    )
  }
  
  D  <- ncol(x[[1]])
  p  <- if (isTRUE(params$includeDiag)) D * D else D * (D - 1) / 2
  # `vapply` enforces fixed output length and catches shape regressions early.
  Fe <- t(vapply(
    x,
    function(tr) atmTrialFeature(
      trial = tr, zfit = zfit, zThresh = params$z, bin = binSel,
      minLen = params$minLen, useAbs = params$useAbs, includeDiag = params$includeDiag
    ),
    FUN.VALUE = numeric(p)
  ))
  
  if (exists("perf_log_event", mode = "function")) perf_log_event("ATM_train_end")
  
  list(
    features = Fe,
    object = list(
      atmType       = "ATM",
      zfit          = zfit,
      z             = params$z,
      bin           = binSel,
      minLen        = params$minLen,
      useAbs        = params$useAbs,
      includeDiag   = params$includeDiag,
      binCandidates = params$binCandidates
    )
  )
}

#' Transform test trials to ATM features using a fitted ATM object.
#'
#' @param x Trial list or single trial matrix.
#' @param object Fitted ATM object from `atmTrain()`.
#' @return Feature matrix.
#' @export
atmTest <- function(x, object) {
  if (!is.list(x)) x <- list(x)
  D <- ncol(x[[1]])
  p <- if (isTRUE(object$includeDiag)) D * D else D * (D - 1) / 2
  
  if (exists("perf_log_event", mode = "function")) perf_log_event("ATM_test_start")
  
  Fe <- t(vapply(
    x,
    function(tr) atmTrialFeature(
      trial = tr, zfit = object$zfit, zThresh = object$z, bin = object$bin,
      minLen = object$minLen, useAbs = object$useAbs, includeDiag = object$includeDiag
    ),
    FUN.VALUE = numeric(p)
  ))
  
  if (exists("perf_log_event", mode = "function")) perf_log_event("ATM_test_end")
  
  Fe
}
