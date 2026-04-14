
# Module: CSP/FBCSP/FBCSSP spatial filtering utilities.
# Scope: Learns class-discriminative spatial filters and extracts log-variance
# features after one-stage or two-stage filter-bank projections.

#' Train one-vs-one multiclass CSP filters.
#'
#' @param x List of trial matrices (samples x channels).
#' @param labels Class labels for each trial.
#' @param ncomps Number of components kept per class pair.
#' @param mean_type Class covariance mean type: `euclidean` or `riemannian`.
#' @param max_iterations Maximum iterations for Riemannian mean.
#' @param epsilon Numerical floor for Riemannian computations.
#' @param metric Mean metric used when `mean_type = "riemannian"`; forwarded to
#'   `riemannian_mean()`. Defaults to `"euclid"` for speed.
#' @return List with flattened `filter`, projected `pattern`, eigen `var_ratio`,
#'   and class `pairs`.
#' @export
multiclass_csp <- function(x, labels,
                           ncomps = 6,
                           mean_type = "euclidean",
                           max_iterations = 100,
                           epsilon = 1e-6,
                           metric = c("euclid", "logeuclid", "riemann")) {
  metric <- match.arg(metric)
  stopifnot(is.list(x))
  stopifnot(length(x) == length(labels))
  labels <- as.integer(factor(labels))
  nclass <- length(unique(labels))
  ntrials <- length(labels)
  nchans <- sapply(x, ncol)
  if (any(nchans != nchans[1])) {
    stop("All matrices in 'x' must have the same number of columns (channels).")
  }
  nchans <- nchans[1]
  ncomps <- min(ncomps, nchans)
  x <- lapply(x, scale, center = TRUE, scale = FALSE)
  calc_mean <- function(data_idx) {
    if (mean_type == "riemannian") {
      cov_matrices <- lapply(x[data_idx], crossprod)
      riemannian_mean(
        array(unlist(cov_matrices),
              dim = c(nchans, nchans, length(cov_matrices))),
        max_iterations = max_iterations,
        epsilon = epsilon,
        metric = metric
      )
    } else if (mean_type == "euclidean") {
      sum_matrix <- matrix(0, nchans, nchans)
      for (i in data_idx) {
        sum_matrix <- sum_matrix + crossprod(x[[i]])
      }
      sum_matrix / length(data_idx)
    } else {
      stop("Invalid mean_type. Use 'riemannian' or 'euclidean'.")
    }
  }
  pairs <- combn(nclass, 2)
  npairs <- ncol(pairs)
  spatial_filters <- array(dim = c(nchans, ncomps, npairs))
  var_ratio <- matrix(NA, ncomps, npairs)
  # Solve each OvO generalized eigenproblem via Cholesky whitening of cov2:
  # cov1 v = lambda cov2 v   ⇔   (U')^{-1} cov1 U^{-1} w = lambda w, v = U^{-1} w
  # where U is upper-tri with U'U = cov2. Avoids geigen() which is notably
  # slower on moderate p and drags in a heavy dependency on the CSP hot path.
  for (k in 1:npairs) {
    i1 <- pairs[1, k]
    i2 <- pairs[2, k]
    idx1 <- which(labels == i1)
    idx2 <- which(labels == i2)
    cov1 <- calc_mean(idx1)
    cov2 <- calc_mean(idx2)
    C1 <- (cov1 + t(cov1)) / 2
    C2 <- (cov2 + t(cov2)) / 2
    U <- tryCatch(chol(C2), error = function(e) {
      diag(C2) <- diag(C2) + epsilon * max(1, mean(diag(C2)))
      chol(C2)
    })
    invU <- backsolve(U, diag(nrow(U)))
    M <- crossprod(invU, C1) %*% invU
    M <- (M + t(M)) / 2
    eig <- eigen(M, symmetric = TRUE)
    vals <- eig$values
    # Rank by deviation from 1 (discriminative contrast), matching original logic.
    abs_vals <- pmax(abs(vals), 1e-12)
    ordv <- pmax(abs_vals, 1 / abs_vals)
    idx <- order(ordv, decreasing = TRUE)[seq_len(ncomps)]
    spatial_filters[, , k] <- invU %*% eig$vectors[, idx, drop = FALSE]
    var_ratio[, k] <- vals[idx]
  }
  # Concatenate pair-specific filters into a single projection matrix.
  spatial_filters <- aperm(spatial_filters, c(1, 3, 2))
  dim(spatial_filters) <- c(nchans, npairs * ncomps)
  filtered_data <- lapply(x, "%*%", y = spatial_filters)
  list(
    filter    = spatial_filters,
    pattern   = filtered_data,
    var_ratio = var_ratio,
    pairs     = pairs
  )
}

#' Apply per-band CSP filters to each trial and concatenate outputs.
#'
#' @param data List of trial matrices.
#' @param filters List of projection matrices, typically one per band.
#' @return List of projected trial matrices.
#' @export
apply_csp_filters <- function(data, filters) {
  lapply(data, function(trial) {
    filtered <- lapply(filters, function(f) as.matrix(trial) %*% f)
    do.call(cbind, filtered)
  })
}

#' Extract log-variance features from projected trial matrices.
#'
#' @param filtered_data List of projected trial matrices.
#' @return Numeric feature matrix (`n_trials` x `n_features`).
#' @export
extract_logvar_features <- function(filtered_data) {
  features <- lapply(filtered_data, function(trial_matrix) {
    if (is.complex(trial_matrix)) trial_matrix <- Mod(trial_matrix)
    log(pmax(.col_vars(trial_matrix), 1e-16))
  })
  do.call(rbind, features)
}

#' Train and package CSP-family filters for train/test pipelines.
#'
#' @param data_train List of training trials.
#' @param labels_train Training labels.
#' @param fs Sampling rate.
#' @param channels Optional subset of channels.
#' @param csp_type One of `FBCSP` or `FBCSSP`.
#' @param frequency_bands List of band ranges.
#' @param ncomps Components per stage.
#' @return List with projected training patterns and trained filters.
#' @export
prepare_csp_filters <- function(data_train, labels_train,
                                fs, channels, csp_type,
                                frequency_bands = list(c(8, 30), c(11, 20),
                                                       c(21, 30), c(31, 40)),
                                ncomps = c(3, 3)) {
  # Stage-1 FBCSP is needed in both paths; FBCSSP consumes its output directly
  # (its stage-detection branch reuses the stage-1 result without recomputing).
  stage1 <- FBCSP(data_list = data_train,
                  labels = labels_train,
                  fs = fs,
                  frequency_bands = frequency_bands,
                  channels = channels,
                  ncomps = ncomps)
  if (csp_type == "FBCSP") {
    list(pattern = stage1$pattern, filter = stage1$filter)
  } else if (csp_type == "FBCSSP") {
    stage2 <- FBCSSP(stage1,
                     labels_train,
                     fs = fs,
                     frequency_bands = frequency_bands,
                     channels = channels,
                     ncomps1 = ncomps,
                     ncomps2 = ncomps)
    list(pattern = stage2$pattern,
         filter  = stage2$filter,
         filter2 = stage2$filter2)
  } else {
    stop("Invalid csp_type. Choose either 'FBCSP' or 'FBCSSP'.")
  }
}

#' Train Filter Bank CSP and return per-band features.
#'
#' @param data_list List of trial matrices.
#' @param labels Trial labels.
#' @param fs Sampling rate.
#' @param frequency_bands List of pass bands.
#' @param channels Optional channel subset.
#' @param ncomps Components per band.
#' @return List with per-band filters, projected patterns, features, and eigen ratios.
#' @export
FBCSP <- function(data_list, labels, fs,
                  frequency_bands,
                  channels = NULL, ncomps = 5) {
  stopifnot(is.list(data_list))
  ntrials <- length(data_list)
  nbands <- length(frequency_bands)
  nchannels <- ncol(data_list[[1]])
  if (!is.null(channels) && !setequal(channels, 1:nchannels)) {
    data_list <- lapply(data_list, function(x) x[, channels, drop = FALSE])
  }
  # Design Butterworth band filters once and reuse them across all trials; the
  # coefficients depend only on (fs, order, band), so the previous per-trial
  # design inside `freqBank` was pure overhead on the FBCSP hot path.
  band_filters <- lapply(
    frequency_bands,
    function(band) gsignal::butter(3, band / (fs / 2), type = "pass")
  )
  filtered_data_list <- lapply(data_list, function(trial) {
    lapply(band_filters, gsignal::filtfilt, x = trial)
  })
  filtered_data_list <- unlist(filtered_data_list, recursive = FALSE)
  dim(filtered_data_list) <- c(nbands, ntrials)
  spatial_filter <- feature <- var_ratio <- vector("list", nbands)
  pattern <- vector("list", ntrials * nbands)
  dim(pattern) <- c(ntrials, nbands)

  for (b in 1:nbands) {
    # Train CSP independently in each frequency band.
    CSP <- multiclass_csp(filtered_data_list[b, ], labels, ncomps)
    spatial_filter[[b]] <- CSP$filter
    pattern[, b] <- CSP$pattern
    var_csp <- sapply(CSP$pattern, function(x) {
      if (is.complex(x)) x <- Mod(x)
      storage.mode(x) <- "double"
      pmax(colMeans(x^2), 1e-8)
    })
    feature[[b]] <- log(t(var_csp))
    var_ratio[[b]] <- CSP$var_ratio
  }

  feature <- do.call(cbind, feature)

  spatial_pattern <- vector("list", ntrials)
  for (trial in 1:ntrials) {
    spatial_pattern[[trial]] <- do.call(cbind, pattern[trial, ])
  }

  return(list(
    filter  = spatial_filter,
    pattern = spatial_pattern,
    feature = feature,
    var_ratio = var_ratio
  ))
}

#' Train two-stage Filter Bank Common Sparse Spatial Pattern filters.
#'
#' @param data_list Raw trial list or a precomputed first-stage FBCSP object.
#' @param labels Trial labels.
#' @param fs Sampling rate when raw trials are provided.
#' @param frequency_bands Frequency bands used in first stage.
#' @param channels Optional channel subset.
#' @param ncomps1 Components for stage-1 (per band).
#' @param ncomps2 Components for stage-2 (across concatenated stage-1 outputs).
#' @return List with combined filters, second-stage filter, projected patterns, and features.
#' @export
FBCSSP <- function(data_list, labels, fs = NULL, frequency_bands,
                   channels = NULL, ncomps1 = 2, ncomps2 = 2)
{
  names1 <- c("filter", "pattern", "feature", "var_ratio")
  stage1 <- if (identical(names(data_list), names1)) {
    data_list
  } else {
    FBCSP(data_list, labels, fs, frequency_bands, channels, ncomps1)
  }
  stage2 <- multiclass_csp(stage1$pattern, labels, ncomps2)
  ncomps1 <- ncol(stage1$filter[[1]])
  ncomps2 <- ncol(stage2$filter)
  nbands <- length(stage1$filter)
  combined_filter <- vector("list", nbands)

  for (b in 1:nbands) {
    # Compose first-stage band filter with its corresponding stage-2 block.
    idx <- seq((b - 1) * ncomps1 + 1, b * ncomps1)
    combined_filter[[b]] <- stage1$filter[[b]] %*% stage2$filter[idx, , drop = FALSE]
  }

  combined_pattern <- stage2$pattern
  var_csp <- sapply(
    combined_pattern,
    function(mat) {
      if (is.complex(mat)) mat <- Mod(mat)
      storage.mode(mat) <- "double"
      pmax(colMeans(mat^2), 1e-8)
    }
  )

  combined_feature <- log(t(var_csp))

  return(list(
    filter    = combined_filter,
    filter2   = stage2$filter,
    pattern   = combined_pattern,
    feature   = combined_feature,
    var_ratio = stage2$var_ratio
  ))
}
