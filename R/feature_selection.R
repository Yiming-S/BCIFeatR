# Module: Feature-selection utilities.
# Scope: Implements MI-based selection, PCA-based selection, and PCA+Fisher
# ranking over session-wise feature matrices.

#' Select informative features using mutual information (MIBIF-style).
#'
#' @param data_list List of feature matrices.
#' @param labels Session-level labels aligned with `data_list`.
#' @param k Fraction of features to keep.
#' @param method `combine` (use all sessions) or `random` (subset sessions).
#' @param num_random Number of sessions sampled when `method = 'random'`.
#' @return List with selected feature frame and original column indices.
#' @export
mibif <- function(data_list, labels, k = 0.3, method = 'combine', num_random = NULL) {
  if (!requireNamespace("infotheo", quietly = TRUE))
    stop("Package 'infotheo' is required for mibif(). Please install it.")
  method <- match.arg(method, c("combine", "random"))
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("`k` must be a positive numeric scalar.")
  }

  if (length(data_list) != length(labels)) {
    stop("The length of data_list and labels must be the same.")
  }

  if (method == 'random') {
    if (is.null(num_random)) {
      stop("num_random must be provided for the 'random' method.")
    }
    if (!is.numeric(num_random) || length(num_random) != 1L || !is.finite(num_random) ||
        num_random < 1 || as.integer(num_random) != num_random) {
      stop("num_random must be a positive integer.")
    }
    num_random <- as.integer(num_random)
    if (num_random > length(data_list)) {
      stop("num_random cannot be greater than the length of data_list.")
    }

    selected_indices <- sample(1:length(data_list), num_random)
    data_list <- data_list[selected_indices]
    labels <- labels[selected_indices]
  }
  class_vector <- unlist(mapply(rep, labels, sapply(data_list, nrow)))
  combined_data <- do.call(rbind, data_list)

  k <- ceiling(ncol(combined_data) * k)
  k <- min(k, ncol(combined_data))
  compute_mi <- function(feature, class) {
    return(infotheo::mutinformation(infotheo::discretize(feature), infotheo::discretize(class)))
  }
  feature_names <- colnames(combined_data)
  if (is.null(feature_names)) {
    feature_names <- paste0("V", 1:ncol(combined_data))
  }
  all_indices <- seq_len(ncol(combined_data))
  S <- list()
  selected_set <- logical(length(all_indices))
  selected_indices <- integer(0)
  while (length(S) < k) {
    # Greedy forward selection by maximal mutual information.
    mi_values <- vapply(all_indices, function(i) {
      if (selected_set[i]) -Inf
      else compute_mi(combined_data[, i], class_vector)
    }, numeric(1))
    best_idx <- which.max(mi_values)
    best_name <- feature_names[best_idx]

    S[[best_name]] <- combined_data[, best_idx]
    selected_set[best_idx] <- TRUE
    selected_indices <- c(selected_indices, best_idx)
  }
  return(list(selected_features = as.data.frame(S),
              selected_indices = selected_indices))
}

#' Select principal components by explained-variance rank.
#'
#' @param data_list List of feature matrices.
#' @param k Fraction of components to retain.
#' @param method `combine` or `random`.
#' @param num_random Number of sessions sampled when `method = 'random'`.
#' @return List containing reduced features, selected indices, and loadings.
#' @export
pca <- function(data_list, k = 0.9, method = 'combine', num_random = NULL) {
  method <- match.arg(method, c("combine", "random"))
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0 || k > 1) {
    stop("`k` must be a numeric scalar in (0, 1].")
  }
  if (method == 'random') {
    if (is.null(num_random)) {
      stop("num_random must be provided for the 'random' method.")
    }
    if (!is.numeric(num_random) || length(num_random) != 1L || !is.finite(num_random) ||
        num_random < 1 || as.integer(num_random) != num_random) {
      stop("num_random must be a positive integer.")
    }
    num_random <- as.integer(num_random)
    if (num_random > length(data_list)) {
      stop("num_random cannot be greater than the length of data_list.")
    }

    selected_indices <- sample(1:length(data_list), num_random)
    data_list <- data_list[selected_indices]
  }
  combined_data <- do.call(rbind, data_list)
  pca_result <- prcomp(combined_data, center = TRUE, scale. = TRUE)
  num_components <- min(ncol(pca_result$x), max(1L, ceiling(ncol(pca_result$x) * k)))
  variances <- pca_result$sdev^2
  sorted_indices <- order(variances, decreasing = TRUE)
  keep <- sorted_indices[1:num_components]
  selected_features <- pca_result$x[, keep, drop = FALSE]
  selected_indices <- sorted_indices[1:num_components]
  principal_components <- pca_result$rotation[, keep, drop = FALSE]
  channel_contribution <- rowSums(abs(principal_components))
  return(list(selected_features = as.data.frame(selected_features),
              selected_indices = selected_indices,
              principal_components = principal_components,
              channel_contribution = channel_contribution))
}

#' Rank PCA components with Fisher score and keep top subset.
#'
#' @param data_list List of feature matrices.
#' @param labels Session-level or sample-level labels.
#' @param k Fraction of principal components to retain.
#' @param method `combine` or `random`.
#' @param num_random Number of sessions sampled when `method = 'random'`.
#' @param center Whether to center features before PCA.
#' @param scale Whether to scale features before PCA.
#' @return List with Fisher ranking, PCA model, selected indices, and variance stats.
#' @export
fisher <- function(data_list, labels, k = 0.9, method = 'combine',
                   num_random = NULL, center = TRUE, scale = TRUE) {
  method <- match.arg(method, c("combine", "random"))
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0 || k > 1) {
    stop("`k` must be a numeric scalar in (0, 1].")
  }

  if (!is.list(data_list) || length(data_list) == 0L) {
    stop("`data_list` must be a non-empty list of matrices.")
  }
  if (!all(vapply(data_list, is.matrix, logical(1)))) {
    stop("All elements of `data_list` must be matrices.")
  }

  n_sessions_raw <- length(data_list)
  rows_per_session <- vapply(data_list, nrow, integer(1))
  total_rows_raw <- sum(rows_per_session)
  if (method == 'random') {
    if (is.null(num_random)) {
      stop("num_random must be provided for the 'random' method.")
    }
    if (!is.numeric(num_random) || length(num_random) != 1L || !is.finite(num_random) ||
        num_random < 1 || as.integer(num_random) != num_random) {
      stop("num_random must be a positive integer.")
    }
    num_random <- as.integer(num_random)
    if (num_random > length(data_list)) {
      stop("num_random cannot be greater than the length of data_list.")
    }

    selected_indices <- sample(1:length(data_list), num_random)
    if (length(labels) == n_sessions_raw) {
      labels <- labels[selected_indices]
    } else if (length(labels) == total_rows_raw) {
      ends <- cumsum(rows_per_session)
      starts <- c(1L, head(ends, -1L) + 1L)
      keep_rows <- unlist(Map(
        function(s, e) seq.int(s, e),
        starts[selected_indices],
        ends[selected_indices]
      ), use.names = FALSE)
      labels <- labels[keep_rows]
    }
    data_list <- data_list[selected_indices]
  }
  data <- do.call(rbind, data_list)

  total_rows <- nrow(data)
  if (length(labels) == length(data_list)) {
    labels_expanded <- unlist(mapply(
      function(lbl, mat) rep(lbl, nrow(mat)),
      labels, data_list,
      SIMPLIFY = FALSE
    ), use.names = FALSE)
  } else if (length(labels) == total_rows) {
    labels_expanded <- labels
  } else {
    stop("`labels` must be either session-level (length(data_list)) or sample-level (sum(nrow(data_list))).")
  }
  data_standardized <- scale(data, center, scale)
  pca_result <- prcomp(data_standardized)
  num_components <- min(ncol(pca_result$x), max(1L, ceiling(ncol(pca_result$x) * k)))
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(variance_explained)
  pca_data <- pca_result$x
  compute_fisher_score <- function(data, labels) {
    n_classes <- length(unique(labels))
    overall_mean <- colMeans(data)
    fisher_scores <- numeric(ncol(data))

    for (j in 1:ncol(data)) {
      between_class_var <- 0
      within_class_var <- 0

      for (class in unique(labels)) {
        class_data <- data[labels == class, j, drop = FALSE]
        class_mean <- mean(class_data)
        class_size <- length(class_data)

        between_class_var <- between_class_var + class_size * (class_mean - overall_mean[j])^2
        within_class_var <- within_class_var + sum((class_data - class_mean)^2)
      }

      # Fisher criterion: maximize between-class variance / within-class variance.
      fisher_scores[j] <- between_class_var / max(within_class_var, 1e-12)
    }

    return(fisher_scores)
  }

  fisher_scores <- compute_fisher_score(pca_data, labels_expanded)
  fisher_df <- data.frame(Principal_Component = 1:ncol(pca_data), Fisher_Score = fisher_scores)
  fisher_df <- fisher_df[order(fisher_df$Fisher_Score, decreasing = TRUE), ]
  fisher_df <- fisher_df[1:num_components, ]

  return(list(fisher_df = fisher_df, pca_result = pca_result,
              selected_indices = fisher_df[,1],
              variance_explained = variance_explained,
              cumulative_variance = cumulative_variance))
}
