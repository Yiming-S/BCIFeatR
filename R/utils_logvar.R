# Module: Log-variance feature utilities.
# Scope: Computes log-variance features directly and with optional PCA projection.

#' Compute column-wise sample variance without `apply()` overhead.
#'
#' `apply(X, 2, var)` is the idiomatic choice but ~5-10x slower than this
#' closed-form on matrices with many columns, which matters on FBCSP /
#' logvar hot paths that call it once per trial per band.
#'
#' @param X Numeric matrix with observations in rows.
#' @return Numeric vector of unbiased column variances (NA when n < 2).
#' @keywords internal
.col_vars <- function(X) {
  n <- nrow(X)
  if (n < 2L) return(rep(NA_real_, ncol(X)))
  m <- colMeans(X)
  (colMeans(X * X) - m * m) * (n / (n - 1))
}

#' Compute log-variance features for vector, matrix, or trial-list inputs.
#'
#' @param data Numeric vector, numeric matrix, or list of numeric matrices.
#' @return Numeric scalar/vector/matrix of log-variance values.
#' @export
log_var <- function(data) {
  if (is.numeric(data) && is.vector(data)) {
    out <- log(var(data))
  } else if (is.matrix(data)) {
    out <- log(.col_vars(data))
    p <- ncol(data)
    names(out) <- if (is.null(colnames(data))) paste0("x", 1:p) else colnames(data)
  } else {
    n <- length(data)
    p <- ncol(data[[1]])
    out <- matrix(, n, p)
    colnames(out) <- if (is.null(colnames(data[[1]]))) {
      paste0("x", 1:p)
    } else {
      colnames(data[[1]])
    }
    for (i in 1:n) {
      out[i, ] <- log(.col_vars(data[[i]]))
    }
  }
  # Clamp `-Inf` from zero-variance channels to a finite floor.
  infinite <- is.infinite(out)
  if (any(infinite)) out[infinite] <- log(1e-16)
  out
}

#' Compute log-variance after optional PCA fit/apply step.
#'
#' @param data Numeric matrix where rows are observations.
#' @param pca_model Optional trained PCA model when `apply_pca = TRUE`.
#' @param apply_pca Whether to apply a provided PCA model.
#' @param n_components Optional number of principal components to keep.
#' @return List with `logvar` vector and possibly fitted `pca_model`.
#' @export
logvar_transform <- function(data, 
                             pca_model = NULL, 
                             apply_pca = FALSE, 
                             n_components = NULL) {
  if (apply_pca) {
    if (is.null(pca_model)) {
      stop("A trained pca_model must be provided to apply PCA transformation to test data.")
    }
    data <- scale(data, center = pca_model$center, scale = pca_model$scale)
    data <- data %*% pca_model$rotation
    if (!is.null(n_components)) {
      data <- data[, 1:n_components, drop = FALSE]
    }
  } else {
    if (!is.null(n_components)) {
      # Fit PCA on training matrix and keep leading components.
      pca_model <- prcomp(data, center = TRUE, scale. = TRUE)
      data <- pca_model$x[, 1:n_components, drop = FALSE]
    }
  }
  variances <- apply(data, 2, var)
  log_variances <- log(variances + 1e-10)
  
  return(list(logvar = log_variances, pca_model = pca_model))
}
