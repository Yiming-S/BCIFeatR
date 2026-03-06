# Module: Riemannian geometry utilities on SPD matrices.
# Scope: Mean/log-exp maps, tangent-space conversion, geodesic filtering,
# and manifold-side variable selection used by TS/ACM/Riemannian features.

#' Compute the affine-invariant Riemannian mean of SPD matrices.
#'
#' @param cov_matrices 3D array (p x p x n) or list of SPD matrices.
#' @param max_iterations Maximum number of fixed-point updates.
#' @param epsilon Convergence threshold and eigenvalue floor in log-domain steps.
#' @return SPD matrix representing the manifold mean.
#' @export
riemannian_mean <- function(cov_matrices, max_iterations = 500, epsilon = 1e-5) {
  if (is.list(cov_matrices)) cov_matrices <- simplify2array(cov_matrices)
  m <- dim(cov_matrices)[3]
  n <- dim(cov_matrices)[1]
  P_omega <- rowMeans(cov_matrices, dims = 2L)
  
  for (iteration in 1:max_iterations) {
    # Move to local tangent space at current estimate.
    eig <- eigen(P_omega, symmetric = TRUE)
    safe_vals <- pmax(eig$values, epsilon)
    sqrt_P_omega     <- eig$vectors %*% (t(eig$vectors) * sqrt(safe_vals))
    inv_sqrt_P_omega <- eig$vectors %*% (t(eig$vectors) / sqrt(safe_vals))
    S <- matrix(0, n, n)
    for (i in 1:m) {
      eig_i <- eigen(inv_sqrt_P_omega %*% cov_matrices[, , i] %*% inv_sqrt_P_omega,
                     symmetric = TRUE)
      S <- S + (eig_i$vectors %*% (t(eig_i$vectors) * log(pmax(eig_i$values, epsilon))))
    }
    S <- S / m
    PS <- P_omega %*% S
    # Exponentiate averaged tangent direction back to SPD manifold.
    eigS <- eigen(S, symmetric = TRUE)
    P_omega <- sqrt_P_omega %*% eigS$vectors %*%
      (t(eigS$vectors) * exp(eigS$values)) %*% sqrt_P_omega
    
    if (sqrt(sum(PS * t(PS))) < epsilon) break
  }
  P_omega
}

#' Apply the affine-invariant logarithmic map to an SPD matrix.
#'
#' @param P_omega Reference SPD matrix or list containing `inv_sqrt`.
#' @param P_i SPD matrix to be mapped to tangent space at `P_omega`.
#' @param eps Eigenvalue floor for numerical stability.
#' @return Symmetric tangent matrix.
#' @export
log_map <- function(P_omega, P_i, eps = 1e-12) {
  if (is.list(P_omega)) {
    inv_sqrt_P <- P_omega$inv_sqrt
  } else {
    eigP <- eigen((P_omega + t(P_omega)) / 2, symmetric = TRUE)
    lamP <- pmax(eigP$values, eps)
    VP   <- eigP$vectors
    inv_sqrt_P <- VP %*% (t(VP) / sqrt(lamP))
  }
  Pi   <- (P_i + t(P_i)) / 2
  W    <- inv_sqrt_P %*% Pi %*% inv_sqrt_P
  eigW <- eigen((W + t(W)) / 2, symmetric = TRUE)
  lamW <- pmax(eigW$values, eps)
  VW   <- eigW$vectors
  L <- VW %*% (t(VW) * log(lamW))
  (L + t(L)) / 2
}

#' Apply the affine-invariant exponential map from tangent to SPD space.
#'
#' @param P_omega Reference SPD matrix or list with `sqrt` and `inv_sqrt`.
#' @param S_i Tangent matrix at `P_omega`.
#' @return SPD matrix on the manifold.
#' @export
exp_map <- function(P_omega, S_i) {
  if (is.list(P_omega)) {
    sqrt_P_omega     <- P_omega[["sqrt"]]
    inv_sqrt_P_omega <- P_omega[["inv_sqrt"]]
  } else {
    eig <- eigen(P_omega, symmetric = TRUE)
    safe_vals <- pmax(eig$values, 1e-12)
    sqrt_P_omega     <- eig$vectors %*% (t(eig$vectors) * sqrt(safe_vals))
    inv_sqrt_P_omega <- eig$vectors %*% (t(eig$vectors) / sqrt(safe_vals))
  }
  w     <- inv_sqrt_P_omega %*% S_i %*% inv_sqrt_P_omega
  eig_w <- eigen(w, symmetric = TRUE)
  exp_w <- eig_w$vectors %*% (t(eig_w$vectors) * exp(eig_w$values))
  sqrt_P_omega %*% exp_w %*% sqrt_P_omega
}

#' Estimate FGDA filter weights in tangent space.
#'
#' @param cov_matrices List of SPD matrices.
#' @param labels Class labels aligned with `cov_matrices`.
#' @param epsilon Numerical floor.
#' @return List with `W_tilde` selector weights and manifold reference `P_omega`.
#' @export
fgda_filters <- function(cov_matrices, labels, epsilon = 1e-6) {
  
  ntrials <- length(labels)
  nchans  <- ncol(cov_matrices[[1]])
  nfeats  <- choose(nchans + 1, 2)
  P_omega <- riemannian_mean(simplify2array(cov_matrices))
  eig     <- eigen(P_omega, TRUE)
  sqrt_eigvals <- sqrt(pmax(eig$values, epsilon))
  eigvecs <- eig$vectors
  P_omega <- list(
    sqrt     = eigvecs %*% (t(eigvecs) * sqrt_eigvals),
    inv_sqrt = eigvecs %*% (t(eigvecs) / sqrt_eigvals)
  )
  S_tilde <- matrix(0, nfeats, ntrials)
  mask    <- upper.tri(diag(nchans), diag = TRUE)
  
  for (i in 1:ntrials) {
    # Flatten upper triangle of tangent matrices as feature coordinates.
    S_i <- log_map(P_omega, cov_matrices[[i]])
    S_tilde[, i] <- S_i[mask]
  }
  unique_labels <- unique(labels)
  nclass <- length(unique_labels)
  class_mean <- matrix(0, nfeats, nclass)
  for (k in 1:nclass) {
    class_mean[, k] <- rowMeans(S_tilde[, labels == unique_labels[k], drop = FALSE])
  }
  # Generalized eigenproblem between between-class and within-class scatter.
  Sigma_b <- cov(t(class_mean))
  S_centered <- S_tilde - class_mean[, match(labels, unique_labels)]
  Sigma_w <- tcrossprod(S_centered) / ntrials
  # Apply Ledoit-Wolf-style shrinkage directly to the scatter matrix.
  trSw <- sum(diag(Sigma_w))
  mu_w <- trSw / nfeats
  Sigma_w <- (1 - 1 / ntrials) * Sigma_w
  diag(Sigma_w) <- diag(Sigma_w) + (1 / ntrials) * mu_w
  eig <- geigen(Sigma_b, Sigma_w, TRUE)
  ncomps <- min(nclass - 1, nfeats)
  keep <- head(order(eig$values, decreasing = TRUE), ncomps)
  Wtilde <- eig$vectors[, keep, drop = FALSE]
  Wtilde <- sweep(Wtilde, 2, sqrt(colSums(Wtilde^2)), "/")
  
  list(W_tilde = Wtilde, P_omega = P_omega)
}

#' Apply geodesic feature filtering and map back to SPD.
#'
#' @param P_x Input SPD matrix.
#' @param P_omega Reference manifold point.
#' @param W_tilde Tangent-space selector weights.
#' @param epsilon Unused placeholder for interface compatibility.
#' @return Filtered SPD matrix.
#' @export
geodesic_filtering <- function(P_x, P_omega, W_tilde, epsilon = 1e-6) {
  S_x <- log_map(P_omega, P_x)
  mask <- upper.tri(S_x, diag = TRUE)
  vec_S_x <- S_x[mask]
  if (is.matrix(W_tilde)) {
    if (nrow(W_tilde) != length(vec_S_x)) {
      stop("W_tilde rows must match upper-tri length.")
    }
    W_tilde <- rowSums(W_tilde^2)
    max_w <- max(W_tilde)
    if (!is.finite(max_w) || max_w <= 0) {
      W_tilde[] <- 0
    } else {
      W_tilde <- W_tilde / max_w
    }
  } else if (length(W_tilde) != length(vec_S_x)) {
    stop("W_tilde length must match upper-tri length.")
  }
  W_tilde[!is.finite(W_tilde)] <- 0
  S_tilde_x <- vec_S_x * W_tilde
  n <- ncol(P_x)
  S_tilde_x_matrix <- matrix(0, n, n)
  S_tilde_x_matrix[mask] <- S_tilde_x
  S_tilde_x_matrix <- S_tilde_x_matrix + t(S_tilde_x_matrix)
  diag(S_tilde_x_matrix) <- diag(S_tilde_x_matrix) / 2
  frob <- sqrt(sum(S_tilde_x_matrix^2))
  if (frob > 0) S_tilde_x_matrix <- S_tilde_x_matrix / frob
  exp_map(P_omega, S_tilde_x_matrix)
}

#' Compute affine-invariant Riemannian distance between two SPD matrices.
#'
#' @param A First SPD matrix.
#' @param B Second SPD matrix.
#' @return Scalar geodesic distance.
#' @export
riemannian_distance <- function(A, B) {
	vals <- geigen(A, B, TRUE, TRUE)$values
	sqrt(sum(log(vals[vals > 0])^2))
}

#' Map SPD matrices to upper-triangular tangent-space vectors.
#'
#' @param cov_matrices List of SPD matrices.
#' @param epsilon Eigenvalue floor.
#' @param P_omega Optional manifold reference; if NULL it is estimated from data.
#' @return List with tangent coordinates `S` (d x n) and reference `P_omega`.
#' @export
map_2_tangent_space <- function(cov_matrices, epsilon = 1e-12, P_omega = NULL) {
  stopifnot(is.list(cov_matrices))
  n_trials <- length(cov_matrices)
  p <- ncol(cov_matrices[[1]])
  d <- p * (p + 1) / 2
  
  .triu_idx_rowmajor <- function(p) {
    idx <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
    idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
  }
  
  if (is.null(P_omega)) {
    P_omega <- riemannian_mean(simplify2array(cov_matrices), epsilon = epsilon)
  }
  eigP <- eigen((P_omega + t(P_omega)) / 2, symmetric = TRUE)
  lamP <- pmax(eigP$values, epsilon)
  VP   <- eigP$vectors
  Pstr <- list(
    sqrt     = VP %*% (t(VP) * sqrt(lamP)),
    inv_sqrt = VP %*% (t(VP) / sqrt(lamP))
  )
  
  S  <- matrix(NA_real_, nrow = d, ncol = n_trials)
  ij <- .triu_idx_rowmajor(p)
  
  for (i in seq_len(n_trials)) {
    L <- log_map(Pstr, cov_matrices[[i]], eps = epsilon)
    # Off-diagonal coordinates are scaled by sqrt(2) under Frobenius inner product.
    diag_mask <- row(L) == col(L)
    L[!diag_mask] <- L[!diag_mask] * sqrt(2)
    S[, i] <- L[ij]
  }
  
  list(S = S, P_omega = P_omega)
}

#' Perform SVD-weighted BH selection on tangent coordinates.
#'
#' @param S Tangent feature matrix (d x n).
#' @param fdr_level Benjamini-Hochberg threshold.
#' @return List with filtered `S` and selected feature `indices`.
#' @export
variable_selection <- function(S, fdr_level = .05) {
  svd_result <- svd(S)
  U <- svd_result$u
  Lambda <- svd_result$d
  
  S_o <- t(U) %*% S
  p_values <- apply(S_o, 1, function(x) {
    tryCatch(t.test(x)$p.value, error = function(e) 1)
  })
  
  weights <- Lambda / sum(Lambda)
  adjusted_p <- p.adjust(p_values, method = "BH") * weights
  
  selected_indices <- which(adjusted_p < fdr_level)
  list(S = S[selected_indices, , drop = FALSE], indices = selected_indices)
}

#' Select manifold variables with adaptive FDR escalation.
#'
#' @param cov_matrices List of SPD matrices.
#' @param fdr_level Initial FDR level.
#' @param max_fdr Maximum FDR allowed during escalation.
#' @param min_features Minimum number of selected features.
#' @param epsilon Numerical floor.
#' @return List with reference `P_omega` and binary selector `W_tilde`.
#' @export
riem_var_select <- function(cov_matrices, fdr_level = .05,
                            max_fdr = 0.9, min_features = 3,
                            epsilon = 1e-6) {
  map_result <- map_2_tangent_space(cov_matrices, epsilon = epsilon)
  S <- map_result$S
  P_omega <- map_result$P_omega
  
  repeat {
    var_select_result <- variable_selection(S, fdr_level)
    selected_indices  <- var_select_result$indices
    # Increase FDR until enough coordinates survive or max_fdr is reached.
    if (length(selected_indices) >= min_features) break
    if (fdr_level >= max_fdr) {
      warning("Not enough features selected at max_fdr; using all features.")
      selected_indices <- seq_len(nrow(S))
      break
    }
    fdr_level <- min(fdr_level * 2, max_fdr)
    warning(sprintf("Not enough features; increasing fdr_level to %.3f and retrying.", fdr_level))
  }
  
  W_tilde <- numeric(nrow(S))
  W_tilde[selected_indices] <- 1
  
  list(P_omega = P_omega, W_tilde = W_tilde)
}
