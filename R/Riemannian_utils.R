# Module: Riemannian geometry utilities on SPD matrices.
# Scope: Mean/log-exp maps, tangent-space conversion, geodesic filtering,
# and manifold-side variable selection used by TS/ACM/Riemannian features.

#' Compute a mean of SPD matrices under a chosen Riemannian metric.
#'
#' @param cov_matrices 3D array (p x p x n) or list of SPD matrices.
#' @param max_iterations Maximum number of fixed-point updates (affine-invariant only).
#' @param epsilon Convergence threshold and eigenvalue floor in log-domain steps.
#' @param metric One of `"euclid"` (arithmetic mean, fastest — default),
#'   `"logeuclid"` (matrix-log arithmetic mean, closed form), or
#'   `"riemann"` (affine-invariant Fréchet mean, iterative).
#' @return SPD matrix representing the manifold mean.
#' @export
riemannian_mean <- function(cov_matrices, max_iterations = 30, epsilon = 1e-5,
                            metric = c("euclid", "logeuclid", "riemann")) {
  metric <- match.arg(metric)
  if (is.list(cov_matrices)) cov_matrices <- simplify2array(cov_matrices)
  m <- dim(cov_matrices)[3]
  n <- dim(cov_matrices)[1]

  if (metric == "euclid") {
    P <- rowMeans(cov_matrices, dims = 2L)
    return((P + t(P)) / 2)
  }

  if (metric == "logeuclid") {
    # log-Euclidean: closed-form, one eigen per trial, no iteration.
    Lsum <- matrix(0, n, n)
    for (i in seq_len(m)) {
      Ci  <- (cov_matrices[, , i] + t(cov_matrices[, , i])) / 2
      ei  <- eigen(Ci, symmetric = TRUE)
      lam <- log(pmax(ei$values, epsilon))
      Lsum <- Lsum + ei$vectors %*% (t(ei$vectors) * lam)
    }
    Lmean <- Lsum / m
    Lmean <- (Lmean + t(Lmean)) / 2
    eL <- eigen(Lmean, symmetric = TRUE)
    out <- eL$vectors %*% (t(eL$vectors) * exp(pmin(eL$values, 50)))
    return((out + t(out)) / 2)
  }

  # metric == "riemann": affine-invariant Fréchet mean.
  P_omega <- rowMeans(cov_matrices, dims = 2L)
  for (iteration in seq_len(max_iterations)) {
    eig <- eigen(P_omega, symmetric = TRUE)
    safe_vals <- pmax(eig$values, epsilon)
    sqrt_P_omega     <- eig$vectors %*% (t(eig$vectors) * sqrt(safe_vals))
    inv_sqrt_P_omega <- eig$vectors %*% (t(eig$vectors) / sqrt(safe_vals))
    S <- matrix(0, n, n)
    for (i in seq_len(m)) {
      eig_i <- eigen(inv_sqrt_P_omega %*% cov_matrices[, , i] %*% inv_sqrt_P_omega,
                     symmetric = TRUE)
      S <- S + (eig_i$vectors %*% (t(eig_i$vectors) * log(pmax(eig_i$values, epsilon))))
    }
    S <- S / m
    # Convergence gauged directly by tangent-step Frobenius norm.
    if (sqrt(sum(S * S)) < epsilon) break
    eigS <- eigen(S, symmetric = TRUE)
    # Cap exp argument to prevent overflow on pathological trials.
    capped <- pmin(eigS$values, 50)
    P_omega <- sqrt_P_omega %*% eigS$vectors %*%
      (t(eigS$vectors) * exp(capped)) %*% sqrt_P_omega
  }
  (P_omega + t(P_omega)) / 2
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
  # Replace geigen(Sigma_b, Sigma_w) with a Cholesky-whitened standard eigen.
  # Solves Sigma_b v = lambda Sigma_w v by forming M = (U')^{-1} Sigma_b U^{-1}
  # where U is the upper-triangular Cholesky factor of Sigma_w (U' U = Sigma_w),
  # recovering eigenvectors v = U^{-1} w. Substantially faster than geigen for
  # the sizes produced by ACM/TS and keeps the hot path free of geigen deps.
  Sigma_w_sym <- (Sigma_w + t(Sigma_w)) / 2
  U <- tryCatch(chol(Sigma_w_sym), error = function(e) {
    diag(Sigma_w_sym) <- diag(Sigma_w_sym) + epsilon
    chol(Sigma_w_sym)
  })
  invU <- backsolve(U, diag(nrow(U)))
  M    <- crossprod(invU, Sigma_b) %*% invU
  M    <- (M + t(M)) / 2
  eig  <- eigen(M, symmetric = TRUE)
  ncomps <- min(nclass - 1, nfeats)
  keep   <- head(order(eig$values, decreasing = TRUE), ncomps)
  Wtilde <- invU %*% eig$vectors[, keep, drop = FALSE]
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
#' @param eps Eigenvalue floor for numerical stability.
#' @return Scalar geodesic distance.
#' @export
riemannian_distance <- function(A, B, eps = 1e-12) {
  # Whitening-based formulation: d(A,B)^2 = sum(log eig(B^{-1/2} A B^{-1/2})^2).
  # Avoids the generalized eigenproblem (geigen) which is markedly slower on
  # moderate-sized SPD matrices and often unavailable in lean runtimes.
  Bs   <- (B + t(B)) / 2
  eigB <- eigen(Bs, symmetric = TRUE)
  lamB <- pmax(eigB$values, eps)
  inv_sqrt_B <- eigB$vectors %*% (t(eigB$vectors) / sqrt(lamB))
  W <- inv_sqrt_B %*% ((A + t(A)) / 2) %*% inv_sqrt_B
  eW <- eigen((W + t(W)) / 2, symmetric = TRUE)
  sqrt(sum(log(pmax(eW$values, eps))^2))
}

#' Map SPD matrices to upper-triangular tangent-space vectors.
#'
#' @param cov_matrices List of SPD matrices.
#' @param epsilon Eigenvalue floor.
#' @param P_omega Optional manifold reference; if NULL it is estimated from data.
#' @param metric Mean metric passed to `riemannian_mean()` when `P_omega` is NULL.
#'   Defaults to `"euclid"` (arithmetic mean), which is substantially faster
#'   than the affine-invariant Fréchet mean with negligible accuracy loss for
#'   tangent-space classification. See `riemannian_mean()` for alternatives.
#' @return List with tangent coordinates `S` (d x n) and reference `P_omega`.
#' @export
map_2_tangent_space <- function(cov_matrices, epsilon = 1e-12, P_omega = NULL,
                                metric = c("euclid", "logeuclid", "riemann")) {
  metric <- match.arg(metric)
  stopifnot(is.list(cov_matrices))
  n_trials <- length(cov_matrices)
  p <- ncol(cov_matrices[[1]])
  d <- p * (p + 1) / 2

  if (is.null(P_omega)) {
    P_omega <- riemannian_mean(simplify2array(cov_matrices),
                               epsilon = epsilon, metric = metric)
  }
  eigP <- eigen((P_omega + t(P_omega)) / 2, symmetric = TRUE)
  lamP <- pmax(eigP$values, epsilon)
  VP   <- eigP$vectors
  Pstr <- list(
    sqrt     = VP %*% (t(VP) * sqrt(lamP)),
    inv_sqrt = VP %*% (t(VP) / sqrt(lamP))
  )

  # Pre-compute row-major upper-tri index and diagonal linear positions once.
  ij_all <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  ij <- ij_all[order(ij_all[, "row"], ij_all[, "col"]), , drop = FALSE]
  diag_pos <- which(ij[, "row"] == ij[, "col"])

  # Batched two-sided whitening: compute W_i = inv_sqrt %*% C_i %*% inv_sqrt
  # across all trials using two large GEMM calls, then eigen-decompose each
  # whitened slice individually. This collapses 2 x n_trials small GEMMs
  # into 2 large ones, which BLAS handles much faster on ACM-sized matrices.
  inv_sqrt <- Pstr$inv_sqrt
  C_arr <- simplify2array(cov_matrices)
  C_left <- inv_sqrt %*% matrix(C_arr, p, p * n_trials)       # (p, p*n)
  C_left_arr <- array(C_left, c(p, p, n_trials))
  # Right multiply by inv_sqrt: reshape so each trial's rows become columns,
  # then another single GEMM completes the two-sided whitening.
  C_left_perm <- aperm(C_left_arr, c(2, 1, 3))                # (p, p, n)
  rhs_mat <- inv_sqrt %*% matrix(C_left_perm, p, p * n_trials)
  W_arr <- aperm(array(rhs_mat, c(p, p, n_trials)), c(2, 1, 3))

  S <- matrix(NA_real_, nrow = d, ncol = n_trials)
  for (i in seq_len(n_trials)) {
    Wi <- (W_arr[, , i] + t(W_arr[, , i])) / 2
    eigW <- eigen(Wi, symmetric = TRUE)
    lam  <- log(pmax(eigW$values, epsilon))
    L    <- eigW$vectors %*% (t(eigW$vectors) * lam)
    vec  <- L[ij]
    # Off-diagonal coordinates are scaled by sqrt(2) under Frobenius inner product.
    vec[-diag_pos] <- vec[-diag_pos] * sqrt(2)
    S[, i] <- vec
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
  # Right singular vectors are unused downstream; skip them to save work on
  # tall/wide tangent feature matrices (e.g. ACM_TS produces very wide S).
  svd_result <- svd(S, nv = 0)
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
