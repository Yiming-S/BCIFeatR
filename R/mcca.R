# Module: Regularized SUMCOR multiset CCA for feature-view fusion.
# Scope: Learns per-view linear projections that maximize summed pairwise
# correlation across views, then projects new views into aligned canonical
# scores. Operates on already-extracted feature matrices (one per view).

#' Train a regularized SUMCOR multiset CCA (MCCA) over feature views.
#'
#' Each view is a feature matrix on the same trials (rows aligned across views).
#' MCCA finds per-view linear projections whose projected scores are maximally
#' correlated across views. Under the SUMCOR formulation this is the generalized
#' eigenproblem `C w = lambda B w`, where `C` is the full stacked covariance and
#' `B` is its block-diagonal (within-view) part with ridge regularization. It is
#' solved by Cholesky-whitening `B` (the same device used elsewhere in the
#' package) and taking the leading eigenvectors.
#'
#' @param views List of numeric matrices, one per view, each
#'   `n_trials x d_view`, all sharing the same `n_trials`.
#' @param ncomp Number of canonical components to keep.
#' @param reg Ridge added to each within-view block, relative to its mean
#'   diagonal, for numerical stability.
#' @param center,scale Whether to center/scale each view; the statistics are
#'   stored for [mcca_transform()].
#' @return An object of class `"mcca"` with per-view `weights`, the stored
#'   `center`/`scale`, view `dims`, `ncomp`, `reg`, kept `eigenvalues`, and
#'   `n_views`.
#' @seealso [mcca_transform()]
#' @export
mcca_train <- function(views, ncomp = 1L, reg = 1e-3,
                       center = TRUE, scale = TRUE) {
  if (!is.list(views) || length(views) < 2L) {
    stop("`views` must be a list of at least two matrices.")
  }
  views <- lapply(views, as.matrix)
  n <- nrow(views[[1]])
  if (any(vapply(views, nrow, integer(1)) != n)) {
    stop("All views must have the same number of rows (aligned trials).")
  }
  V <- length(views)
  dims <- vapply(views, ncol, integer(1))
  ncomp <- as.integer(ncomp)
  if (!is.finite(ncomp) || ncomp < 1L) stop("`ncomp` must be a positive integer.")
  ncomp <- min(ncomp, min(dims))

  centers <- vector("list", V)
  scales  <- vector("list", V)
  Zs      <- vector("list", V)
  for (v in seq_len(V)) {
    Xv <- views[[v]]
    cv <- if (isTRUE(center)) colMeans(Xv) else rep(0, dims[v])
    Xc <- sweep(Xv, 2, cv, "-")
    sv <- if (isTRUE(scale)) apply(Xc, 2, stats::sd) else rep(1, dims[v])
    sv[!is.finite(sv) | sv <= 0] <- 1
    Zs[[v]] <- sweep(Xc, 2, sv, "/")
    centers[[v]] <- cv
    scales[[v]]  <- sv
  }

  Z <- do.call(cbind, Zs)
  D <- ncol(Z)
  C <- crossprod(Z) / n
  C <- (C + t(C)) / 2

  # Block-diagonal within-view covariance with ridge regularization.
  offs <- c(0L, cumsum(dims))
  B <- matrix(0, D, D)
  for (v in seq_len(V)) {
    rng <- (offs[v] + 1L):offs[v + 1L]
    blk <- C[rng, rng, drop = FALSE]
    ridge <- reg * mean(diag(blk))
    if (!is.finite(ridge) || ridge <= 0) ridge <- reg
    diag(blk) <- diag(blk) + ridge
    B[rng, rng] <- blk
  }

  # Solve C w = lambda B w via Cholesky whitening of B.
  U <- chol((B + t(B)) / 2)
  invU <- backsolve(U, diag(D))
  M <- crossprod(invU, C) %*% invU
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  keep <- seq_len(ncomp)
  W <- invU %*% eig$vectors[, keep, drop = FALSE]    # D x ncomp stacked weights

  weights <- lapply(seq_len(V), function(v) {
    rng <- (offs[v] + 1L):offs[v + 1L]
    W[rng, , drop = FALSE]
  })

  structure(
    list(weights = weights, center = centers, scale = scales,
         dims = dims, ncomp = ncomp, reg = reg,
         eigenvalues = eig$values[keep], n_views = V),
    class = "mcca"
  )
}

#' Project feature views into aligned MCCA canonical scores.
#'
#' @param object An `"mcca"` object from [mcca_train()].
#' @param views List of view matrices matching the trained views' dimensions.
#' @return List with per-view `scores` (each `n_trials x ncomp`), the cross-view
#'   `consensus` (mean of per-view scores, `n_trials x ncomp`), and `concat`
#'   (column-bound per-view scores, `n_trials x (n_views * ncomp)`).
#' @seealso [mcca_train()]
#' @export
mcca_transform <- function(object, views) {
  if (!inherits(object, "mcca")) {
    stop("`object` must be an 'mcca' object from mcca_train().")
  }
  if (!is.list(views) || length(views) != object$n_views) {
    stop("`views` must be a list matching the number of trained views.")
  }
  V <- object$n_views
  scores <- vector("list", V)
  for (v in seq_len(V)) {
    Xv <- as.matrix(views[[v]])
    if (ncol(Xv) != object$dims[v]) {
      stop(sprintf("View %d has %d columns; expected %d.",
                   v, ncol(Xv), object$dims[v]))
    }
    Zc <- sweep(sweep(Xv, 2, object$center[[v]], "-"), 2, object$scale[[v]], "/")
    scores[[v]] <- Zc %*% object$weights[[v]]
  }
  consensus <- Reduce(`+`, scores) / V
  list(scores = scores, consensus = consensus, concat = do.call(cbind, scores))
}
