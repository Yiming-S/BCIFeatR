# Module: Markov-Switching VAR (MSVAR) estimation engine (vendored).
# Scope: EM estimation (Hamilton/Kim switching Kalman filter-smoother E-step,
# closed-form weighted-least-squares M-step) of Markov-switching vector
# autoregressions, plus initialization and Viterbi decoding. Powers the
# "MSVAR" feature and the msvar_train()/msvar_predict() classifier.
#
# Provenance: adapted from the Degras Lab 'msvar' source by David Degras,
# Chee-Ming Ting, and Hernando Ombao. Method reference:
#   D. Degras, C.-M. Ting, H. Ombao (2022). "Markov-switching state-space
#   models with applications to neuroimaging." Computational Statistics &
#   Data Analysis.
# Vendored changes: MASS::ginv replaced by an internal SVD pseudo-inverse
# (.ginv); the optional foreach/%dopar% parallel path removed (sequential
# only); an uninitialized `LLflag` initialized to 0; a stray `w` -> `weights`
# fix in init_msvar; control parameters extracted explicitly instead of via
# assign(). All functions here are package-internal except fit_var(), which is
# exported.

# Moore-Penrose pseudo-inverse via SVD (drop-in for MASS::ginv).
.ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  s <- svd(X)
  pos <- s$d > max(tol * s$d[1], 0)
  if (!any(pos)) return(matrix(0, ncol(X), nrow(X)))
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

# Companion-matrix stability test for a VAR coefficient block.
check_stable_var <- function(A, p) {
  n <- ncol(A)
  companion_matrix <- matrix(0, n * p, n * p)
  companion_matrix[1:n, ] <- A
  if (p > 1) {
    companion_matrix[(n + 1):(n * p), 1:(n * (p - 1))] <- diag(n * (p - 1))
  }
  if (any(is.infinite(companion_matrix) | is.na(companion_matrix))) return(FALSE)
  eigenvalues <- eigen(companion_matrix, only.values = TRUE)$values
  all(Mod(eigenvalues) < 1)
}

# Random contiguous segments [start, end] within 1:tt (length-biased sampling).
rand_segments <- function(tt, minlen = 1, maxlen = tt, nsegments = 100) {
  tt <- as.integer(tt)
  minlen <- as.integer(minlen)
  maxlen <- as.integer(maxlen)
  prob <- (tt - minlen + 1):(tt - maxlen + 1)
  len <- if (minlen == maxlen) minlen else sample(minlen:maxlen, nsegments, TRUE, prob)
  start <- ceiling(runif(nsegments) * (tt - len + 1))
  start[start == 0] <- 1
  start <- as.integer(start)
  end <- start + len - 1
  cbind(start, end)
}

#' Fit a least-squares vector autoregression VAR(p).
#'
#' Ordinary-least-squares VAR(p) estimate of the stacked coefficient matrix,
#' innovation covariance, and mean of a multivariate series, with a
#' pseudo-inverse fallback for rank-deficient designs and an eigenvalue floor on
#' the innovation covariance.
#'
#' @param y Numeric matrix (`samples x channels`); a vector is treated as a
#'   single channel.
#' @param p Autoregressive order; a positive integer with `p < nrow(y)`.
#' @param subset Optional integer index of rows to use as regression targets;
#'   defaults to all rows after the first `p`.
#' @return List with `theta` (`mu` length-`n`, `A` `n x (n*p)` stacked
#'   coefficients, `Q` `n x n` innovation covariance), the log-likelihood `LL`,
#'   and pseudo-`R2`; with attributes `model = "var"` and `length = nrow(y)`.
#'   Returns `NULL` when the usable target set is empty.
#' @references D. Degras, C.-M. Ting, H. Ombao (2022). Markov-switching
#'   state-space models with applications to neuroimaging. CSDA.
#' @export
fit_var <- function(y, p, subset = NULL) {
  p <- as.integer(p)
  stopifnot(p > 0 & p < NROW(y))
  if (!is.matrix(y)) dim(y) <- c(length(y), 1)
  dimy <- dim(y)
  n <- dimy[2]

  if (is.null(subset)) {
    subset0 <- 1:p
    subset <- (p + 1):nrow(y)
  } else {
    subset <- as.integer(subset)
    subset <- subset[(subset > 0) & (subset <= nrow(y))]
    subset0 <- subset[subset <= p]
    subset <- subset[subset > p]
  }
  tt0 <- length(subset0)
  tt <- length(subset)
  if (tt == 0) return(NULL)

  x <- matrix(0, tt, n * p)
  for (l in 1:p) x[, ((l - 1) * n + 1):(l * n)] <- y[subset - l, ]
  y0 <- y[subset0, , drop = FALSE]
  y <- y[subset, , drop = FALSE]
  ybar <- colMeans(y)
  xbar <- colMeans(x)

  yty <- crossprod(y) / tt - tcrossprod(ybar)
  xtx <- crossprod(x) / tt - tcrossprod(xbar)
  xty <- crossprod(x, y) / tt - tcrossprod(xbar, ybar)
  A <- tryCatch(t(solve(xtx, xty)), error = function(e) NULL)
  if (is.null(A)) {
    eig <- eigen(xtx, TRUE)
    eigvals <- eig$values
    lb <- max(1e-14, 1e-8 * eigvals[1])
    keep <- which(eigvals > lb)
    if (length(keep) == 0) {
      A <- matrix(0, n, n * p)
    } else {
      eigvals <- eigvals[keep]
      eigvecs <- eig$vectors[, keep]
      A <- crossprod(xty, eigvecs) %*% (t(eigvecs) / eigvals)
    }
  }

  Q <- yty + A %*% tcrossprod(xtx, A) - (A %*% xty) - t(A %*% xty)
  eig <- eigen(Q, TRUE)
  eigvals <- eig$values
  lb <- max(1e-14, 1e-8 * eigvals[1])
  if (eigvals[n] < lb) {
    eigvals <- pmax(eigvals, lb)
    eigvecs <- eig$vectors
    Q <- eigvecs %*% (eigvals * t(eigvecs))
  }

  mu <- ybar - A %*% xbar
  logdetQ <- sum(log(eigvals))
  ss0 <- if (tt0 == 0) 0 else sum(y0^2)
  nLL <- (n * (tt + tt0) * log(2 * pi)) + ss0 + (tt * logdetQ) + (n * tt)
  sse <- ss0 + tt * sum(diag(Q))
  if (tt0 == 0) {
    sst <- sum(y^2) - tt * sum(ybar^2)
  } else {
    ybar <- (tt0 / (tt0 + tt)) * colMeans(y0) + (tt / (tt0 + tt)) * ybar
    sst <- sum(y0^2, y^2) - (tt0 + tt) * sum(ybar^2)
  }

  out <- list(theta = list(mu = mu, A = A, Q = Q), LL = -nLL / 2, R2 = 1 - sse / sst)
  attr(out, "model") <- "var"
  attr(out, "length") <- dimy[1]
  out
}

# M-step VAR stability regularizer: shrink coefficient eigenvalues toward the
# unit circle, keeping the candidate that lowers the weighted-LS objective.
regularize_A <- function(A, eigvals, sum_MPb, sum_MCP, A_old, Q_old) {
  n <- nrow(A)
  p <- ncol(A) / n
  A_best <- A_old
  invQA <- solve(Q_old, A_old)
  objective_best <- sum(invQA * (A_old %*% sum_MPb)) - 2 * sum(invQA * sum_MCP)

  cc <- .99 / max(abs(eigvals))
  Aj <- A
  dim(Aj) <- c(n, n, p)
  for (l in 1:p) Aj[, , l] <- cc^l * Aj[, , l]
  dim(Aj) <- c(n, n * p)
  invQA <- solve(Q_old, Aj)
  objective <- sum(invQA * (Aj %*% sum_MPb)) - 2 * sum(invQA * sum_MCP)
  if (objective < objective_best) {
    objective_best <- objective
    A_best <- Aj
  }

  if (p == 1) {
    eigA <- eigen(A)
    nrm <- abs(eigA$values)
    shrink <- (nrm > .99)
    eigA$values[shrink] <- eigA$values[shrink] * (.99 / nrm[shrink])
    Aj <- eigA$vectors %*% diag(eigA$values) %*% solve(eigA$vectors)
    Aj <- Re(Aj)
    invQA <- solve(Q_old, Aj)
    objective <- sum(invQA * (Aj %*% sum_MPb)) - 2 * sum(invQA * sum_MCP)
    if (objective < objective_best) {
      objective_best <- objective
      A_best <- Aj
    }
  }
  A_best
}

# M-step: closed-form regime-wise updates of (mu, A, Q, Pi, Z) from the
# weighted sufficient statistics accumulated in the E-step.
M_var <- function(Ms, sum_MCP, sum_MP, sum_MPb, sum_Ms2, xbar, ybar,
                  theta_prev = NULL, fixed, init, verbose) {
  M <- ncol(Ms)
  tt <- nrow(Ms)
  n <- nrow(sum_MP)
  p <- ncol(sum_MCP) / n

  if (is.null(theta_prev)) {
    Abig <- array(0, dim = c(n, n, p, M))
    Abig[, , 1, ] <- rep(diag(.99, n), M)
    theta_prev <- list(A = Abig,
                       Q = array(rep(diag(n), M), c(n, n, M)),
                       mu = matrix(0, n, M), Pi = rep(1 / M, M),
                       Z = matrix(1 / M, M, M))
  }

  if (is.null(fixed$A)) {
    A <- array(0, dim = c(n, n * p, M))
    for (j in 1:M) {
      Aj <- tryCatch(t(solve(sum_MPb[, , j], t(sum_MCP[, , j]))),
                     error = function(e) NULL)
      if (is.null(Aj) || any(is.nan(Aj) | is.infinite(Aj)))
        Aj <- sum_MCP[, , j] %*% .ginv(sum_MPb[, , j])
      A[, , j] <- Aj
    }
    Abig <- matrix(0, n * p, n * p)
    if (p > 1) Abig[cbind((n + 1):(n * p), 1:(n * (p - 1)))] <- 1
    for (j in 1:M) {
      Abig[1:n, ] <- A[, , j]
      eigvals <- abs(eigen(Abig, FALSE, TRUE)$values)
      if (any(abs(eigvals) > .99)) {
        if (verbose) {
          warning(paste0("Eigenvalues of A", j, " greater than 0.99. Regularizing."),
                  immediate. = TRUE)
        }
        A_old <- matrix(theta_prev$A[, , , j], n, n * p)
        Q_old <- theta_prev$Q[, , j]
        A[, , j] <- regularize_A(A[, , j], eigvals, sum_MPb[, , j],
                                 sum_MCP[, , j], A_old, Q_old)
      }
    }
  } else {
    A <- fixed$A
  }

  if (is.null(fixed$Q)) {
    Q <- array(0, dim = c(n, n, M))
    sum_M <- colSums(Ms[(p + 1):tt, ])
    for (j in 1:M) {
      if (sum_M[j] == 0) { Q[, , j] <- diag(n); next }
      Aj <- A[, , j]
      Qj <- (sum_MP[, , j] - tcrossprod(sum_MCP[, , j], Aj) -
               tcrossprod(Aj, sum_MCP[, , j]) +
               Aj %*% sum_MPb[, , j] %*% t(Aj)) / sum_M[j]
      Q[, , j] <- 0.5 * (Qj + t(Qj))
    }
    for (j in 1:M) {
      eig <- eigen(Q[, , j], TRUE)
      eigvals <- abs(eig$values)
      lb <- 1e-8 * max(1, eigvals)
      if (min(eigvals) < lb) {
        eigvals[eigvals < lb] <- lb
        Q[, , j] <- eig$vectors %*% (eigvals * t(eig$vectors))
      }
    }
  } else {
    Q <- fixed$Q
  }

  if (is.null(fixed$mu)) {
    mu <- matrix(0, n, M)
    for (j in 1:M) mu[, j] <- ybar[, j] - A[, , j] %*% xbar[, j]
  } else {
    mu <- fixed$mu
  }

  t0 <- switch(init, mvn = 1, conditional = p + 1)
  Pi <- if (is.null(fixed$Pi)) Ms[t0, ] else fixed$Pi

  if (is.null(fixed$Z)) {
    Z <- sum_Ms2 / rowSums(sum_Ms2)
    Z[is.nan(Z)] <- 1 / M
  } else {
    Z <- fixed$Z
  }

  dim(A) <- c(n, n, p, M)
  list(mu = mu, A = A, Q = Q, Pi = Pi, Z = Z)
}

# Initialization: random-segment VAR fits clustered into M regimes, then one
# M-step to assemble a starting theta. RNG-driven (rand_segments + clustering).
init_msvar <- function(y, M, p, minlen = NULL, maxlen = NULL,
                       nseg = 100L, minclustsize = 1L, cluster = c("mu", "A", "Q"),
                       rmv_outliers = TRUE) {
  cluster <- match.arg(cluster, several.ok = TRUE)
  if (!is.matrix(y)) dim(y) <- c(length(y), 1)
  tt <- nrow(y)
  n <- ncol(y)

  if (is.null(minlen)) minlen <- min(5 * n * p, tt)
  if (is.null(maxlen)) maxlen <- tt
  stopifnot((minlen > 0) && (minlen <= maxlen) && (maxlen <= nrow(y)))

  nseg <- max(nseg, 2 * M)
  segments <- rand_segments(tt, minlen, maxlen, nseg)

  mu <- A <- Q <- matrix(, nseg, n)
  mask <- seq(0, by = n, len = n) + (1:n)
  stable <- logical(nseg)
  for (k in 1:nseg) {
    idx <- segments[k, 1]:segments[k, 2]
    thetak <- fit_var(y, p, idx)$theta
    if (rmv_outliers) stable[k] <- check_stable_var(thetak$A, 1)
    mu[k, ] <- thetak$mu
    A[k, ] <- thetak$A[mask]
    Q[k, ] <- thetak$Q[mask]
  }

  if (rmv_outliers) {
    err <- paste("Too few (less than 'M') segments to cluster",
                 "after removing those with unstable or outlying estimates.",
                 "\nConsider increasing 'minlen', decreasing 'M' or 'p',",
                 "or setting rmv_outliers = FALSE.")
    keep <- which(stable)
    if (length(keep) < M) stop(err)
    Amax <- apply(abs(A[keep, ]), 1, max)
    keep <- keep[Amax <= 3 * median(Amax)]
    if (length(keep) < M) stop(err)
    Qmax <- apply(abs(Q[keep, ]), 1, max)
    keep <- keep[Qmax <= 3 * median(Qmax)]
    if (length(keep) < M) stop(err)
    mu <- mu[keep, ]; A <- A[keep, ]; Q <- Q[keep, ]
    segments <- segments[keep, ]
  }

  xclust <- NULL
  if ("mu" %in% cluster) xclust <- cbind(xclust, mu)
  if ("A" %in% cluster) xclust <- cbind(xclust, A)
  if ("Q" %in% cluster) xclust <- cbind(xclust, Q)
  diss <- dist(xclust, "manhattan")
  dendro <- hclust(diss)
  success <- FALSE
  idx <- NULL
  for (k in M:nrow(segments)) {
    cluster <- cutree(dendro, k)
    freq <- table(cluster)
    idx <- which(freq >= minclustsize)
    if (length(idx) == M) {
      cluster[!is.element(cluster, idx)] <- NA
      success <- TRUE
      break
    }
  }
  if (!success) {
    warning("Could not obtain 'M' clusters with size at least 'minclustsize' ",
            "by hierarchical clustering.")
    cluster <- cutree(dendro, M)
    cluster[!is.element(cluster, idx)] <- NA
  }

  segments <- segments[!is.na(cluster), ]
  cluster <- cluster[!is.na(cluster)]
  cluster <- match(cluster, unique(cluster))
  weights <- matrix(0, tt, M)
  for (j in 1:M) {
    idxj <- which(cluster == j)
    idxt <- mapply(seq, from = segments[idxj, 1], to = segments[idxj, 2])
    idxt <- factor(unlist(idxt), levels = 1:tt)
    weights[, j] <- as.numeric(table(idxt))
  }
  weights <- weights / rowSums(weights)
  weights[is.na(weights)] <- 0
  wpos <- (rowSums(weights[1:(p + 1), ]) > 0)
  if (!all(wpos)) {
    t1 <- min(segments)
    if (t1 > p + 10) {
      weights[1:(p + 1), ] <- 1 / M
    } else if (t1 > p) {
      weights[1:t1, ] <- matrix(weights[t1, ], t1, M, byrow = TRUE)
    } else {
      zero <- which(!wpos)
      wpos <- which(wpos)
      weights[zero, ] <- matrix(colMeans(weights[wpos, , drop = FALSE]),
                                length(zero), M, byrow = TRUE)
    }
  }

  sum_Ms2 <- crossprod(weights[1:(tt - p), ], weights[(p + 1):tt, ])
  xmat <- matrix(0, tt - p, n * p)
  for (k in 1:p) xmat[, ((k - 1) * n + 1):(k * n)] <- y[(p + 1 - k):(tt - k), ]
  yy <- y[(p + 1):tt, ]
  yty <- array(dim = c(n, n, M))
  xtx <- array(dim = c(n * p, n * p, M))
  ytx <- array(dim = c(n, n * p, M))
  xbar <- matrix(, n * p, M)
  ybar <- matrix(, n, M)
  for (j in 1:M) {
    wj <- weights[(p + 1):tt, j]
    xbar[, j] <- crossprod(xmat, wj) / sum(wj)
    ybar[, j] <- crossprod(yy, wj) / sum(wj)
    xj <- (xmat - matrix(xbar[, j], tt - p, n * p, byrow = TRUE)) * sqrt(wj)
    yj <- (yy - matrix(ybar[, j], tt - p, n, byrow = TRUE)) * sqrt(wj)
    yty[, , j] <- crossprod(yj)
    xtx[, , j] <- crossprod(xj)
    ytx[, , j] <- crossprod(yj, xj)
  }
  theta <- M_var(Ms = weights, sum_MCP = ytx, sum_MP = yty,
                 sum_MPb = xtx, sum_Ms2 = sum_Ms2, xbar = xbar, ybar = ybar,
                 fixed = NULL, init = "mvn", verbose = FALSE)
  dim(theta$A) <- c(n, n, p, M)
  c(theta, list(weights = weights, clustsize = table(cluster)))
}

# E-step: Kim/Hamilton switching Kalman filter + smoother for one series.
skfs_var <- function(y, M, p, theta, init) {
  tt <- nrow(y)
  n <- ncol(y)
  mu <- theta$mu; A <- theta$A; Q <- theta$Q; Pi <- theta$Pi; Z <- theta$Z
  Zt <- t(Z)

  Mf <- matrix(0, tt, M)
  Ms <- matrix(0, tt, M)
  sum_Ms2 <- matrix(0, M, M)
  Lp <- matrix(0, tt, M)
  LL <- numeric(tt)

  invSqrtQ <- apply(Q, 3, function(x) solve(chol(x)))
  dim(invSqrtQ) <- c(n, n, M)
  logDetInvSqrtQ <- apply(invSqrtQ, 3, function(x) sum(log(diag(x))))
  logDetInvSqrtQ <- logDetInvSqrtQ - (n / 2) * log(2 * pi)

  for (j in 1:M) {
    ej <- y[1:p, , drop = FALSE]
    Lp[1:p, j] <- sqrt((2 * pi)^n) * exp(-0.5 * rowSums(ej^2))
    ej <- y[(p + 1):tt, ] - matrix(mu[, j], tt - p, n, byrow = TRUE)
    for (k in 1:p) ej <- ej - tcrossprod(y[(p + 1 - k):(tt - k), ], A[, , k, j])
    ej <- ej %*% invSqrtQ[, , j]
    Lp[(p + 1):tt, j] <- exp(-0.5 * rowSums(ej^2) + logDetInvSqrtQ[j])
  }

  t0 <- switch(init, mvn = 1L, conditional = p + 1L)
  Acc <- Pi * Lp[t0, ]
  LL[t0] <- log(sum(Acc))
  if (all(Acc == 0)) Acc <- rep(1, M)
  Mf[t0, ] <- Acc / sum(Acc)

  for (t in (t0 + 1):tt) {
    Acc <- Lp[t, ] * (Zt %*% Mf[t - 1, ])
    LL[t] <- log(sum(Acc))
    if (all(Acc == 0)) Acc <- rep(1, M)
    Mf[t, ] <- Acc / sum(Acc)
  }

  test <- is.infinite(LL)
  if (all(test)) LL <- -Inf else if (any(test)) LL[test] <- min(LL[!test])

  Ms[tt, ] <- Mf[tt, ]
  for (t in (tt - 1):t0) {
    U <- Mf[t, ] * Z
    U <- U / rep(colSums(U), each = M)
    U[is.nan(U)] <- 0
    Ms2 <- U * rep(Ms[t + 1, ], each = M)
    if (all(Ms2 == 0)) Ms2 <- matrix(1 / M^2, M, M)
    Ms2 <- Ms2 / sum(Ms2)
    sum_Ms2 <- sum_Ms2 + Ms2
    Ms[t, ] <- rowSums(Ms2)
  }

  list(Ms = Ms, LL = sum(LL), sum_Ms2 = sum_Ms2, Lp = Lp)
}

# Viterbi decode of the most-likely regime path (probability-space max-product).
Viterbi <- function(Pi, Lp, Z) {
  T <- nrow(Lp)
  M <- ncol(Lp)
  delta <- matrix(0, T, M)
  psi <- matrix(0, T, M)
  delta[1, ] <- Pi * Lp[1, ]
  for (t in 2:T) {
    for (j in 1:M) {
      max_val <- -Inf; max_state <- 0
      for (i in 1:M) {
        val <- delta[t - 1, i] * Z[i, j] * Lp[t, j]
        if (is.na(val) || is.nan(val)) val <- -Inf
        if (val > max_val) { max_val <- val; max_state <- i }
      }
      delta[t, j] <- max_val
      psi[t, j] <- max_state
    }
  }
  states <- numeric(T)
  states[T] <- which.max(delta[T, ])
  for (t in (T - 1):1) states[t] <- psi[t + 1, states[t + 1]]
  states
}

# Single-series MSVAR EM. Returns the best-LL iterate's parameters, smoothed
# regime posteriors Ms, Viterbi path x, and log-likelihood LL.
fit_msvar <- function(y, M, p = 1, theta = NULL, control = list()) {
  if (!is.matrix(y)) dim(y) <- c(length(y), 1)
  M <- as.integer(M); stopifnot(M > 0)
  p <- as.integer(p); stopifnot(p > 0); stopifnot(p < nrow(y))

  con <- list(maxit = 1000L, tol = 1e-6, verbose = FALSE, fixed = NULL, init = "mvn")
  idx <- intersect(names(control), names(con))
  if (length(idx) > 0) con[idx] <- control[idx]
  maxit <- con$maxit; tol <- con$tol; verbose <- con$verbose
  fixed <- con$fixed; init <- con$init
  choices <- c("mvn", "conditional")
  idx <- pmatch(init, choices)
  if (is.na(idx)) stop("Please specify 'init' as 'mvn' or 'conditional'.")
  init <- choices[idx]

  tt <- nrow(y); n <- ncol(y)
  LL_best <- LL <- -Inf

  X <- matrix(0, tt - p, p * n)
  for (k in 1:p) X[, ((k - 1) * n + 1):(k * n)] <- y[(p + 1 - k):(tt - k), ]
  sum_MP <- array(0, c(n, n, M))
  sum_MPb <- array(0, c(n * p, n * p, M))
  sum_MCP <- array(0, c(n, n * p, M))

  parnames <- c("mu", "A", "Q", "Pi", "Z")
  if (is.null(theta) && !(all(parnames %in% names(fixed)))) {
    args <- c("minlen", "maxlen", "nseg", "minclustsize", "cluster", "rmv_outliers")
    args <- intersect(args, names(control))
    theta <- if (length(args) == 0) {
      init_msvar(y, M, p, rmv_outliers = FALSE)
    } else {
      do.call(init_msvar, c(list(y = y, M = M, p = p), control[args]))
    }
  }
  parnames <- intersect(parnames, names(fixed))
  if (length(parnames) > 0) theta[parnames] <- fixed[parnames]
  theta_best <- theta
  Ms_best <- NULL
  Lp <- NULL

  LLflag <- 0L
  for (it in 1:maxit) {
    LL_prev <- LL
    filtsmooth <- skfs_var(y, M, p, theta, init)
    Ms <- filtsmooth$Ms
    LL <- filtsmooth$LL
    sum_Ms2 <- filtsmooth$sum_Ms2
    if (verbose) cat("Iteration", it, "LL = ", LL, "\n")

    if (LL > LL_best) {
      LL_best <- LL
      theta_best <- theta
      Ms_best <- Ms
      Lp <- filtsmooth$Lp
    }

    progress <- (LL - LL_prev) >= tol * abs(LL_prev)
    LLflag <- ifelse(progress, 0, LLflag + 1)
    if (LLflag == 5) break

    sum_Ms <- colSums(Ms[(p + 1):tt, ])
    xbar <- crossprod(X, Ms[(p + 1):tt, ]) / matrix(sum_Ms, n * p, M, byrow = TRUE)
    ybar <- crossprod(y[(p + 1):tt, ], Ms[(p + 1):tt, ]) / matrix(sum_Ms, n, M, byrow = TRUE)
    for (j in 1:M) {
      sqrtwj <- sqrt(Ms[(p + 1):tt, j])
      Xj <- X * sqrtwj
      yj <- y[(p + 1):tt, ] * sqrtwj
      sum_MP[, , j] <- crossprod(yj) - sum_Ms[j] * tcrossprod(ybar[, j])
      sum_MPb[, , j] <- crossprod(Xj) - sum_Ms[j] * tcrossprod(xbar[, j])
      sum_MCP[, , j] <- crossprod(yj, Xj) - sum_Ms[j] * tcrossprod(ybar[, j], xbar[, j])
    }
    theta <- M_var(Ms, sum_MCP, sum_MP, sum_MPb, sum_Ms2, xbar, ybar,
                   theta, fixed, init, verbose)
  }

  if (is.null(Ms_best)) Ms_best <- Ms
  if (is.null(Lp)) Lp <- skfs_var(y, M, p, theta_best, init)$Lp
  dim(theta_best$A) <- c(n, n, p, M)
  x <- Viterbi(theta_best$Pi, Lp, theta_best$Z)
  out <- list(theta = theta_best, x = x, Ms = Ms_best, LL = LL_best, iters = it)
  attr(out, "model") <- "msvar"
  attr(out, "length") <- tt
  out
}

# Multi-series MSVAR EM with a shared model across trials (sequential E-step).
fit_msvar_multi <- function(y, M, p, theta, control = list(),
                            center = TRUE, scale = TRUE) {
  stopifnot(is.numeric(y) || is.list(y))
  if (is.list(y)) {
    test <- sapply(y, is.vector) | sapply(y, is.matrix)
    if (!all(test)) stop("All components of list 'y' must be vectors/matrices.")
    test <- sapply(y, NCOL)
    if (!all(test == test[1])) stop("All components of 'y' must share NCOL.")
  }
  M <- as.integer(M); p <- as.integer(p)
  stopifnot(M > 0); stopifnot(p > 0)

  if (is.list(y) && length(y) == 1) return(fit_msvar(y[[1]], M, p, theta, control))
  if (is.numeric(y) && NCOL(y) == 1) return(fit_msvar(y, M, p, theta, control))

  if (is.numeric(y)) {
    ndims <- length(dim(y))
    y <- apply(y, ndims, as.matrix, simplify = FALSE)
  } else {
    test <- sapply(y, is.matrix)
    if (!all(test)) y[!test] <- lapply(y[!test], as.matrix)
  }

  nsubs <- length(y)
  nvars <- NCOL(y[[1]])
  nvols <- sapply(y, NROW)

  if (center) for (i in 1:nsubs) y[[i]] <- sweep(y[[i]], 2, colMeans(y[[i]]), "-")
  if (scale)  for (i in 1:nsubs) y[[i]] <- sweep(y[[i]], 2, sqrt(colMeans(y[[i]]^2)), "/")

  con <- list(maxit = 1000L, tol = 1e-6, verbose = FALSE, fixed = NULL, init = "mvn")
  idx <- intersect(names(control), names(con))
  if (length(idx) > 0) con[idx] <- control[idx]
  maxit <- con$maxit; tol <- con$tol; verbose <- con$verbose
  fixed <- con$fixed; init <- con$init
  choices <- c("mvn", "conditional")
  idx <- pmatch(init, choices)
  if (is.na(idx)) stop("Please specify 'init' as 'mvn' or 'conditional'.")
  init <- choices[idx]

  LLbest <- LL <- -Inf
  LLsub <- numeric(nsubs)
  singlePiZ <- (length(theta$Pi) == M)
  t0 <- switch(init, mvn = 1, conditional = p + 1)

  parnames <- intersect(c("A", "Q", "mu", "Pi", "Z"), names(fixed))
  if (length(parnames) > 0) theta[parnames] <- fixed[parnames]
  thetabest <- theta
  Msbest <- NULL; LLsubbest <- LLsub; Lp <- NULL

  LLflag <- 0L
  for (it in 1:maxit) {
    LLprev <- LL
    theta_prev <- theta

    filtsmooth <- vector("list", nsubs)
    for (i in 1:nsubs) {
      thetai <- theta
      if (!singlePiZ) { thetai$Pi <- thetai$Pi[, i]; thetai$Z <- thetai$Z[, , i] }
      filtsmooth[[i]] <- skfs_var(y[[i]], M, p, thetai, init)
    }
    LLsub <- sapply(filtsmooth, "[[", "LL")
    LL <- sum(LLsub)
    Ms <- lapply(filtsmooth, "[[", "Ms")
    if (verbose) cat("\nIteration", it, "LL = ", LL)

    if (LL > LLbest) {
      LLsubbest <- LLsub; LLbest <- LL; thetabest <- theta
      Msbest <- Ms; Lp <- lapply(filtsmooth, "[[", "Lp")
    }

    progress <- (LL - LLprev) >= tol * abs(LLprev)
    LLflag <- ifelse(progress, 0, LLflag + 1)
    if (LLflag == 5 || it == maxit) break

    sum_Ms <- numeric(M)
    sum_MP <- array(0, c(nvars, nvars, M))
    sum_MPb <- array(0, c(nvars * p, nvars * p, M))
    sum_MCP <- array(0, c(nvars, nvars * p, M))
    if (singlePiZ) { Pi <- numeric(M); Z <- matrix(0, M, M) }
    else { Pi <- matrix(0, M, nsubs); Z <- array(0, c(M, M, nsubs)) }
    xbar <- matrix(0, nvars * p, M)
    ybar <- matrix(0, nvars, M)

    for (i in 1:nsubs) {
      Msi <- Ms[[i]]
      sum_Ms2i <- filtsmooth[[i]]$sum_Ms2
      if (singlePiZ) { Pi <- Pi + Msi[t0, ]; Z <- Z + sum_Ms2i }
      else { Pi[, i] <- Msi[t0, ]; Z[, , i] <- sum_Ms2i / rowSums(sum_Ms2i) }

      tt <- nrow(y[[i]])
      X <- matrix(0, tt - p, p * nvars)
      for (l in 1:p) X[, ((l - 1) * nvars + 1):(l * nvars)] <- y[[i]][(p + 1 - l):(tt - l), ]
      Msi <- Msi[-(1:p), ]
      sum_Ms <- sum_Ms + colSums(Msi)
      xbar <- xbar + crossprod(X, Msi)
      ybar <- ybar + crossprod(y[[i]][-(1:p), ], Msi)
      for (j in 1:M) {
        yj <- y[[i]][-(1:p), ] * sqrt(Msi[, j])
        sum_MP[, , j] <- sum_MP[, , j] + crossprod(yj)
        Xj <- X * sqrt(Msi[, j])
        sum_MPb[, , j] <- sum_MPb[, , j] + crossprod(Xj)
        sum_MCP[, , j] <- sum_MCP[, , j] + crossprod(yj, Xj)
      }
    }
    for (j in 1:M) {
      xbar[, j] <- xbar[, j] / sum_Ms[j]
      ybar[, j] <- ybar[, j] / sum_Ms[j]
      sum_MP[, , j] <- sum_MP[, , j] - sum_Ms[j] * tcrossprod(ybar[, j])
      sum_MPb[, , j] <- sum_MPb[, , j] - sum_Ms[j] * tcrossprod(xbar[, j])
      sum_MCP[, , j] <- sum_MCP[, , j] - sum_Ms[j] * tcrossprod(ybar[, j], xbar[, j])
    }

    Ms <- matrix(0, p + 2, M)
    Ms[p + 1, ] <- sum_Ms
    sum_Ms2 <- matrix(0, M, M)
    theta <- M_var(Ms, sum_MCP, sum_MP, sum_MPb, sum_Ms2, xbar, ybar,
                   theta_prev, fixed, init, verbose)

    if (singlePiZ) { Pi <- Pi / sum(Pi); Z <- Z / rowSums(Z) }
    Z[is.nan(Z)] <- 1 / M
    theta$Pi <- Pi
    theta$Z <- Z
  }

  if (is.null(Msbest)) Msbest <- Ms
  x <- vector("list", nsubs)
  if (!is.null(Lp)) {
    for (i in 1:nsubs) {
      if (singlePiZ) { Pi <- thetabest$Pi; Z <- thetabest$Z }
      else { Pi <- thetabest$Pi[, i]; Z <- thetabest$Z[, , i] }
      x[[i]] <- Viterbi(Pi, Lp[[i]], Z)
    }
  }
  out <- list(theta = thetabest, x = x, Ms = Msbest, LL = LLbest,
              LLsub = LLsubbest, iters = it)
  attr(out, "model") <- "msvar"
  attr(out, "length") <- nvols
  out
}

# Per-class trainer: per-trial MSVAR fits -> k-means cluster the AR matrices to
# seed a shared theta0 -> fit one group MSVAR across the class's trials.
.train_msvar <- function(train_data, M, p, control, update_fit0 = TRUE) {
  train_data <- lapply(train_data, as.matrix)
  n_channels <- ncol(train_data[[1]])
  ni <- length(train_data)

  fit0 <- lapply(train_data, fit_msvar, M = M, p = p, control = control)
  A_tmp <- unlist(lapply(fit0, function(z) z$theta$A))
  dim(A_tmp) <- c(n_channels, n_channels, p, M * ni)
  Q_tmp <- unlist(lapply(fit0, function(z) z$theta$Q))
  dim(Q_tmp) <- c(n_channels, n_channels, M * ni)
  mu_tmp <- unlist(lapply(fit0, function(z) z$theta$mu))
  dim(mu_tmp) <- c(n_channels, M * ni)

  A_reshaped <- matrix(A_tmp, nrow = n_channels * n_channels * p, ncol = M * ni)
  km <- kmeans(t(A_reshaped), centers = M)

  A_clustered <- array(0, c(n_channels, n_channels, p, M))
  Q_clustered <- array(0, c(n_channels, n_channels, M))
  mu_clustered <- matrix(0, n_channels, M)
  for (m in 1:M) {
    ci <- which(km$cluster == m)
    A_clustered[, , , m] <- apply(A_tmp[, , , ci, drop = FALSE], c(1, 2, 3), mean)
    Q_clustered[, , m] <- apply(Q_tmp[, , ci, drop = FALSE], c(1, 2), mean)
    mu_clustered[, m] <- rowMeans(mu_tmp[, ci, drop = FALSE])
  }

  if (update_fit0) {
    con2 <- control
    con2$fixed <- list(A = A_clustered, Q = Q_clustered, mu = mu_clustered)
    fit0b <- lapply(train_data, fit_msvar, M = M, p = p, control = con2)
    Pi_tmp <- rowMeans(sapply(fit0b, function(z) z$theta$Pi))
    Z_tmp <- matrix(apply(sapply(fit0b, function(z) as.vector(z$theta$Z)), 1, mean), M, M)
  } else {
    Pi_tmp <- rep(1 / M, M)
    Z_tmp <- matrix(1 / M, M, M)
  }

  theta0 <- list(A = A_clustered, Q = Q_clustered, mu = mu_clustered,
                 Pi = Pi_tmp, Z = Z_tmp)
  fit_msvar_multi(train_data, M, p, theta0, control = control)
}
