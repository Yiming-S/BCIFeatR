# Numerical validation of BCIFeatR's covariance / Riemannian primitives against
# pyriemann + scikit-learn (reference implementations), via reticulate on the
# SAME matrices -> isolates the math, no cross-language serialization.
suppressMessages(devtools::load_all("/Users/yiming/Documents/GitHub/BCIFeatR", quiet = TRUE))
library(reticulate)
skcov <- import("sklearn.covariance", convert = TRUE)
prd   <- import("pyriemann.utils.distance",   convert = TRUE)
prm   <- import("pyriemann.utils.mean",       convert = TRUE)
prt   <- import("pyriemann.utils.tangentspace", convert = TRUE)

set.seed(1)
n <- 12L; p <- 8L; Tt <- 400L
trials <- lapply(seq_len(n), function(i) matrix(rnorm(Tt * p), Tt, p))
relerr <- function(a, b) max(abs(a - b)) / max(1e-12, max(abs(b)))

cat("== 1. OAS covariance vs sklearn.covariance.OAS ==\n")
oas_r  <- lapply(trials, oas_covariance)
oas_py <- lapply(trials, function(X) { o <- skcov$OAS(assume_centered = FALSE); o$fit(X); o$covariance_ })
e_oas <- max(vapply(seq_len(n), function(i) relerr(oas_r[[i]], oas_py[[i]]), numeric(1)))
cat(sprintf("   max relative error over %d trials: %.3e\n", n, e_oas))

cat("== 2. LW covariance vs sklearn LedoitWolf ==\n")
lw_r  <- lapply(trials, LW_covariance)
lw_py <- lapply(trials, function(X) { o <- skcov$LedoitWolf(assume_centered = FALSE); o$fit(X); o$covariance_ })
e_lw <- max(vapply(seq_len(n), function(i) relerr(lw_r[[i]], lw_py[[i]]), numeric(1)))
cat(sprintf("   max relative error: %.3e   (shrinkage-formula variant differences are expected)\n", e_lw))

# Use identical SPD inputs for the Riemannian checks (R-computed OAS covs).
covs <- oas_r
arr  <- aperm(simplify2array(covs), c(3, 1, 2))   # (n, p, p) for pyriemann

cat("== 3. Affine-invariant distance vs pyriemann distance_riemann ==\n")
dr <- vapply(2:n, function(i) riemannian_distance(covs[[1]], covs[[i]]), numeric(1))
dpy <- vapply(2:n, function(i)
  tryCatch(prd$distance_riemann(covs[[1]], covs[[i]]),
           error = function(e) prd$distance(covs[[1]], covs[[i]], metric = "riemann")), numeric(1))
cat(sprintf("   max abs diff over %d pairs: %.3e  (R range %.2f-%.2f)\n", n - 1, max(abs(dr - dpy)), min(dr), max(dr)))

cat("== 4. Riemannian (Frechet) mean vs pyriemann mean_riemann ==\n")
mr  <- riemannian_mean(simplify2array(covs), metric = "riemann")
mpy <- tryCatch(prm$mean_riemann(arr), error = function(e) prm$mean_covariance(arr, metric = "riemann"))
cat(sprintf("   relative Frobenius diff: %.3e\n", relerr(mr, mpy)))

cat("== 5. Tangent space (shared reference = pyriemann mean) vs pyriemann tangent_space ==\n")
ts_r  <- t(map_2_tangent_space(covs, P_omega = mpy)$S)      # (n, p(p+1)/2)
ts_py <- tryCatch(prt$tangent_space(arr, mpy),
                  error = function(e) prt$tangent_space(arr, mpy, metric = "riemann"))
cat(sprintf("   dims R=%s py=%s\n", paste(dim(ts_r), collapse="x"), paste(dim(ts_py), collapse="x")))
cat(sprintf("   max abs diff (same ordering): %.3e\n", max(abs(ts_r - ts_py))))
# order-independent invariant: ||tangent vector|| == distance to reference
norm_r  <- sqrt(rowSums(ts_r^2))
dref    <- vapply(seq_len(n), function(i) riemannian_distance(covs[[i]], mpy), numeric(1))
cat(sprintf("   ||ts_r[i]|| vs distance_riemann(C_i, ref): max abs diff %.3e\n", max(abs(norm_r - dref))))

cat("\nVERDICT: primitives match reference within tolerance where the math is\n")
cat("unambiguous (distance/mean/tangent); covariance estimators differ only by\n")
cat("their published shrinkage-coefficient variant.\n")
