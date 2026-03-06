# Module: Augmented covariance matrix construction for time-delay embeddings.
# Scope: Expands each trial with delayed channel copies and returns covariance
# matrices with optional shrinkage estimators.

#' Compute augmented covariance matrices from delayed embeddings.
#'
#' @param x Trial matrix or list of trial matrices (`samples x channels`).
#' @param order Embedding order (number of delayed blocks).
#' @param delay Delay between blocks in samples.
#' @param shrinkage Covariance estimator: `no`, `LW`, or `oas`.
#' @return SPD covariance matrix or list of matrices matching input shape.
#' @export
compute_ACM <- function(x, order = 2, delay = 1, shrinkage = c("no", "LW", "oas")) {
	order <- as.integer(order)
	delay <- as.integer(delay)
	if (!is.finite(order) || length(order) != 1L || order < 1L) {
		stop("`order` must be an integer >= 1.")
	}
	if (!is.finite(delay) || length(delay) != 1L || delay < 1L) {
		stop("`delay` must be an integer >= 1.")
	}
	ntrials <- if (is.list(x)) length(x) else 1
	acm <- vector("list", ntrials)
	shrinkage <- match.arg(shrinkage)
	cov_fun <- switch(shrinkage, no = cov, LW = LW_covariance, oas = oas_covariance)
	
	for (i in 1:ntrials) {
		xi <- if (is.list(x)) x[[i]] else x
		if (order == 1) {
			augx <- xi
		} else {
			nsamp <- nrow(xi)
			nchan <- ncol(xi)
			eff_len <- nsamp - (order - 1L) * delay
			if (eff_len < 2L)
				stop(sprintf("Insufficient time points (%d) for order=%d, delay=%d", nsamp, order, delay))
			idx <- lapply(1:order, function(o) seq(1+(o-1)*delay, nsamp - (order-o)*delay))
			augx <- matrix(xi[unlist(idx),], eff_len, nchan * order)
		}
		acm[[i]] <- cov_fun(augx)
	}
	if (!is.list(x)) acm <- acm[[1]]
	return(acm)
}
