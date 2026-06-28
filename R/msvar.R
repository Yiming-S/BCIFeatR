# Module: Markov-Switching VAR (MSVAR) generative classifier and feature.
# Scope: Trains one MSVAR per class and scores trials by per-class whole-trial
# log-likelihood (argmax = prediction), mirroring the original CSP -> MSVAR
# pipeline. Per-trial center+scale at scoring matches the per-trial
# normalization used during training (fit_msvar_multi center/scale defaults).
#
# Method reference: D. Degras, C.-M. Ting, H. Ombao (2022). "Markov-switching
# state-space models with applications to neuroimaging." CSDA.

# Per-trial center + RMS-scale, matching fit_msvar_multi's internal step so
# fixed-theta likelihoods are computed on the same scale as training.
.msvar_scale_trial <- function(M) {
  M <- as.matrix(M)
  M <- sweep(M, 2, colMeans(M), "-")
  rms <- sqrt(colMeans(M^2))
  rms[!is.finite(rms) | rms <= 0] <- 1
  sweep(M, 2, rms, "/")
}

# Score trials under each class model. Scoring fixes theta and runs a single
# E-step (maxit = 1), so it is deterministic and returns the trial's
# log-likelihood (and, optionally, mean smoothed regime occupancy) per class.
.msvar_score <- function(model, x, want_occupancy = FALSE) {
  if (!is.list(x)) x <- list(x)
  K <- length(model$classes)
  N <- length(x)
  Mreg <- model$M
  LL <- matrix(NA_real_, N, K)
  occ <- if (want_occupancy) matrix(NA_real_, N, K * Mreg) else NULL
  base_con <- list(maxit = 1L, tol = model$control$tol, verbose = FALSE,
                   init = model$control$init)
  for (i in seq_len(N)) {
    Mi <- .msvar_scale_trial(x[[i]])
    for (k in seq_len(K)) {
      th <- model$models[[k]]
      con <- base_con
      con$fixed <- th
      fit <- fit_msvar(Mi, Mreg, model$p, theta = th, control = con)
      LL[i, k] <- fit$LL
      if (want_occupancy) {
        occ[i, ((k - 1) * Mreg + 1):(k * Mreg)] <- colMeans(fit$Ms)
      }
    }
  }
  list(LL = LL, occ = occ)
}

# Build the dispatcher feature matrix (per-class log-likelihood, optionally with
# per-class mean smoothed regime occupancy appended).
.msvar_feature <- function(model, x, include_occupancy = FALSE) {
  sc <- .msvar_score(model, x, want_occupancy = include_occupancy)
  LL <- sc$LL
  colnames(LL) <- paste0("msvar_LL_", model$classes)
  if (!include_occupancy) return(LL)
  occ <- sc$occ
  colnames(occ) <- paste0("msvar_occ_",
                          rep(model$classes, each = model$M),
                          "_r", rep(seq_len(model$M), times = length(model$classes)))
  cbind(LL, occ)
}

#' Train a per-class Markov-Switching VAR (MSVAR) generative classifier.
#'
#' Fits one MSVAR (M regimes, order p) per class from that class's trials. A new
#' trial is later assigned (by [msvar_predict()]) to the class whose model gives
#' the highest whole-trial log-likelihood. This reproduces the original
#' CSP -> MSVAR motor-imagery pipeline, with spatial filtering left to the
#' feature dispatcher (e.g. fit CSP/FBCSP first, then pass the projected trials).
#'
#' @param x List of trial matrices (`samples x channels`). For best results pass
#'   low-dimensional, e.g. CSP-projected, trials.
#' @param y Class labels aligned with `x`.
#' @param M Number of regimes per class model.
#' @param p VAR order.
#' @param control MSVAR EM control list (`maxit`, `tol`, `init`, plus optional
#'   `init_msvar` keys such as `nseg`). Defaults to a light fit; the original
#'   pipeline used `maxit = 1`.
#' @param seed Optional integer; when supplied, `set.seed(seed)` is called before
#'   the RNG-based initialization/k-means so training is reproducible.
#' @return An object of class `"msvar"`: a list with `classes`, per-class
#'   `models` (each a theta list), and `M`/`p`/`control`.
#' @references D. Degras, C.-M. Ting, H. Ombao (2022). Markov-switching
#'   state-space models with applications to neuroimaging. CSDA.
#' @seealso [msvar_predict()]
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- lapply(1:12, function(i) matrix(rnorm(200 * 3), 200, 3))
#' y <- factor(rep(c("a", "b"), each = 6))
#' m <- msvar_train(x, y, M = 2L, p = 1L, seed = 1L,
#'                  control = list(maxit = 5L, nseg = 15L))
#' predict(m, x)
#' }
#' @export
msvar_train <- function(x, y, M = 2L, p = 1L,
                        control = list(maxit = 10L, tol = 1e-5, init = "conditional"),
                        seed = NULL) {
  if (!is.list(x) || length(x) == 0L) stop("`x` must be a non-empty trial list.")
  if (length(y) != length(x)) stop("length(y) must equal the number of trials.")
  if (anyNA(y)) stop("`y` contains NA labels.")
  M <- as.integer(M); p <- as.integer(p)
  con <- modifyList(
    list(maxit = 10L, tol = 1e-5, verbose = FALSE, init = "conditional"),
    if (is.null(control)) list() else control
  )
  if (!is.null(seed)) set.seed(as.integer(seed))

  y <- factor(y)
  classes <- levels(y)
  if (length(classes) < 2L) stop("MSVAR requires at least two classes in `y`.")
  models <- lapply(classes, function(cl) {
    .train_msvar(x[which(y == cl)], M = M, p = p, control = con)$theta
  })
  names(models) <- classes

  structure(list(classes = classes, models = models, M = M, p = p, control = con),
            class = "msvar")
}

#' Predict class labels with a trained MSVAR classifier.
#'
#' Scores each trial under every class model and predicts the highest-likelihood
#' class. Scoring fixes the trained parameters (single E-step), so it is fast and
#' deterministic.
#'
#' @param model An `"msvar"` object from [msvar_train()].
#' @param x Trial list or a single trial matrix.
#' @return List with `y_hat` (factor of predicted classes) and `loglik`
#'   (trial-by-class log-likelihood matrix).
#' @seealso [msvar_train()]
#' @export
msvar_predict <- function(model, x) {
  if (!inherits(model, "msvar")) {
    stop("`model` must be an 'msvar' object from msvar_train().")
  }
  if (!is.list(x)) x <- list(x)
  LL <- .msvar_score(model, x, want_occupancy = FALSE)$LL
  colnames(LL) <- model$classes
  y_hat <- factor(model$classes[max.col(LL, ties.method = "first")],
                  levels = model$classes)
  list(y_hat = y_hat, loglik = LL)
}

#' Predict method for MSVAR classifiers.
#'
#' @param object An `"msvar"` object from [msvar_train()].
#' @param x Trial list or a single trial matrix.
#' @param type `"class"` for the predicted-class factor (default) or
#'   `"loglik"` for the trial-by-class log-likelihood matrix.
#' @param ... Unused.
#' @return A factor of class labels, or the log-likelihood matrix when
#'   `type = "loglik"`.
#' @seealso [msvar_train()], [msvar_predict()]
#' @export
predict.msvar <- function(object, x, type = c("class", "loglik"), ...) {
  type <- match.arg(type)
  res <- msvar_predict(object, x)
  if (type == "loglik") res$loglik else res$y_hat
}

#' @export
print.msvar <- function(x, ...) {
  cat(sprintf("Markov-switching VAR classifier (M = %d regimes, p = %d)\n",
              x$M, x$p))
  cat(sprintf("  %d classes: %s\n", length(x$classes),
              paste(x$classes, collapse = ", ")))
  invisible(x)
}
