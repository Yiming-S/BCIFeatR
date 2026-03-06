# Module: Feature extraction dispatchers and shared preprocessing/validation helpers.
# Scope: Converts trial/session inputs into feature matrices for all supported methods
# and preserves train-time objects required for deterministic test-time transforms.

#' Normalize scalar-like parameters from list wrappers.
#'
#' @param x A scalar value or a single-element list containing a scalar.
#' @return Scalar value used by downstream parameter readers.
.as_scalar_param <- function(x) {
  if (is.list(x)) x[[1]] else x
}

#' Apply second-stage CSP projection when provided.
#'
#' @param filtered_trials List of trial matrices after first-stage filtering.
#' @param filter2 Optional second-stage projection matrix.
#' @return List of projected trial matrices.
.apply_optional_filter2 <- function(filtered_trials, filter2) {
  if (is.null(filter2)) return(filtered_trials)
  lapply(filtered_trials, function(M) {
    if (ncol(M) == nrow(filter2)) {
      M %*% filter2
    } else if (ncol(M) == ncol(filter2)) {
      M %*% t(filter2)
    } else {
      stop("filter2 dimensions incompatible with trial matrix")
    }
  })
}

#' Vectorize upper-triangular entries from a covariance-like matrix list.
#'
#' @param mats List of square matrices.
#' @param mask Optional logical upper-triangular mask reused across train/test.
#' @return List with `features` matrix and the `mask` used for extraction.
.extract_upper_tri_features <- function(mats, mask = NULL) {
  stopifnot(is.list(mats), length(mats) > 0L)
  if (is.null(mask)) mask <- upper.tri(mats[[1]], diag = TRUE)
  feats <- t(vapply(mats, function(M) M[mask], numeric(sum(mask))))
  list(features = feats, mask = mask)
}

#' Resolve covariance estimator function from type string.
#'
#' @param cov_type One of `"oas"` or `"lw"`.
#' @return Covariance estimator function.
.resolve_cov_fun <- function(cov_type) {
  switch(tolower(cov_type),
         oas = oas_covariance,
         lw  = LW_covariance,
         stop("`cov_type` must be 'oas' or 'lw'."))
}

.is_positive_intish_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && (as.integer(x) == x)
}

.validate_frequency_bands <- function(frequency_bands, fs) {
  if (!is.list(frequency_bands) || length(frequency_bands) == 0L) {
    stop("`params$frequency_bands` must be a non-empty list of c(low, high) bands.")
  }
  for (i in seq_along(frequency_bands)) {
    band <- frequency_bands[[i]]
    if (!is.numeric(band) || length(band) != 2L || any(!is.finite(band))) {
      stop(sprintf("`params$frequency_bands[[%d]]` must be numeric length-2.", i))
    }
    lo <- band[1]
    hi <- band[2]
    if (lo <= 0 || hi <= 0 || lo >= hi) {
      stop(sprintf("`params$frequency_bands[[%d]]` must satisfy 0 < low < high.", i))
    }
    if (!is.null(fs) && hi >= fs / 2) {
      stop(sprintf("`params$frequency_bands[[%d]]$high` must be < fs/2.", i))
    }
  }
  frequency_bands
}

.validate_feature_params <- function(feature, params) {
  if (is.null(params)) params <- list()
  if (!is.list(params)) stop("`params` must be a list.")

  if (feature == "logvar_pca") {
    if (!.is_positive_intish_scalar(params$ncomps)) {
      stop("`params$ncomps` must be a positive integer for logvar_pca.")
    }
    params$ncomps <- as.integer(params$ncomps)
  } else if (feature == "CSP") {
    if (is.null(params$ncomps)) params$ncomps <- 6L
    if (!.is_positive_intish_scalar(params$ncomps)) {
      stop("`params$ncomps` must be a positive integer for CSP.")
    }
    params$ncomps <- as.integer(params$ncomps)
  } else if (feature %in% c("FBCSP", "FBCSSP")) {
    if (!is.numeric(params$fs) || length(params$fs) != 1L || !is.finite(params$fs) || params$fs <= 0) {
      stop("`params$fs` must be a positive scalar for FBCSP/FBCSSP.")
    }
    params$frequency_bands <- .validate_frequency_bands(params$frequency_bands, fs = params$fs)
    if (is.null(params$ncomps)) {
      params$ncomps <- if (feature == "FBCSSP") c(3L, 3L) else 5L
    }
    ncomps <- .as_scalar_param(params$ncomps)
    if (!is.numeric(ncomps) || any(!is.finite(ncomps)) || any(ncomps <= 0)) {
      stop("`params$ncomps` must contain positive values for FBCSP/FBCSSP.")
    }
    params$ncomps <- ncomps
    if (!is.null(params$channels)) {
      channels <- .as_scalar_param(params$channels)
      if (!is.numeric(channels) || any(!is.finite(channels)) ||
          any(channels < 1) || any(as.integer(channels) != channels)) {
        stop("`params$channels` must contain positive integer channel indices.")
      }
      params$channels <- as.integer(channels)
    }
  } else if (feature == "TS") {
    if (is.null(params$cov_type)) params$cov_type <- "oas"
    .resolve_cov_fun(params$cov_type)
  } else if (feature == "ACM_TS") {
    if (is.null(params$order)) params$order <- 2L
    if (is.null(params$delay)) params$delay <- 1L
    if (is.null(params$shrinkage)) params$shrinkage <- "oas"
    if (!.is_positive_intish_scalar(params$order)) {
      stop("`params$order` must be a positive integer for ACM_TS.")
    }
    if (!.is_positive_intish_scalar(params$delay)) {
      stop("`params$delay` must be a positive integer for ACM_TS.")
    }
    params$order <- as.integer(params$order)
    params$delay <- as.integer(params$delay)
    params$shrinkage <- match.arg(params$shrinkage, c("no", "LW", "oas"))
  } else if (feature == "Riemannian") {
    if (is.null(params$cov_type)) params$cov_type <- "oas"
    .resolve_cov_fun(params$cov_type)
    if (!is.null(params$jitter_sd) &&
        (!is.numeric(params$jitter_sd) || length(params$jitter_sd) != 1L ||
         !is.finite(params$jitter_sd) || params$jitter_sd < 0)) {
      stop("`params$jitter_sd` must be a non-negative scalar for Riemannian.")
    }
  }
  params
}

#' Add deterministic jitter to break exact ties in flattened SPD features.
#'
#' @param feats Numeric feature matrix.
#' @param jitter_sd Non-negative jitter scale.
#' @return Feature matrix with deterministic perturbation.
.add_deterministic_jitter <- function(feats, jitter_sd = 0) {
  if (is.null(jitter_sd) || !is.finite(jitter_sd) || jitter_sd <= 0) return(feats)
  # Build a deterministic index grid so the same input always yields the same perturbation.
  idx <- matrix(seq_len(length(feats)), nrow = nrow(feats), ncol = ncol(feats))
  centered <- idx - mean(idx)
  scale_den <- max(abs(centered))
  if (scale_den == 0) return(feats)
  feats + (centered / scale_den) * jitter_sd
}

#' Normalize preprocess configuration to a canonical method list.
#'
#' @param preprocess NULL, scalar method name, or list with `method`.
#' @return List with normalized `method`.
.normalize_preprocess_spec <- function(preprocess) {
  if (is.null(preprocess)) return(list(method = "none"))
  if (is.character(preprocess) && length(preprocess) == 1L) {
    method <- tolower(preprocess)
    if (method == "zscore") method <- "scale"
    if (!method %in% c("none", "center", "scale", "whiten")) {
      stop("`params$preprocess` must be one of: none, center, scale, whiten.")
    }
    return(list(method = method))
  }
  if (is.list(preprocess)) {
    method <- if (is.null(preprocess$method)) "none" else tolower(preprocess$method)
    if (method == "zscore") method <- "scale"
    if (!method %in% c("none", "center", "scale", "whiten")) {
      stop("`params$preprocess$method` must be one of: none, center, scale, whiten.")
    }
    return(list(method = method))
  }
  stop("`params$preprocess` must be NULL, character scalar, or list.")
}

#' Fit channel-wise preprocessing statistics on trial data.
#'
#' @param trials List of trial matrices.
#' @param preprocess Normalized preprocess specification.
#' @param epsilon Lower bound for whitening eigenvalues.
#' @return Preprocessor object consumed by `.apply_trial_preprocessor()`.
.fit_trial_preprocessor <- function(trials, preprocess, epsilon = 1e-6) {
  method <- preprocess$method
  if (method == "none") return(list(method = "none"))

  Xall <- do.call(rbind, lapply(trials, as.matrix))
  center <- colMeans(Xall)
  scalev <- apply(Xall, 2, stats::sd)
  scalev[!is.finite(scalev) | scalev <= 0] <- 1

  if (method == "center") {
    return(list(method = method, center = center))
  }
  if (method == "scale") {
    return(list(method = method, center = center, scale = scalev))
  }

  # Whitening = standardize, estimate covariance, then apply inverse square root.
  Z <- sweep(sweep(Xall, 2, center, "-"), 2, scalev, "/")
  S <- stats::cov(Z)
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  lam <- pmax(eig$values, epsilon)
  V <- eig$vectors
  whitener <- V %*% (t(V) / sqrt(lam))
  list(method = method, center = center, scale = scalev, whitener = whitener)
}

#' Apply a fitted preprocessor to each trial matrix.
#'
#' @param trials List of trial matrices.
#' @param preproc Object returned by `.fit_trial_preprocessor()`.
#' @return List of transformed trial matrices.
.apply_trial_preprocessor <- function(trials, preproc) {
  if (is.null(preproc) || is.null(preproc$method) || identical(preproc$method, "none")) {
    return(lapply(trials, as.matrix))
  }

  lapply(trials, function(tr) {
    M <- as.matrix(tr)
    if (identical(preproc$method, "center")) {
      return(sweep(M, 2, preproc$center, "-"))
    }
    if (identical(preproc$method, "scale")) {
      return(sweep(sweep(M, 2, preproc$center, "-"), 2, preproc$scale, "/"))
    }
    if (identical(preproc$method, "whiten")) {
      Z <- sweep(sweep(M, 2, preproc$center, "-"), 2, preproc$scale, "/")
      return(Z %*% preproc$whitener)
    }
    stop("Unknown preprocess method in object.")
  })
}

#' Validate a single trial matrix shape and numeric finiteness.
#'
#' @param trial Candidate trial matrix.
#' @param path Name used in error messages.
#' @return Invisibly returns NULL or raises an error.
.assert_trial_matrix <- function(trial, path) {
  if (!is.matrix(trial)) stop(sprintf("%s must be a matrix.", path))
  if (!is.numeric(trial)) stop(sprintf("%s must be numeric.", path))
  if (nrow(trial) < 2L) stop(sprintf("%s must have at least 2 rows.", path))
  if (ncol(trial) < 1L) stop(sprintf("%s must have at least 1 column.", path))
  if (any(!is.finite(trial))) stop(sprintf("%s contains non-finite values.", path))
}

#' Validate trial-list contract and channel consistency.
#'
#' @param trials List of trial matrices.
#' @param path Name used in error messages.
#' @return Invisibly returns NULL or raises an error.
.validate_trial_list <- function(trials, path = "x") {
  if (!is.list(trials) || length(trials) == 0L) {
    stop(sprintf("%s must be a non-empty trial list.", path))
  }
  for (i in seq_along(trials)) {
    .assert_trial_matrix(trials[[i]], sprintf("%s[[%d]]", path, i))
  }
  nch <- vapply(trials, ncol, integer(1))
  if (any(nch != nch[1])) {
    stop(sprintf("%s trials must have the same number of channels.", path))
  }
}

#' Validate label vector/list against trial count.
#'
#' @param labels Labels to validate.
#' @param n_trials Expected number of trials.
#' @param path Name used in error messages.
#' @param allow_null Whether NULL labels are accepted.
#' @return Invisibly returns NULL or raises an error.
.validate_labels <- function(labels, n_trials, path = "y", allow_null = FALSE) {
  if (allow_null && is.null(labels)) return(invisible(NULL))
  if (is.null(labels)) stop(sprintf("%s must not be NULL.", path))
  if (length(labels) != n_trials) stop(sprintf("length(%s) must equal number of trials.", path))
  if (any(is.na(labels))) stop(sprintf("%s contains NA labels.", path))
}

#' Validate training input contract for trial-list and session-list modes.
#'
#' @param x Trial list or session list.
#' @param y Label vector or label-list aligned with `x`.
#' @param feature Feature mode to enforce feature-specific constraints.
#' @return Invisibly returns NULL or raises an error.
.validate_train_contract <- function(x, y, feature) {
  if (!is.list(x) || length(x) == 0L) stop("`x` must be a non-empty list.")
  input_is_session_list <- is.list(x[[1]])

  if (input_is_session_list) {
    if (feature == "ATM" && is.null(y)) {
      for (i in seq_along(x)) {
        .validate_trial_list(x[[i]], sprintf("x[[%d]]", i))
      }
    } else {
      if (!is.list(y) || length(y) != length(x)) {
        stop("For session-list input, `y` must be a list with the same length as `x`.")
      }
      for (i in seq_along(x)) {
        .validate_trial_list(x[[i]], sprintf("x[[%d]]", i))
        .validate_labels(y[[i]], length(x[[i]]), sprintf("y[[%d]]", i), allow_null = (feature == "ATM"))
        if (feature %in% c("CSP", "FBCSP", "FBCSSP") &&
            !is.null(y[[i]]) && length(unique(y[[i]])) < 2L) {
          stop(sprintf("%s requires at least two classes in y[[%d]].", feature, i))
        }
      }
    }
  } else {
    .validate_trial_list(x, "x")
    .validate_labels(y, length(x), "y", allow_null = (feature == "ATM"))
    if (feature %in% c("CSP", "FBCSP", "FBCSSP") && !is.null(y) && length(unique(y)) < 2L) {
      stop(sprintf("%s requires at least two classes in y.", feature))
    }
  }
}

#' Validate test-time contract for selected feature mode.
#'
#' @param x Trial list.
#' @param object Train-time object produced by `featEx4Train`.
#' @param feature Feature mode name.
#' @return Invisibly returns NULL or raises an error.
.validate_test_contract <- function(x, object, feature) {
  .validate_trial_list(x, "x")
  if (feature != "logvar" && is.null(object)) {
    stop("`object` must be provided for this feature in featEx4Test().")
  }
  if (feature %in% c("FBCSP", "FBCSSP") && is.null(object$filter)) {
    stop("`object$filter` is required for FBCSP/FBCSSP test transform.")
  }
}

#' Build standardized metadata for train outputs.
#'
#' @param feature Feature mode name.
#' @param feats Feature matrix.
#' @param preprocess_method Preprocess method string.
#' @param params Raw parameter list used during training.
#' @return Metadata list stored in `object$metadata`.
.build_metadata <- function(feature, feats, preprocess_method, params = NULL) {
  list(
    feature = feature,
    n_trials = nrow(feats),
    n_features = ncol(feats),
    preprocess = preprocess_method,
    version = "0.2",
    params = params
  )
}


#' Train feature extractors and compute feature matrix.
#'
#' Supports both trial-list input (`x` is list of matrices) and session-list input
#' (`x` is list of trial lists). Returns train-time state required to reproduce the
#' same transform in `featEx4Test`.
#'
#' @param x Trial list or session list.
#' @param y Labels aligned with `x`; can be NULL only for ATM.
#' @param feature One of logvar/logvar_pca/CSP/FBCSP/FBCSSP/TS/ACM_TS/Riemannian/ATM.
#' @param params Feature-specific parameter list.
#' @param epsilon Numerical floor used in SPD operations.
#' @param simplify Whether to unwrap a single-session return to a single object.
#' @return List with `features` matrix and `object`; or list of such results for sessions.
#' @export
featEx4Train <- function(x, y, feature,
                         params = list(), epsilon = 1e-6, simplify = TRUE) {

  valid_features <- c("logvar", "logvar_pca", "CSP", "FBCSP", "FBCSSP",
                      "TS", "ACM_TS", "Riemannian", "ATM")
  feature <- match.arg(feature, valid_features)
  params <- .validate_feature_params(feature, params)
  .validate_train_contract(x, y, feature)

  input_is_session_list <- is.list(x[[1]])
  n_sess <- if (input_is_session_list) length(x) else 1L

  result <- vector("list", n_sess)
  for (i in seq_len(n_sess)) {
    # Unify input handling so each iteration works on a single session.
    if (input_is_session_list) {
      session_x <- x[[i]]
      session_y <- if (is.null(y)) NULL else y[[i]]
    } else {
      session_x <- x
      session_y <- y
    }

    preproc_spec <- .normalize_preprocess_spec(params$preprocess)
    preproc_obj  <- .fit_trial_preprocessor(session_x, preproc_spec, epsilon = epsilon)
    # All feature branches operate on the same preprocessed trials.
    session_x_proc <- .apply_trial_preprocessor(session_x, preproc_obj)

    if (feature == "logvar") {
      feats <- log_var(session_x_proc)
      obj <- list()

    } else if (feature == "logvar_pca") {
      ncomps <- params$ncomps
      Xall <- do.call(rbind, session_x_proc)
      pca_res <- logvar_transform(Xall, NULL, FALSE, ncomps)
      pca_mod <- pca_res$pca_model
      attr(pca_mod, "ncomps") <- ncomps
      trial_feats <- lapply(
        session_x_proc,
        function(tr) logvar_transform(tr, pca_mod, TRUE, ncomps)$logvar
      )
      feats <- do.call(rbind, trial_feats)
      obj <- list(pca_mod = pca_mod, ncomps = ncomps)

    } else if (feature == "CSP") {
      ncomps <- params$ncomps
      csp_obj <- multiclass_csp(session_x_proc, session_y, ncomps)
      proj <- lapply(session_x_proc, `%*%`, csp_obj$filter)
      feats <- log_var(proj)
      obj <- csp_obj

    } else if (feature %in% c("FBCSP", "FBCSSP")) {
      ncomps <- .as_scalar_param(params$ncomps)
      fs <- params$fs
      ch <- .as_scalar_param(params$channels)
      frequency_bands <- params$frequency_bands
      csp_info <- prepare_csp_filters(
        data_train = session_x_proc,
        labels_train = session_y,
        fs = fs,
        channels = ch,
        frequency_bands = frequency_bands,
        csp_type = feature,
        ncomps = ncomps
      )
      filt <- apply_csp_filters(session_x_proc, csp_info$filter)
      filt <- .apply_optional_filter2(filt, csp_info$filter2)
      feats <- extract_logvar_features(filt)
      obj <- csp_info

    } else if (feature == "TS") {
      cov_type_ts <- if (is.null(params$cov_type)) "oas" else tolower(params$cov_type)
      cov_fun <- .resolve_cov_fun(cov_type_ts)
      cov_list <- lapply(session_x_proc, function(tr) cov_fun(as.matrix(tr)))
      # Map SPD covariances to a Euclidean tangent space for linear models.
      tg <- map_2_tangent_space(cov_list, epsilon = epsilon)
      feats <- t(tg$S)
      obj <- list(tangent_ref = tg$P_omega, cov_type = cov_type_ts, varsel = NULL)

    } else if (feature == "ACM_TS") {
      # Build augmented covariance matrices from delayed channel embeddings.
      cov_list <- compute_ACM(
        session_x_proc,
        params$order,
        params$delay,
        params$shrinkage
      )
      extra_obj <- NULL
      if (isTRUE(params$use_filter)) {
        vs <- riem_var_select(cov_list, epsilon = epsilon)
        # Filter in tangent space, then map back to SPD before final vectorization.
        cov_list <- lapply(
          cov_list, geodesic_filtering,
          P_omega = vs$P_omega, W_tilde = vs$W_tilde, epsilon = epsilon
        )
        extra_obj <- vs
      }
      tg <- map_2_tangent_space(cov_list, epsilon = epsilon)
      feats <- t(tg$S)
      obj <- list(
        tangent_ref = tg$P_omega,
        varsel = extra_obj,
        order = params$order,
        delay = params$delay,
        shrinkage = params$shrinkage
      )

    } else if (feature == "Riemannian") {
      cov_type_riem <- if (is.null(params$cov_type)) "oas" else tolower(params$cov_type)
      cov_fun <- .resolve_cov_fun(cov_type_riem)
      cov_list <- lapply(session_x_proc, function(tr) cov_fun(as.matrix(tr)))
      jitter_sd <- if (is.null(params$jitter_sd)) 1e-8 else params$jitter_sd

      if (isTRUE(params$use_filter)) {
        vs <- riem_var_select(cov_list, epsilon = epsilon)
        # Geodesic filtering keeps the result on the SPD manifold.
        cov_list <- lapply(
          cov_list, geodesic_filtering,
          P_omega = vs$P_omega, W_tilde = vs$W_tilde, epsilon = epsilon
        )
        flat <- .extract_upper_tri_features(cov_list)
        feats <- flat$features
        obj <- list(
          varsel = vs,
          mask = flat$mask,
          cov_type = cov_type_riem,
          jitter_sd = 0,
          use_filter = TRUE
        )
      } else {
        flat <- .extract_upper_tri_features(cov_list)
        # Deterministic jitter avoids exact duplicate columns without RNG side effects.
        feats <- .add_deterministic_jitter(flat$features, jitter_sd = jitter_sd)
        obj <- list(
          mean_cov = riemannian_mean(simplify2array(cov_list)),
          mask = flat$mask,
          cov_type = cov_type_riem,
          jitter_sd = jitter_sd,
          use_filter = FALSE
        )
      }

    } else if (feature == "ATM") {
      # Note: ATM applies its own z-score normalization internally.
      # Setting params$preprocess may cause double normalization.
      atm_res <- atmTrain(x = session_x_proc, params = params)
      feats <- atm_res$features
      obj <- atm_res$object

    } else {
      stop("Unsupported feature.")
    }

    obj$preprocess <- preproc_obj
    obj$metadata <- .build_metadata(
      feature = feature,
      feats = feats,
      preprocess_method = preproc_obj$method,
      params = params
    )
    result[[i]] <- list(features = feats, object = obj)
  }

  if (n_sess == 1L && !input_is_session_list && isTRUE(simplify)) {
    result[[1]]
  } else {
    result
  }
}

#' Transform test trials using a trained feature object.
#'
#' @param x Trial list or single trial matrix.
#' @param object Train-time object returned by `featEx4Train`.
#' @param feature Feature mode, must match training branch.
#' @param epsilon Numerical floor used in SPD operations.
#' @return Feature matrix with one row per test trial.
#' @export
featEx4Test <- function(x, object,
                        feature = c("logvar", "logvar_pca",
                                    "CSP", "FBCSP", "FBCSSP",
                                    "TS", "ACM_TS", "Riemannian", "ATM"),
                        epsilon = 1e-6) {

  if (!is.list(x)) x <- list(x)
  feature <- match.arg(feature)
  .validate_test_contract(x, object, feature)

  preproc_obj <- if (is.list(object) && !is.null(object$preprocess)) object$preprocess else list(method = "none")
  x_proc <- .apply_trial_preprocessor(x, preproc_obj)

  if (feature == "logvar") {
    feats <- log_var(x_proc)

  } else if (feature == "logvar_pca") {
    if (is.list(object) && !inherits(object, "prcomp")) {
      pca_mod <- object$pca_mod
      ncomps <- object$ncomps
    } else {
      pca_mod <- object
      ncomps <- attr(pca_mod, "ncomps")
    }
    feats <- lapply(x_proc, logvar_transform,
                    pca_model = pca_mod, apply_pca = TRUE, n_components = ncomps)
    feats <- lapply(feats, "[[", "logvar")
    feats <- do.call(rbind, feats)

  } else if (feature == "CSP") {
    W <- if (!is.null(object$filter)) object$filter else object$W
    proj <- lapply(x_proc, `%*%`, W)
    feats <- log_var(proj)

  } else if (feature == "FBCSP" || feature == "FBCSSP") {
    sig <- apply_csp_filters(x_proc, object$filter)
    sig <- .apply_optional_filter2(sig, object$filter2)
    feats <- extract_logvar_features(sig)

  } else if (feature == "TS") {
    cov_fun <- .resolve_cov_fun(if (is.null(object$cov_type)) "oas" else object$cov_type)
    cov_ls <- lapply(x_proc, function(tr) cov_fun(as.matrix(tr)))
    # Reuse training tangent reference to preserve coordinate system.
    tg <- map_2_tangent_space(cov_ls, epsilon = epsilon, P_omega = object$tangent_ref)
    feats <- t(tg$S)

  } else if (feature == "ACM_TS") {
    cov_ls <- compute_ACM(x_proc, object$order, object$delay, object$shrinkage)
    if (!is.null(object$varsel)) {
      # Apply the same manifold filter learned during training.
      cov_ls <- lapply(
        cov_ls, geodesic_filtering,
        P_omega = object$varsel$P_omega,
        W_tilde = object$varsel$W_tilde,
        epsilon = epsilon
      )
    }
    tg <- map_2_tangent_space(cov_ls, epsilon = epsilon, P_omega = object$tangent_ref)
    feats <- t(tg$S)

  } else if (feature == "Riemannian") {
    cov_type_riem <- if (is.null(object$cov_type)) "oas" else tolower(object$cov_type)
    cov_fun <- .resolve_cov_fun(cov_type_riem)
    mask <- object$mask
    cov_ls <- lapply(x_proc, function(tr) cov_fun(as.matrix(tr)))
    jitter_sd <- if (is.null(object$jitter_sd)) 0 else object$jitter_sd
    if (!is.null(object$varsel)) {
      cov_ls <- lapply(
        cov_ls, geodesic_filtering,
        P_omega = object$varsel$P_omega,
        W_tilde = object$varsel$W_tilde,
        epsilon = epsilon
      )
    }
    flat <- .extract_upper_tri_features(cov_ls, mask = mask)
    # Keep deterministic tie-breaking behavior aligned with training.
    feats <- .add_deterministic_jitter(flat$features, jitter_sd = jitter_sd)

  } else if (feature == "ATM") {
    feats <- atmTest(x = x_proc, object = object)
  }

  if (is.list(object) && is.list(object$metadata) && !is.null(object$metadata$n_features)) {
    if (ncol(feats) != object$metadata$n_features) {
      stop("Test feature dimension does not match training metadata.")
    }
  }
  feats
}
