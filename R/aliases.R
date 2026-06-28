# Module: snake_case aliases for public API.
# Scope: Provides consistent snake_case entry points alongside the original
# mixed-case names. Each alias is a direct binding to its camelCase counterpart
# (not a re-typed wrapper), so argument defaults can never drift between the two.
# The legacy names remain fully supported so existing scripts keep working.

# Direct bindings require the targets to exist when this file is sourced, so pin
# the source order via @include (roxygen writes a matching Collate field).
#' @include FeatureExtraction.R Bandpass_Filtering.R covariance.R compute_ACM.R
NULL

#' @rdname featEx4Train
#' @export
feat_ex_train <- featEx4Train

#' @rdname featEx4Test
#' @export
feat_ex_test <- featEx4Test

#' @rdname LW_covariance
#' @export
lw_covariance <- LW_covariance

#' @rdname compute_ACM
#' @export
compute_acm <- compute_ACM

#' @rdname freqBank
#' @export
freq_bank <- freqBank

#' @rdname freqBandSelect
#' @export
freq_band_select <- freqBandSelect

#' @rdname splitTimeRange
#' @export
split_time_range <- splitTimeRange
