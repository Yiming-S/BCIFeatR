# Module: snake_case aliases for public API.
# Scope: Provides consistent snake_case entry points alongside the original
# mixed-case names. The legacy names remain fully supported so existing
# downstream scripts keep working without modification.

#' @rdname featEx4Train
#' @export
feat_ex_train <- function(x, y, feature, params = list(), epsilon = 1e-6,
                          simplify = TRUE, skip_validation = FALSE) {
  featEx4Train(x = x, y = y, feature = feature, params = params,
               epsilon = epsilon, simplify = simplify,
               skip_validation = skip_validation)
}

#' @rdname featEx4Test
#' @export
feat_ex_test <- function(x, object, feature, epsilon = 1e-6,
                         skip_validation = FALSE) {
  featEx4Test(x = x, object = object, feature = feature, epsilon = epsilon,
              skip_validation = skip_validation)
}

#' @rdname LW_covariance
#' @export
lw_covariance <- function(x) LW_covariance(x)

#' @rdname compute_ACM
#' @export
compute_acm <- function(x, order = 2, delay = 1,
                        shrinkage = c("no", "LW", "oas")) {
  compute_ACM(x = x, order = order, delay = delay, shrinkage = shrinkage)
}

#' @rdname freqBank
#' @export
freq_bank <- function(eeg, fs, order = 3,
                      bands = list(c(8, 13), c(13, 18), c(18, 23), c(23, 28)),
                      filters = NULL) {
  freqBank(eeg = eeg, fs = fs, order = order, bands = bands, filters = filters)
}

#' @rdname freqBandSelect
#' @export
freq_band_select <- function(x, y, fs, flo = 1, fhi = fs / 2,
                             winlen = NULL, overlap = 0.5,
                             by_label = TRUE, trials = NULL) {
  freqBandSelect(x = x, y = y, fs = fs, flo = flo, fhi = fhi,
                 winlen = winlen, overlap = overlap,
                 by_label = by_label, trials = trials)
}
