#' @keywords internal
"_PACKAGE"

#' @importFrom stats cov var sd prcomp p.adjust t.test cor median mad
#' @importFrom stats dist hclust cutree kmeans runif
#' @importFrom utils combn head modifyList
#' @importFrom gsignal butter filtfilt pwelch
NULL

# Optional perf hook resolved at runtime via exists() — not part of this package.
utils::globalVariables("perf_log_event")
