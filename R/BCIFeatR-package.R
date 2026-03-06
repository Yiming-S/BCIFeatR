#' @keywords internal
"_PACKAGE"

#' @importFrom stats cov var sd prcomp p.adjust t.test cor
#' @importFrom utils combn head
#' @importFrom geigen geigen
#' @importFrom gsignal butter filtfilt pwelch
NULL

# Functions checked at runtime via exists() — not part of this package.
utils::globalVariables(c("perf_log_event", "splitTimeRange"))
