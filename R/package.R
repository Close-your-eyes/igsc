#' tcrdistR
#'
#' Calculate T-cell receptor distances using the Python tcrdist3
#' package.
#'
#' @keywords internal
"_PACKAGE"


.onLoad <- function(libname, pkgname) {
  reticulate::py_require(c(
    "pandas",
    "tcrdist3"
  ))
}
