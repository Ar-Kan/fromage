load_fromage <- function(path = ".") {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required to load fromage from the source tree.", call. = FALSE)
  }

  pkgload::load_all(path, quiet = TRUE)
}
