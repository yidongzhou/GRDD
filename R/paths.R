# Project root for file paths. Set environment variable GRDD_ROOT to the repository
# root if you cannot run scripts with working directory = GRDD (e.g. some IDEs).
grdd_root <- function() {
  r <- Sys.getenv("GRDD_ROOT", "")
  if (nzchar(r)) {
    return(normalizePath(r, winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

grdd_path <- function(...) {
  file.path(grdd_root(), ...)
}
