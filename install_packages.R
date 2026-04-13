# One-shot install of CRAN dependencies for this repository.
# Usage (from repository root):  Rscript install_packages.R

pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "lubridate",
  "rdrobust",
  "truncnorm",
  "foreach",
  "doSNOW",
  "haven",
  "purrr",
  "plotly",
  "htmlwidgets",
  "fdapace",
  "pracma",
  "manifold",
  "Matrix",
  "osqp",
  "trust"
)

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

message("All packages available: ", paste(pkgs, collapse = ", "))
