# List of packages required
packages <- c(
  "corrplot",
  "ggplot2",
  "dplyr",
  "Matrix",
  "shru",
  "GenomicSEM",
  "readr",
  "qqman",
  "data.table",
  "utils",
  "reshape2",
  "devtools"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

if (!require("GenomicSEM", character.only = TRUE)) {
  devtools::install_github("GenomicSEM/GenomicSEM")
}

if (!require("shru", character.only = TRUE)) {
  devtools::install_github("GenomicsMEL/shru")
}
