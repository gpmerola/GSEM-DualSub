install.packages(c("GenomicSEM", "readr", "dplyr", "ggplot2", "qqman", "data.table", "utils", "reshape2", "corrplot", "forcats", "tidyr", "Matrix"))

# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("johanzvrskovec/shru")
