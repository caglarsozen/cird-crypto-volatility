# Install the packages used by cird_analysis.R.
packages <- c("quantmod", "xts", "zoo", "dplyr", "tidyr", "purrr",
              "ggplot2", "readr", "tibble", "curl", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
cat("All required R packages are available.\n")
