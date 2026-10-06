# Compare a regenerated analysis against the supplied numerical reference.
# Usage: Rscript verify_results.R results reported
# PNG byte equality is not expected across operating systems or font devices.
# All eight plot-data tables and PNG dimensions are checked instead.
args <- commandArgs(trailingOnly = TRUE)
generated <- if (length(args) >= 1L) args[1L] else "results"
reference <- if (length(args) >= 2L) args[2L] else "reported"
stopifnot(dir.exists(generated), dir.exists(reference))
files <- list.files(file.path(reference, "tables"), pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No reference CSV tables found.")
checks <- list()
stable_keys <- c("symbol", "method", "ablation_id", "weight_id", "dgp", "rep",
                 "sample", "episode_id", "sensitivity", "H", "B", "cal_frac",
                 "alpha_day", "h", "block_id")
normalize_table <- function(x) {
  keys <- intersect(stable_keys, names(x))
  if (length(keys)) {
    x <- x[do.call(order, c(unname(x[keys]), list(na.last = TRUE))), , drop = FALSE]
  }
  rownames(x) <- NULL
  for (j in seq_along(x)) {
    if (is.character(x[[j]])) x[[j]][!is.na(x[[j]]) & !nzchar(x[[j]])] <- NA_character_
  }
  x
}
for (path in files) {
  name <- basename(path)
  other <- file.path(generated, "tables", name)
  if (!file.exists(other)) stop("Missing regenerated table: ", other)
  a <- normalize_table(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  b <- normalize_table(read.csv(other, check.names = FALSE, stringsAsFactors = FALSE))
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) stop("Table structure differs: ", name)
  for (j in names(a)) {
    same <- isTRUE(all.equal(a[[j]], b[[j]], tolerance = 1e-9, check.attributes = FALSE))
    if (!same) stop("Numerical/reference mismatch: ", name, " / ", j)
  }
  checks[[length(checks) + 1L]] <- data.frame(check = paste0("table:", name), passed = TRUE)
}
figure_data <- c("cird_reliability_grid.csv", "cird_component_ablation_summary.csv",
                 "simulation_main_summary.csv", "benchmark_comparison.csv",
                 "cird_weight_sensitivity.csv", "cird_H_B_calibration_sensitivity.csv",
                 "simulation_parameter_sensitivity.csv", "cird_component_ablation.csv")
sizes <- list(c(6000L, 3720L), c(6000L, 2880L), c(6000L, 2880L), c(6000L, 3720L),
              c(6000L, 3720L), c(6000L, 3720L), c(5280L, 5640L), c(6000L, 3720L))
for (i in seq_len(8L)) {
  figure <- file.path(generated, "figures", paste0("Fig", i, ".png"))
  if (!file.exists(figure) || file.info(figure)$size < 1000) stop("Missing/empty figure: ", figure)
  con <- file(figure, "rb")
  magic <- readBin(con, "raw", n = 8L)
  if (!identical(as.integer(magic), c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))) stop("Invalid PNG: ", figure)
  seek(con, 16L)
  dimensions <- readBin(con, "integer", n = 2L, size = 4L, endian = "big")
  close(con)
  if (!identical(dimensions, sizes[[i]])) stop("Unexpected PNG dimensions: ", figure)
  if (!file.exists(file.path(reference, "tables", figure_data[i]))) stop("Missing plot-data reference: ", figure_data[i])
  checks[[length(checks) + 1L]] <- data.frame(check = paste0("figure:", i, ":data-and-dimensions"), passed = TRUE)
}
manifest <- read.csv(file.path("data", "checksums_md5.csv"), stringsAsFactors = FALSE)
actual <- unname(tools::md5sum(file.path("data", manifest$file)))
if (!identical(actual, manifest$md5)) stop("Frozen-input checksums differ.")
checks[[length(checks) + 1L]] <- data.frame(check = "frozen-input-checksums", passed = TRUE)
out <- do.call(rbind, checks)
dir.create(file.path(generated, "audit"), showWarnings = FALSE)
write.csv(out, file.path(generated, "audit", "clean_run_comparison.csv"), row.names = FALSE)
cat("PASS:", length(files), "numerical CSV tables, eight figure-data sets and PNG dimensions, and four frozen-input checksums.\n")
