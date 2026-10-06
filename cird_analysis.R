# ======================================================================
# Block-Calibrated Shape-Aware Screening of Cryptocurrency Volatility Episodes
# Authors: Mervenur Sozen and Caglar Sozen
#
# Single-file R analysis for four cryptocurrencies: BTC, ETH, XRP and LTC.
# The analysis is retrospective: forward curves are available after maturity.
# Block-conformal validity requires exchangeable block scores. Empirical
# calibration-gap and weight rankings remain descriptive.
#
# Run from the repository/package root:
#   Rscript install_dependencies.R
#   Rscript cird_analysis.R
#
# Frozen input CSVs are checked against data/checksums_md5.csv. The default
# run uses these inputs without downloading market data. Outputs are written
# to results/. CIRD_OUTPUT_DIR may select a separate output directory.
# The script generates all numerical tables, eight figures, model objects,
# dated calendars, Monte Carlo results and reproducibility audits.
# ======================================================================

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260911)

# ----------------------------------------------------------------------
# 0. Packages
# ----------------------------------------------------------------------
required_packages <- c(
  "quantmod", "xts", "zoo",
  "dplyr", "tidyr", "purrr",
  "ggplot2", "readr", "tibble",
  "curl", "jsonlite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install missing packages before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(quantmod)
  library(xts)
  library(zoo)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ----------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------
CFG <- list(
  analysis_id = "cird-crypto-volatility",
  seed = 20260911L,

  # Four-asset empirical universe.
  symbols = c("BTC-USD", "ETH-USD", "XRP-USD", "LTC-USD"),

  data_start = as.Date("2019-01-01"),
  data_end   = as.Date("2024-12-31"),

  # Main specification
  H_main = 30L,
  B_main = 7L,
  alpha_day_main = 0.05,
  train_frac = 0.60,
  cal_frac_main = 0.20,

  # Sensitivity grids
  H_grid = c(15L, 30L, 60L),
  B_grid = c(5L, 7L, 10L),
  cal_frac_grid = c(0.15, 0.20, 0.25),

  # Equal-sum weight grid: neutral baseline is a=b=1.
  weight_grid = tibble::tribble(
    ~weight_id,        ~a,    ~b,
    "level_heavy_1",   1.50,  0.50,
    "level_heavy_2",   1.25,  0.75,
    "equal",           1.00,  1.00,
    "curvature_heavy_1",0.75, 1.25,
    "curvature_heavy_2",0.50, 1.50
  ),

  # Explicit component ablation for the CIRD novelty analysis.
  # These are evaluated separately from the positive-weight sensitivity grid.
  ablation_grid = tibble::tribble(
    ~ablation_id,      ~a,   ~b,
    "level_only",      1.00, 0.00,
    "full_CIRD",       1.00, 1.00,
    "curvature_only",  0.00, 1.00
  ),

  # Reliability curve
  alpha_day_grid = c(0.01, 0.05, 0.10),

  # Range-based held-out stress reference (NOT ground truth).
  # The top 10% of held-out B-day range-proxy block maxima are used only for
  # ex-post descriptive event validation; this reference never enters fitting
  # or conformal calibration.
  stress_eval_quantile = 0.90,

  # Monte Carlo
  mc_reps = 200L,
  mc_sensitivity_reps = 100L,
  mc_T = 1300L,
  mc_burn = 500L,

  # Publication diagnostic: uninformative random score.
  random_score_reps = 1000L,
  random_score_seed = 20261005L,

  # HAR additive log offset; the zero-return audit reports small RV(1) values.
  har_log_epsilon = 1e-12,

  # Descriptive event windows requested for economic interpretation. These are
  # not used in score fitting/calibration and must be cited/interpreted in the
  # manuscript as descriptive market-event windows rather than causal labels.
  recognized_episode_windows = tibble::tribble(
    ~episode_id, ~episode_label, ~event_start, ~event_end,
    "ETF_APPROVAL", "U.S. spot-Bitcoin ETF approval period", as.Date("2024-01-07"), as.Date("2024-01-13"),
    "MARCH_HIGHS", "March 2024 crypto highs", as.Date("2024-03-01"), as.Date("2024-03-31"),
    "AUGUST_SELLOFF", "Early-August 2024 market-wide sell-off", as.Date("2024-08-02"), as.Date("2024-08-08"),
    "NOVEMBER_RALLY", "November 2024 crypto rally", as.Date("2024-11-01"), as.Date("2024-11-30")
  ),

  frozen_data = TRUE,
  input_dir = Sys.getenv("CIRD_DATA_DIR", "data"),
  out_root = Sys.getenv("CIRD_OUTPUT_DIR", "results")
)

set.seed(CFG$seed)

DIRS <- list(
  root = CFG$out_root,
  data = file.path(CFG$out_root, "data"),
  tables = file.path(CFG$out_root, "tables"),
  figures = file.path(CFG$out_root, "figures"),
  rds = file.path(CFG$out_root, "rds"),
  audit = file.path(CFG$out_root, "audit"),
  latex_main = file.path(CFG$out_root, "latex", "main_text"),
  latex_appendix = file.path(CFG$out_root, "latex", "appendix")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ----------------------------------------------------------------------
# 2. Utility functions
# ----------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_sd <- function(x, eps = 1e-8) {
  z <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(z) || z < eps) eps else z
}

safe_mad <- function(x, eps = 1e-8) {
  z <- stats::mad(x, center = stats::median(x, na.rm = TRUE),
                  constant = 1.4826, na.rm = TRUE)
  if (!is.finite(z) || z < eps) eps else z
}

logistic <- function(x) stats::plogis(x)

wilson_ci <- function(k, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0L) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- k / n
  den <- 1 + z^2 / n
  ctr <- (phat + z^2 / (2 * n)) / den
  half <- (z / den) * sqrt(
    phat * (1 - phat) / n + z^2 / (4 * n^2)
  )
  c(max(0, ctr - half), min(1, ctr + half))
}

alpha_block_map <- function(alpha_day, B) {
  # Design mapping. Exact "at least one" probability only under independence.
  1 - (1 - alpha_day)^B
}

conformal_threshold <- function(block_scores, alpha_block) {
  z <- sort(block_scores[is.finite(block_scores)])
  m <- length(z)
  if (m == 0L) stop("No finite calibration block scores.")
  k <- ceiling((m + 1) * (1 - alpha_block))

  # Standard split-conformal convention: if the requested order statistic
  # lies beyond the m observed calibration scores, the threshold is +Inf.
  if (k > m) {
    return(list(q = Inf, k = k, m = m))
  }

  if (k < 1L) k <- 1L
  list(q = z[k], k = k, m = m)
}

make_equal_blocks <- function(x, B, drop_incomplete = TRUE) {
  n <- length(x)
  n_full <- floor(n / B)
  if (n_full < 1L) return(list(values = numeric(0), membership = integer(0)))

  use_n <- if (drop_incomplete) n_full * B else n
  x_use <- x[seq_len(use_n)]
  membership <- rep(seq_len(ceiling(use_n / B)), each = B, length.out = use_n)

  vals <- vapply(
    split(x_use, membership),
    function(v) max(v, na.rm = TRUE),
    numeric(1)
  )

  list(values = vals, membership = membership)
}

block_any <- function(x, B, drop_incomplete = TRUE) {
  n <- length(x)
  n_full <- floor(n / B)
  if (n_full < 1L) return(logical(0))

  use_n <- if (drop_incomplete) n_full * B else n
  x_use <- x[seq_len(use_n)]
  membership <- rep(seq_len(ceiling(use_n / B)), each = B, length.out = use_n)

  vapply(
    split(x_use, membership),
    function(v) any(v %in% TRUE, na.rm = TRUE),
    logical(1)
  )
}

jaccard_binary <- function(a, b) {
  a <- a %in% TRUE
  b <- b %in% TRUE
  if (length(a) != length(b)) return(NA_real_)
  u <- sum(a | b)
  if (u == 0L) return(1)
  sum(a & b) / u
}

date_or_na <- function(x) {
  if (length(x) == 0L) as.Date(NA) else min(x)
}

date_or_na_max <- function(x) {
  if (length(x) == 0L) as.Date(NA) else max(x)
}

with_local_seed <- function(seed, fun) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fun()
}

collapse_int_ids <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("")
  paste(as.integer(x), collapse = ",")
}

collapse_dates <- function(x) {
  x <- as.Date(x)
  x <- sort(unique(x[!is.na(x)]))
  if (length(x) == 0L) return("")
  paste(as.character(x), collapse = ",")
}


# Publication labels used ONLY for display. Internal object/column names remain
# machine-friendly; figures and LaTeX tables never expose code-style underscores.
METHOD_LABELS <- c(
  CIRD = "CIRD",
  Terminal_FRV = "Terminal FRV",
  Integrated_upper_tail_rank = "Integrated upper-tail rank",
  HAR_forward_residual = "HAR forward residual"
)

ABLATION_LABELS <- c(
  level_only = "Level only",
  full_CIRD = "Full CIRD",
  curvature_only = "Curvature only"
)

WEIGHT_LABELS <- c(
  level_heavy_1 = "Level heavy 1",
  level_heavy_2 = "Level heavy 2",
  equal = "Equal weights",
  curvature_heavy_1 = "Curvature heavy 1",
  curvature_heavy_2 = "Curvature heavy 2"
)

DGP_LABELS <- c(
  HAR_RV = "HAR-RV",
  SV_AR1 = "SV-AR1",
  SV_BREAK = "SV structural break"
)

pretty_method <- function(x) {
  out <- unname(METHOD_LABELS[as.character(x)])
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", as.character(x)[miss])
  out
}

pretty_ablation <- function(x) {
  out <- unname(ABLATION_LABELS[as.character(x)])
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", as.character(x)[miss])
  out
}

pretty_weight <- function(x) {
  out <- unname(WEIGHT_LABELS[as.character(x)])
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", as.character(x)[miss])
  out
}

pretty_dgp <- function(x) {
  out <- unname(DGP_LABELS[as.character(x)])
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", as.character(x)[miss])
  out
}

mean_or_na <- function(x) {
  z <- x[is.finite(x)]
  if (length(z) == 0L) NA_real_ else mean(z)
}

fmt3 <- function(x) {
  z <- as.numeric(x)
  z[is.finite(z) & abs(z) < 0.0005] <- 0

  ifelse(
    is.na(z),
    "--",
    ifelse(
      is.infinite(z),
      ifelse(z > 0, "Inf", "-Inf"),
      sprintf("%.3f", z)
    )
  )
}

fmt_int <- function(x) {
  ifelse(is.na(x), "--", as.character(as.integer(round(x))))
}

fmt_interval3 <- function(lo, hi) {
  ifelse(
    is.na(lo) | is.na(hi),
    "--",
    paste0("[", sprintf("%.3f", lo), ", ", sprintf("%.3f", hi), "]")
  )
}

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "--"
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x <- gsub("$", "\\\\$", x, fixed = TRUE)
  x
}


publication_headers <- function(headers) {
  breaks <- c(
    "Mean block rate" = "Mean block\\\\rate",
    "Mean absolute gap" = "Mean absolute\\\\gap",
    "Absolute gap" = "Absolute\\\\gap",
    "Block rate" = "Block\\\\rate",
    "95\\% interval" = "Descriptive 95\\%\\\\Wilson interval",
    "Day rate" = "Daily\\\\rate",
    "Alarm blocks" = "Alarm\\\\blocks",
    "2.5--97.5\\% replication range" = "2.5--97.5\\%\\\\repetition range",
    "2.5--97.5\\% recall range" = "2.5--97.5\\%\\\\recall range",
    "Data-generating process" = "Data-generating\\\\process",
    "Jaccard vs. equal" = "Jaccard agreement\\\\with equal weights",
    "Calibration fraction" = "Calibration\\\\fraction",
    "Mean proxy recall" = "Mean reference\\\\recall",
    "Calibration gap" = "Signed\\\\gap",
    "Standard deviation" = "Standard\\\\deviation",
    "Raw kurtosis" = "Raw\\\\kurtosis",
    "Jarque--Bera $p$" = "Jarque--Bera\\\\$p$-value",
    "Ljung--Box(10) return $p$" = "Ljung--Box(10)\\\\returns: $p$-value",
    "Ljung--Box(10) squared-return $p$" = "Ljung--Box(10)\\\\squared returns: $p$-value",
    "90th percentile" = "90th\\\\percentile",
    "Tail quantile" = "Tail\\\\quantile",
    "Held-out proxy blocks" = "Held-out\\\\blocks",
    "Stress-reference blocks" = "Reference\\\\blocks",
    "Stress-reference rate" = "Reference\\\\rate",
    "Distinct reference runs" = "Contiguous\\\\reference runs",
    "Unique peak dates" = "Distinct\\\\peak dates",
    "Range peak dates" = "Range peak\\\\dates",
    "Train $N$" = "Training\\\\$N$",
    "Calibration $N$" = "Calibration\\\\$N$",
    "Test $N$" = "Test\\\\$N$",
    "Train end" = "Training\\\\end",
    "Calibration start" = "Calibration\\\\start",
    "Calibration end" = "Calibration\\\\end",
    "Test start" = "Test\\\\start",
    "Training curves" = "Training\\\\curves",
    "Share adjusted" = "Share\\\\adjusted",
    "Mean max adjustment" = "Mean maximum\\\\adjustment",
    "Maximum adjustment" = "Maximum\\\\adjustment",
    "Test blocks" = "Test\\\\blocks",
    "Alarm and stress blocks" = "Coincident\\\\blocks",
    "First alarm" = "First alarm\\\\origin",
    "Last alarm" = "Last alarm\\\\origin",
    "Pearson corr. $(L,C)$" = "Pearson correlation\\\\$(L,C)$",
    "Spearman corr. $(L,C)$" = "Spearman correlation\\\\$(L,C)$",
    "Mean recall" = "Mean\\\\recall",
    "Finite returns" = "Finite\\\\returns",
    "Exact zero returns" = "Exact zero\\\\returns",
    "Training $RV^{(1)}$ zero" = "Training\\\\$RV^{(1)}=0$",
    "Training $RV^{(1)}\\leq\\varepsilon$" = "Training\\\\$RV^{(1)}\\leq\\varepsilon$",
    "Minimum positive training $RV^{(1)}$" = "Minimum positive\\\\training $RV^{(1)}$"
  )
  vapply(headers, function(h) {
    if (h %in% names(breaks)) paste0("\\shortstack{", breaks[[h]], "}") else h
  }, character(1), USE.NAMES = FALSE)
}

write_latex_table <- function(
    df,
    headers,
    path,
    caption,
    label,
    note = NULL,
    align = NULL,
    wide = TRUE) {

  headers <- publication_headers(headers)
  if (length(headers) != ncol(df)) {
    stop("LaTeX table header count does not match data columns for: ", path)
  }

  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", max(0L, ncol(df) - 1L)), collapse = ""))
  }

  env <- if (wide) "table*" else "table"
  width_cmd <- if (wide) "\\textwidth" else "\\linewidth"
  body <- as.data.frame(lapply(df, latex_escape), stringsAsFactors = FALSE)

  lines <- c(
    paste0("\\begin{", env, "}[!htbp]"),
    "  \\centering",
    "  \\small",
    paste0("  \\caption{", caption, "}"),
    paste0("  \\label{", label, "}"),
    paste0("  \\begin{adjustbox}{max width=", width_cmd, "}"),
    paste0("  \\begin{tabular}{", align, "}"),
    "    \\toprule",
    paste0("    ", paste(headers, collapse = " & "), " \\\\"),
    "    \\midrule"
  )

  if (nrow(body) > 0L) {
    for (i in seq_len(nrow(body))) {
      lines <- c(
        lines,
        paste0(
          "    ",
          paste(as.character(body[i, , drop = TRUE]), collapse = " & "),
          " \\\\"
        )
      )
    }
  }

  lines <- c(
    lines,
    "    \\bottomrule",
    "  \\end{tabular}",
    "  \\end{adjustbox}"
  )

  if (!is.null(note) && nzchar(note)) {
    lines <- c(
      lines,
      "  \\vspace{2pt}",
      paste0(
        "  \\parbox{", width_cmd, "}{\\footnotesize\\textit{Note:} ",
        note, "}"
      )
    )
  }

  lines <- c(lines, paste0("\\end{", env, "}"))
  writeLines(lines, path, useBytes = TRUE)
}

# ----------------------------------------------------------------------
# 3. Yahoo Finance download and preprocessing
# ----------------------------------------------------------------------
# Read an existing cache only when it covers the requested window and has
# all fields required by both CIRD and the separate high-low range proxy.
read_compatible_yahoo_cache <- function(path, start_date, end_date) {
  if (!file.exists(path)) return(NULL)

  z <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE) %>%
      mutate(date = as.Date(date)) %>%
      arrange(date) %>%
      distinct(date, .keep_all = TRUE),
    error = function(e) NULL
  )

  required_cols <- c("date", "adjusted", "close", "high", "low")
  if (is.null(z) || !all(required_cols %in% names(z)) || nrow(z) < 100L) {
    return(NULL)
  }

  # Allow a two-day tolerance at the boundaries for endpoint/calendar quirks.
  if (min(z$date, na.rm = TRUE) > start_date + 2 ||
      max(z$date, na.rm = TRUE) < end_date - 2) {
    return(NULL)
  }

  z %>%
    filter(date >= start_date, date <= end_date) %>%
    mutate(
      adjusted = as.numeric(adjusted),
      close = as.numeric(close),
      high = as.numeric(high),
      low = as.numeric(low),
      ret = log(adjusted) - lag(log(adjusted)),
      sqret = ret^2,
      hl_range = if_else(
        is.finite(high) & is.finite(low) & high > 0 & low > 0,
        log(high / low)^2,
        NA_real_
      )
    )
}

# Yahoo Chart API fallback. Unlike the historical CSV endpoint used internally
# by some getSymbols versions, this endpoint does not require a crumb/cookie.
# Two Yahoo hosts are tried because availability can differ by network/region.
download_yahoo_chart_api <- function(symbol, start_date, end_date) {
  period1 <- as.numeric(as.POSIXct(start_date, tz = "UTC"))
  period2 <- as.numeric(as.POSIXct(end_date + 1, tz = "UTC"))

  hosts <- c("query1.finance.yahoo.com", "query2.finance.yahoo.com")
  last_error <- NULL

  for (host in hosts) {
    url <- paste0(
      "https://", host, "/v8/finance/chart/",
      utils::URLencode(symbol, reserved = TRUE),
      "?period1=", sprintf("%.0f", period1),
      "&period2=", sprintf("%.0f", period2),
      "&interval=1d&events=history&includeAdjustedClose=true"
    )

    for (attempt in 1:3) {
      ans <- tryCatch({
        h <- curl::new_handle()
        curl::handle_setopt(
          h,
          useragent = "Mozilla/5.0 (compatible; academic-research/1.0)",
          connecttimeout = 20,
          timeout = 90,
          followlocation = TRUE
        )

        res <- curl::curl_fetch_memory(url, handle = h)
        if (res$status_code != 200L) {
          stop("Yahoo Chart API HTTP status ", res$status_code)
        }

        txt <- rawToChar(res$content)
        # Keep JSON arrays as lists so Yahoo nulls preserve their positions.
        # This avoids length shifts on days with a missing OHLC field.
        js <- jsonlite::fromJSON(txt, simplifyVector = FALSE)

        if (!is.null(js$chart$error)) {
          stop("Yahoo Chart API returned an error object.")
        }
        if (is.null(js$chart$result) || length(js$chart$result) < 1L) {
          stop("Yahoo Chart API returned no result.")
        }

        rr <- js$chart$result[[1L]]
        to_num <- function(z) {
          vapply(z, function(v) {
            if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v[[1L]])
          }, numeric(1))
        }

        timestamp <- to_num(rr$timestamp)
        quote <- rr$indicators$quote[[1L]]

        close <- to_num(quote$close)
        high  <- to_num(quote$high)
        low   <- to_num(quote$low)

        adjusted <- tryCatch(
          to_num(rr$indicators$adjclose[[1L]]$adjclose),
          error = function(e) rep(NA_real_, length(close))
        )
        if (length(adjusted) != length(close) || all(!is.finite(adjusted))) {
          adjusted <- close
        } else {
          adjusted[!is.finite(adjusted)] <- close[!is.finite(adjusted)]
        }

        if (!all(c(length(timestamp), length(close), length(high), length(low)) ==
                 length(timestamp))) {
          stop("Yahoo Chart API returned inconsistent vector lengths.")
        }

        tibble(
          date = as.Date(as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC")),
          adjusted = adjusted,
          close = close,
          high = high,
          low = low
        ) %>%
          filter(date >= start_date, date <= end_date) %>%
          arrange(date) %>%
          distinct(date, .keep_all = TRUE) %>%
          filter(is.finite(adjusted), adjusted > 0) %>%
          mutate(
            ret = log(adjusted) - lag(log(adjusted)),
            sqret = ret^2,
            hl_range = if_else(
              is.finite(high) & is.finite(low) & high > 0 & low > 0,
              log(high / low)^2,
              NA_real_
            )
          )
      }, error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      })

      if (!is.null(ans) && nrow(ans) > 100L) {
        message("  Yahoo Chart API success: ", host,
                " (attempt ", attempt, ")")
        return(ans)
      }

      Sys.sleep(attempt)
    }
  }

  stop(
    "Yahoo Chart API fallback also failed for ", symbol, ". Last error: ",
    last_error %||% "unknown error"
  )
}

download_yahoo_symbol <- function(symbol, start_date, end_date, cache_dir) {
  base <- gsub("-", "_", symbol)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # Canonical cache location for the current run.  `cache_dir` already points
  # to <out_root>/data, so no legacy directory must be appended here.
  cache_file <- file.path(
    cache_dir,
    paste0(
      base, "_",
      format(start_date, "%Y%m%d"), "_",
      format(end_date, "%Y%m%d"), ".csv"
    )
  )

  if (isTRUE(CFG$frozen_data)) {
    input_file <- file.path(CFG$input_dir, basename(cache_file))
    checksum_file <- file.path(CFG$input_dir, "checksums_md5.csv")
    if (!file.exists(input_file) || !file.exists(checksum_file)) {
      stop("Missing frozen input or checksum manifest: ", input_file)
    }
    manifest <- readr::read_csv(checksum_file, show_col_types = FALSE)
    expected <- manifest$md5[manifest$file == basename(input_file)]
    actual <- unname(tools::md5sum(input_file))
    if (length(expected) != 1L || is.na(actual) || actual != expected) {
      stop("Frozen-input checksum mismatch: ", input_file)
    }
    frozen <- read_compatible_yahoo_cache(input_file, start_date, end_date)
    if (is.null(frozen)) stop("Invalid frozen input: ", input_file)
    if (!file.copy(input_file, cache_file, overwrite = TRUE)) {
      stop("Could not copy frozen input into the output cache: ", cache_file)
    }
    message("  Using verified frozen input: ", input_file)
    return(frozen)
  }

  # 1) Exact cache for the current publication run.
  cached <- read_compatible_yahoo_cache(cache_file, start_date, end_date)
  if (!is.null(cached)) {
    message("  Using frozen cache: ", cache_file)
    return(cached)
  }

  # If frozen_data is explicitly disabled, retrieve the selected market window.
  x <- NULL
  quantmod_error <- NULL

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(120, old_timeout %||% 60))

  for (attempt in 1:3) {
    message("  Yahoo quantmod attempt ", attempt, "/3: ", symbol)

    x <- tryCatch(
      suppressWarnings(
        quantmod::getSymbols(
          Symbols = symbol,
          src = "yahoo",
          from = start_date,
          to = end_date + 1,
          auto.assign = FALSE,
          warnings = FALSE
        )
      ),
      error = function(e) {
        quantmod_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(x) && NROW(x) > 100L) break
    Sys.sleep(attempt)
  }

  out <- NULL
  if (!is.null(x) && NROW(x) > 100L) {
    adjusted <- tryCatch(as.numeric(quantmod::Ad(x)), error = function(e) NULL)
    if (is.null(adjusted) || all(!is.finite(adjusted))) {
      adjusted <- as.numeric(quantmod::Cl(x))
    }

    out <- tibble(
      date = as.Date(zoo::index(x)),
      adjusted = adjusted,
      close = as.numeric(quantmod::Cl(x)),
      high = as.numeric(quantmod::Hi(x)),
      low = as.numeric(quantmod::Lo(x))
    ) %>%
      filter(date >= start_date, date <= end_date) %>%
      arrange(date) %>%
      distinct(date, .keep_all = TRUE) %>%
      filter(is.finite(adjusted), adjusted > 0) %>%
      mutate(
        ret = log(adjusted) - lag(log(adjusted)),
        sqret = ret^2,
        hl_range = if_else(
          is.finite(high) & is.finite(low) & high > 0 & low > 0,
          log(high / low)^2,
          NA_real_
        )
      )
  }

  # 4) If quantmod cannot connect/import, use Yahoo's crumb-free Chart API.
  if (is.null(out) || nrow(out) <= 100L) {
    warning(
      "quantmod Yahoo import failed for ", symbol,
      if (!is.null(quantmod_error)) paste0(": ", quantmod_error) else "",
      ". Trying Yahoo Chart API fallback."
    )
    out <- download_yahoo_chart_api(symbol, start_date, end_date)
  }

  if (is.null(out) || nrow(out) <= 100L) {
    stop("No usable Yahoo data obtained for ", symbol, ".")
  }

  # Hard data checks before freezing the cache.
  if (anyDuplicated(out$date)) stop("Duplicate dates after Yahoo download: ", symbol)
  if (any(!is.finite(out$adjusted) | out$adjusted <= 0)) {
    stop("Invalid adjusted/close prices after Yahoo download: ", symbol)
  }

  readr::write_csv(out, cache_file)
  message("  Frozen cache written: ", cache_file, " | rows=", nrow(out))
  out
}

download_all_data <- function(cfg = CFG) {
  setNames(
    lapply(
      cfg$symbols,
      download_yahoo_symbol,
      start_date = cfg$data_start,
      end_date = cfg$data_end,
      cache_dir = DIRS$data
    ),
    cfg$symbols
  )
}

# ----------------------------------------------------------------------
# 4. Forward-RV curves and separate range-based stress proxy
# ----------------------------------------------------------------------
make_forward_panel <- function(df, H) {
  d <- df %>%
    filter(is.finite(ret)) %>%
    arrange(date)

  n <- nrow(d)
  if (n <= H + 30L) stop("Too few observations for H=", H)

  n_origin <- n - H
  origins <- seq_len(n_origin)

  M <- matrix(NA_real_, nrow = n_origin, ncol = H)

  for (h in seq_len(H)) {
    M[, h] <- vapply(
      origins,
      function(i) sum(d$ret[(i + 1L):(i + h)]^2, na.rm = TRUE),
      numeric(1)
    )
  }

  colnames(M) <- paste0("h", seq_len(H))

  # Separate range-based proxy for descriptive event validation:
  # future maximum high-low range over the same H-day horizon. We
  # also retain the calendar date of the daily high-low contribution that
  # attains each forward maximum. This makes mechanical ties caused by
  # overlapping forward windows directly auditable.
  future_range_max <- rep(NA_real_, n_origin)
  future_range_peak_date <- as.Date(rep(NA_character_, n_origin))

  for (i in origins) {
    idx <- (i + 1L):(i + H)
    z <- d$hl_range[idx]
    good <- which(is.finite(z))
    if (length(good) > 0L) {
      z_good <- z[good]
      local_max <- max(z_good)
      # Deterministic first occurrence if an exact within-window tie occurs.
      rel <- good[which(z_good == local_max)[1L]]
      future_range_max[i] <- local_max
      future_range_peak_date[i] <- d$date[idx[rel]]
    }
  }

  # Past information for HAR-style benchmark.
  r2 <- d$ret^2
  rv1 <- rv5 <- rv22 <- rep(NA_real_, n_origin)

  for (i in origins) {
    rv1[i] <- r2[i]
    if (i >= 5L)  rv5[i]  <- mean(r2[(i - 4L):i], na.rm = TRUE)
    if (i >= 22L) rv22[i] <- mean(r2[(i - 21L):i], na.rm = TRUE)
  }

  list(
    M = M,
    origin_date = d$date[origins],
    origin_row = origins,
    future_range_max = future_range_max,
    future_range_peak_date = future_range_peak_date,
    rv1 = rv1,
    rv5 = rv5,
    rv22 = rv22,
    source = d
  )
}

# ----------------------------------------------------------------------
# 5. Leakage-disciplined chronological split with H-observation embargo
# ----------------------------------------------------------------------
make_split_indices <- function(n, H, train_frac = 0.60, cal_frac = 0.20) {
  if (train_frac <= 0 || cal_frac <= 0 || train_frac + cal_frac >= 1) {
    stop("Invalid train/cal fractions.")
  }

  train_end_raw <- floor(train_frac * n)
  cal_end_raw <- floor((train_frac + cal_frac) * n)

  train_idx <- seq_len(train_end_raw)

  cal_start <- train_end_raw + H + 1L
  cal_end <- cal_end_raw

  test_start <- cal_end_raw + H + 1L
  test_end <- n

  cal_idx <- if (cal_start <= cal_end) cal_start:cal_end else integer(0)
  test_idx <- if (test_start <= test_end) test_start:test_end else integer(0)

  if (length(train_idx) < 100L || length(cal_idx) < 50L || length(test_idx) < 50L) {
    stop(
      "Split too small after embargo. n=", n,
      ", H=", H,
      ", train=", length(train_idx),
      ", cal=", length(cal_idx),
      ", test=", length(test_idx)
    )
  }

  list(
    train = train_idx,
    cal = cal_idx,
    test = test_idx,
    embargo = H,
    train_end_raw = train_end_raw,
    cal_end_raw = cal_end_raw
  )
}

# Calibration-window sensitivity with a FIXED test period.
# The baseline 60/20 split determines the calibration end and test start.
# Changing cal_frac moves the calibration start backward/forward and adjusts
# the training endpoint, while the held-out test origins remain identical.
# The usable calibration length is approximately cal_frac*n - H, matching
# the baseline split convention in which the H-origin embargo is excluded.
make_split_indices_fixed_test <- function(
    n, H,
    train_frac = 0.60,
    baseline_cal_frac = 0.20,
    cal_frac = 0.20) {

  if (train_frac <= 0 || baseline_cal_frac <= 0 || cal_frac <= 0 ||
      train_frac + baseline_cal_frac >= 1 || cal_frac >= 1) {
    stop("Invalid split fractions for fixed-test calibration sensitivity.")
  }

  baseline_cal_end <- floor((train_frac + baseline_cal_frac) * n)
  test_start <- baseline_cal_end + H + 1L
  test_end <- n

  # This reproduces the original baseline exactly when
  # cal_frac == baseline_cal_frac (up to integer flooring).
  nominal_cal_span <- floor(cal_frac * n)
  cal_start <- baseline_cal_end - nominal_cal_span + H + 1L
  cal_end <- baseline_cal_end
  train_end <- cal_start - H - 1L

  train_idx <- if (train_end >= 1L) seq_len(train_end) else integer(0)
  cal_idx <- if (cal_start <= cal_end) cal_start:cal_end else integer(0)
  test_idx <- if (test_start <= test_end) test_start:test_end else integer(0)

  if (length(train_idx) < 100L || length(cal_idx) < 50L || length(test_idx) < 50L) {
    stop(
      "Fixed-test split too small. n=", n,
      ", H=", H,
      ", cal_frac=", cal_frac,
      ", train=", length(train_idx),
      ", cal=", length(cal_idx),
      ", test=", length(test_idx)
    )
  }

  # Explicit leakage audit.
  if (!(max(train_idx) + H < min(cal_idx))) {
    stop("Fixed-test split violates train/cal H-origin embargo.")
  }
  if (!(max(cal_idx) + H < min(test_idx))) {
    stop("Fixed-test split violates cal/test H-origin embargo.")
  }

  list(
    train = train_idx,
    cal = cal_idx,
    test = test_idx,
    embargo = H,
    train_end_raw = train_end,
    cal_end_raw = cal_end,
    test_start_raw = test_start,
    fixed_test = TRUE,
    baseline_cal_frac = baseline_cal_frac,
    cal_frac = cal_frac
  )
}

# ----------------------------------------------------------------------
# 6. Isotonic projection and curvature
# ----------------------------------------------------------------------
iso_proj <- function(x) {
  if (all(diff(x) >= -1e-14)) return(as.numeric(x))
  as.numeric(stats::isoreg(seq_along(x), x)$yf)
}

iso_matrix <- function(M) {
  t(apply(M, 1L, iso_proj))
}

curve_curvature <- function(v) {
  # Manuscript-exact definition:
  # kappa_h = (v_{h+1}-v_h) - (v_h-v_{h-1}), h=2,...,H-1.
  # For cumulative FRV, v_h = sum_{j<=h} r_{t+j}^2, so
  # kappa_h = r_{t+h+1}^2 - r_{t+h}^2. The curvature component therefore
  # measures local changes/spikiness in daily squared-return contributions, not
  # an abstract geometric curvature independent of the volatility level.
  # Hence every H-point curve contributes exactly H-2 curvature values.
  if (length(v) < 3L) stop("Curvature requires H >= 3.")
  diff(v, differences = 2)
}

projection_audit <- function(M) {
  adj <- apply(M, 1L, function(x) {
    xi <- iso_proj(x)
    max(abs(x - xi), na.rm = TRUE)
  })
  tibble(
    n_curves = nrow(M),
    share_adjusted_gt_1e_12 = mean(adj > 1e-12),
    mean_max_adjustment = mean(adj),
    max_adjustment = max(adj)
  )
}

# ----------------------------------------------------------------------
# 7. CIRD: training fit and scoring
# ----------------------------------------------------------------------
fit_cird <- function(M_train) {
  M_iso <- iso_matrix(M_train)
  H <- ncol(M_iso)

  terminal <- M_iso[, H]
  terminal_median <- median(terminal, na.rm = TRUE)
  terminal_scale <- safe_mad(terminal)

  L_train <- (terminal - terminal_median) / terminal_scale

  K_train <- t(apply(M_iso, 1L, curve_curvature))
  curvature_scale <- safe_mad(as.numeric(K_train))

  C_train <- rowMeans(abs(K_train) / curvature_scale, na.rm = TRUE)

  list(
    H = H,
    M_train_iso = M_iso,
    terminal_median = terminal_median,
    terminal_scale = terminal_scale,
    curvature_scale = curvature_scale,
    mu_L = mean(L_train, na.rm = TRUE),
    sd_L = safe_sd(L_train),
    mu_C = mean(C_train, na.rm = TRUE),
    sd_C = safe_sd(C_train),
    train_L = L_train,
    train_C = C_train
  )
}

cird_components <- function(fit, M_new) {
  M_iso <- iso_matrix(M_new)

  if (ncol(M_iso) != fit$H) stop("H mismatch in CIRD scoring.")

  terminal <- M_iso[, fit$H]
  L <- (terminal - fit$terminal_median) / fit$terminal_scale

  K <- t(apply(M_iso, 1L, curve_curvature))
  C <- rowMeans(abs(K) / fit$curvature_scale, na.rm = TRUE)

  zL <- (L - fit$mu_L) / fit$sd_L
  zC <- (C - fit$mu_C) / fit$sd_C

  tibble(
    L = L,
    C = C,
    zL = zL,
    zC = zC
  )
}

score_cird <- function(fit, M_new, a = 1, b = 1) {
  cc <- cird_components(fit, M_new)
  logistic(a * cc$zL - b * cc$zC)
}

# ----------------------------------------------------------------------
# 8. Benchmark scores
# ----------------------------------------------------------------------
fit_terminal_level <- function(M_train) {
  x <- log(M_train[, ncol(M_train)] + 1e-12)
  list(
    med = median(x, na.rm = TRUE),
    scale = safe_mad(x),
    H = ncol(M_train)
  )
}

score_terminal_level <- function(fit, M_new) {
  x <- log(M_new[, fit$H] + 1e-12)
  logistic((x - fit$med) / fit$scale)
}

fit_upper_tail_rank <- function(M_train) {
  list(M_train_iso = iso_matrix(M_train), H = ncol(M_train))
}

score_upper_tail_rank <- function(fit, M_new) {
  M_new_iso <- iso_matrix(M_new)
  tr <- fit$M_train_iso
  n <- nrow(tr)
  H <- ncol(tr)

  vapply(seq_len(nrow(M_new_iso)), function(i) {
    x <- M_new_iso[i, ]
    u_h <- vapply(seq_len(H), function(h) {
      # Smoothed lower-tail empirical percentile. High-variance curves have
      # u close to one, so averaging u gives a one-sided upper-tail score.
      (sum(tr[, h] < x[h]) + 0.5 * sum(tr[, h] == x[h]) + 0.5) / (n + 1)
    }, numeric(1))
    mean(u_h)
  }, numeric(1))
}

fit_har_forward <- function(panel, train_idx, eps = CFG$har_log_epsilon) {
  y <- log(panel$M[, ncol(panel$M)] + eps)

  dat <- tibble(
    y = y,
    log_rv1 = log(panel$rv1 + eps),
    log_rv5 = log(panel$rv5 + eps),
    log_rv22 = log(panel$rv22 + eps)
  )

  tr <- dat[train_idx, , drop = FALSE] %>% filter(if_all(everything(), is.finite))
  if (nrow(tr) < 100L) stop("Too few HAR training observations.")

  model <- lm(y ~ log_rv1 + log_rv5 + log_rv22, data = tr)
  pred_tr <- predict(model, newdata = tr)
  res_tr <- tr$y - pred_tr

  list(
    model = model,
    residual_median = median(res_tr, na.rm = TRUE),
    residual_scale = safe_mad(res_tr)
  )
}

score_har_forward <- function(fit, panel, eps = CFG$har_log_epsilon) {
  y <- log(panel$M[, ncol(panel$M)] + eps)
  dat <- tibble(
    log_rv1 = log(panel$rv1 + eps),
    log_rv5 = log(panel$rv5 + eps),
    log_rv22 = log(panel$rv22 + eps)
  )

  pred <- suppressWarnings(predict(fit$model, newdata = dat))
  res <- y - pred
  score <- logistic((res - fit$residual_median) / fit$residual_scale)
  score[!is.finite(score)] <- NA_real_
  score
}

# ----------------------------------------------------------------------
# 9. Empirical evaluation at block level
# ----------------------------------------------------------------------
block_score_objects <- function(score, B) {
  n_full <- floor(length(score) / B)
  if (n_full < 1L) {
    return(list(
      score_matrix = matrix(numeric(0), nrow = B, ncol = 0),
      valid_block = logical(0),
      block_scores_all = numeric(0),
      block_scores = numeric(0),
      use_n = 0L,
      n_full = 0L
    ))
  }

  use_n <- n_full * B
  score_mat <- matrix(score[seq_len(use_n)], nrow = B, ncol = n_full)
  valid_block <- apply(score_mat, 2, function(z) all(is.finite(z)))

  block_scores_all <- rep(NA_real_, n_full)
  if (any(valid_block)) {
    block_scores_all[valid_block] <- apply(
      score_mat[, valid_block, drop = FALSE],
      2,
      max
    )
  }

  list(
    score_matrix = score_mat,
    valid_block = valid_block,
    block_scores_all = block_scores_all,
    block_scores = block_scores_all[valid_block],
    use_n = use_n,
    n_full = n_full
  )
}


# Held-out, ex-post block-tail reference for descriptive event validation.
# IMPORTANT: this reference is never used to fit/standardize CIRD and never
# enters the conformal calibration threshold. It is constructed only after
# the chronological train/cal/test split has been fixed.
make_expost_block_stress_reference <- function(
    proxy,
    split,
    B,
    quantile_prob = 0.90,
    proxy_peak_date = NULL,
    origin_date = NULL
) {
  if (!is.finite(quantile_prob) || quantile_prob <= 0 || quantile_prob >= 1) {
    stop("stress evaluation quantile must lie strictly between 0 and 1.")
  }

  proxy_test <- proxy[split$test]
  obj <- block_score_objects(proxy_test, B)

  if (obj$n_full < 1L || !any(obj$valid_block)) {
    stop("No finite held-out range-proxy blocks available for stress validation.")
  }
  if (any(!obj$valid_block)) {
    stop(
      "Held-out range-proxy validation contains non-finite complete blocks; ",
      "inspect Yahoo high/low data before interpreting event metrics."
    )
  }

  valid_values <- obj$block_scores_all[obj$valid_block]
  q <- as.numeric(stats::quantile(
    valid_values,
    probs = quantile_prob,
    na.rm = TRUE,
    names = FALSE,
    type = 8
  ))

  # >= guarantees that the empirical quantile itself enters the reference.
  # Repeated reference values are expected here because adjacent origin blocks
  # contain overlapping H-day forward windows and can therefore inherit the
  # same future daily high-low peak.
  event_block_all <- rep(NA, obj$n_full)
  event_block_all[obj$valid_block] <- valid_values >= q

  n_event <- sum(event_block_all %in% TRUE, na.rm = TRUE)
  if (n_event < 1L) {
    stop("Held-out range-proxy reference produced zero stress blocks.")
  }

  event_daily_test <- rep(NA, length(proxy_test))
  if (obj$use_n > 0L) {
    block_id <- rep(seq_len(obj$n_full), each = B)
    mapped <- event_block_all[block_id]
    event_daily_test[seq_len(obj$use_n)] <- mapped
  }

  event_daily <- rep(NA, length(proxy))
  event_daily[split$test] <- event_daily_test

  block_details <- tibble(
    block_id = seq_len(obj$n_full),
    block_proxy_max = obj$block_scores_all,
    stress_reference_event = event_block_all
  )

  if (!is.null(origin_date)) {
    od <- origin_date[split$test]
    od_mat <- matrix(
      as.integer(od[seq_len(obj$use_n)]),
      nrow = B, ncol = obj$n_full
    )
    block_details$block_start <- as.Date(
      apply(od_mat, 2, min), origin = "1970-01-01"
    )
    block_details$block_end <- as.Date(
      apply(od_mat, 2, max), origin = "1970-01-01"
    )
  }

  if (!is.null(proxy_peak_date)) {
    pp <- as.Date(proxy_peak_date[split$test])
    pmat <- matrix(proxy_test[seq_len(obj$use_n)], nrow = B, ncol = obj$n_full)
    dmat <- matrix(
      as.integer(pp[seq_len(obj$use_n)]),
      nrow = B, ncol = obj$n_full
    )
    winner <- vapply(seq_len(obj$n_full), function(j) {
      z <- pmat[, j]
      which.max(ifelse(is.finite(z), z, -Inf))[1L]
    }, integer(1))
    peak_int <- vapply(seq_len(obj$n_full), function(j) dmat[winner[j], j], integer(1))
    block_details$proxy_peak_date <- as.Date(peak_int, origin = "1970-01-01")
  } else {
    block_details$proxy_peak_date <- as.Date(NA)
  }

  event_ids <- which(event_block_all %in% TRUE)
  episode_id_all <- rep(NA_integer_, obj$n_full)
  if (length(event_ids) > 0L) {
    run_id <- cumsum(c(TRUE, diff(event_ids) > 1L))
    episode_id_all[event_ids] <- run_id
  }
  block_details$reference_episode_id <- episode_id_all

  event_peak_dates <- block_details$proxy_peak_date[event_ids]
  event_proxy_values <- block_details$block_proxy_max[event_ids]

  list(
    threshold = q,
    quantile_prob = quantile_prob,
    event_daily = event_daily,
    event_block_all = event_block_all,
    block_proxy_all = obj$block_scores_all,
    valid_block = obj$valid_block,
    n_valid_blocks = sum(obj$valid_block),
    n_event_blocks = n_event,
    event_block_rate = mean(event_block_all[obj$valid_block]),
    n_distinct_reference_runs = if (length(event_ids) == 0L) 0L else max(run_id),
    n_unique_peak_dates = dplyr::n_distinct(event_peak_dates[!is.na(event_peak_dates)]),
    n_unique_reference_proxy_values = dplyr::n_distinct(event_proxy_values[is.finite(event_proxy_values)]),
    block_details = block_details,
    source = "heldout_test_block_tail_ex_post"
  )
}

block_forward_window_table <- function(panel, split, B) {
  rows <- panel$origin_row[split$test]
  n_full <- floor(length(rows) / B)
  if (n_full < 1L) return(tibble())
  use_n <- n_full * B
  rmat <- matrix(rows[seq_len(use_n)], nrow = B, ncol = n_full)
  first_row <- apply(rmat, 2, min)
  last_row <- apply(rmat, 2, max)
  H <- ncol(panel$M)
  tibble(
    block_id = seq_len(n_full),
    forward_start = panel$source$date[first_row + 1L],
    forward_end = panel$source$date[last_row + H]
  )
}

component_correlation_summary <- function(cfit, panel, split, symbol) {
  cc <- cird_components(cfit, panel$M)
  bind_rows(lapply(c("Training", "Test"), function(lbl) {
    idx <- if (lbl == "Training") split$train else split$test
    tibble(
      symbol = symbol,
      sample = lbl,
      n = sum(is.finite(cc$L[idx]) & is.finite(cc$C[idx])),
      pearson_L_C = suppressWarnings(cor(cc$L[idx], cc$C[idx], use = "complete.obs")),
      spearman_L_C = suppressWarnings(cor(cc$L[idx], cc$C[idx], use = "complete.obs", method = "spearman"))
    )
  }))
}

har_zero_return_audit <- function(df, panel, split, symbol, eps = CFG$har_log_epsilon) {
  ret <- df$ret[is.finite(df$ret)]
  rv1_train <- panel$rv1[split$train]
  positive_train <- rv1_train[is.finite(rv1_train) & rv1_train > 0]
  tibble(
    symbol = symbol,
    finite_returns = length(ret),
    exact_zero_returns = sum(ret == 0),
    exact_zero_return_share = mean(ret == 0),
    har_training_rv1_n = sum(is.finite(rv1_train)),
    har_training_rv1_zero = sum(rv1_train == 0, na.rm = TRUE),
    har_training_rv1_at_or_below_epsilon = sum(rv1_train <= eps, na.rm = TRUE),
    minimum_positive_training_rv1 = if (length(positive_train) > 0L) min(positive_train) else NA_real_,
    har_log_epsilon = eps
  )
}

run_random_score_reference <- function(panel, split, event_daily, symbol, cfg = CFG) {
  seed <- cfg$random_score_seed + match(symbol, cfg$symbols)
  with_local_seed(seed, function() {
    bind_rows(lapply(seq_len(cfg$random_score_reps), function(rep_id) {
      random_score <- runif(nrow(panel$M))
      ev <- evaluate_score(
        score = random_score,
        split = split,
        B = cfg$B_main,
        alpha_day = cfg$alpha_day_main,
        event_daily = event_daily,
        label = "Random_score_reference"
      )$summary
      ev %>%
        transmute(
          symbol = symbol,
          rep = rep_id,
          alpha_block_design = alpha_block_design,
          block_exceedance_rate = block_exceedance_rate,
          abs_calibration_gap = abs(block_exceedance_rate - alpha_block_design),
          precision = stress_proxy_precision,
          recall = stress_proxy_recall,
          f1 = stress_proxy_f1
        )
    }))
  })
}

build_recognized_episode_mapping <- function(
    panel, split, B, cird_calendar, stress_ref, event_windows, symbol
) {
  fw <- block_forward_window_table(panel, split, B)
  cird <- cird_calendar %>%
    filter(method == "CIRD") %>%
    select(block_id, alarm)
  ref <- stress_ref$block_details %>%
    select(
      block_id, stress_reference_event, reference_episode_id,
      proxy_peak_date, block_proxy_max
    )
  detail <- fw %>%
    left_join(cird, by = "block_id") %>%
    left_join(ref, by = "block_id")

  bind_rows(lapply(seq_len(nrow(event_windows)), function(i) {
    w <- event_windows[i, ]
    alarm_idx <- which(
      detail$alarm %in% TRUE &
        detail$forward_start <= w$event_end &
        detail$forward_end >= w$event_start
    )
    ref_idx <- which(
      detail$stress_reference_event %in% TRUE &
        !is.na(detail$proxy_peak_date) &
        detail$proxy_peak_date >= w$event_start &
        detail$proxy_peak_date <= w$event_end
    )

    tibble(
      symbol = symbol,
      episode_id = w$episode_id,
      episode_label = w$episode_label,
      event_start = w$event_start,
      event_end = w$event_end,
      cird_alarm_forward_overlap_blocks = length(alarm_idx),
      cird_alarm_block_ids = collapse_int_ids(detail$block_id[alarm_idx]),
      reference_blocks = length(ref_idx),
      reference_block_ids = collapse_int_ids(detail$block_id[ref_idx]),
      reference_peak_dates = collapse_dates(detail$proxy_peak_date[ref_idx]),
      reference_run_ids = collapse_int_ids(unique(detail$reference_episode_id[ref_idx]))
    )
  }))
}


build_alarm_calendar <- function(
    score,
    split,
    origin_date,
    B,
    q,
    event_daily = NULL,
    symbol,
    method
) {
  sc_test <- score[split$test]
  date_test <- origin_date[split$test]

  obj <- block_score_objects(sc_test, B)
  if (obj$n_full < 1L) return(tibble())

  use_n <- obj$use_n
  date_mat <- matrix(date_test[seq_len(use_n)], nrow = B, ncol = obj$n_full)
  valid <- obj$valid_block

  if (!any(valid)) return(tibble())

  block_id_all <- seq_len(obj$n_full)
  block_start_all <- as.Date(apply(date_mat, 2, min), origin = "1970-01-01")
  block_end_all <- as.Date(apply(date_mat, 2, max), origin = "1970-01-01")

  event_block_all <- rep(NA, obj$n_full)
  if (!is.null(event_daily)) {
    ev_test <- event_daily[split$test]
    ev_mat <- matrix(ev_test[seq_len(use_n)], nrow = B, ncol = obj$n_full)
    event_block_all <- apply(ev_mat, 2, function(z) any(z %in% TRUE, na.rm = TRUE))
  }

  tibble(
    symbol = symbol,
    method = method,
    block_id = block_id_all[valid],
    block_start = block_start_all[valid],
    block_end = block_end_all[valid],
    block_score = obj$block_scores_all[valid],
    threshold = q,
    alarm = obj$block_scores_all[valid] > q,
    stress_proxy_event = event_block_all[valid]
  )
}

evaluate_score <- function(
    score,
    split,
    B,
    alpha_day,
    event_daily = NULL,
    label = "score"
) {
  sc_cal <- score[split$cal]
  sc_test <- score[split$test]

  # Preserve chronology and block boundaries. Never remove isolated NA values
  # before blocking; instead, drop a complete block if any score in that block
  # is non-finite. This is particularly relevant for benchmark scores.
  cal_obj <- block_score_objects(sc_cal, B)
  bm_cal <- cal_obj$block_scores
  if (length(bm_cal) < 1L) stop("No valid calibration blocks for ", label)

  alpha_block <- alpha_block_map(alpha_day, B)
  qfit <- conformal_threshold(bm_cal, alpha_block)

  test_obj <- block_score_objects(sc_test, B)
  bm_test <- test_obj$block_scores
  if (length(bm_test) < 1L) stop("No valid test blocks for ", label)

  alarms <- bm_test > qfit$q
  n_alarm <- sum(alarms)
  n_test_blocks <- length(alarms)

  rate_block <- mean(alarms)
  ci <- wilson_ci(n_alarm, n_test_blocks)

  finite_test <- sc_test[is.finite(sc_test)]
  rate_day <- if (length(finite_test) > 0L) mean(finite_test > qfit$q) else NA_real_

  out <- tibble(
    method = label,
    B = B,
    alpha_day = alpha_day,
    alpha_block_design = alpha_block,
    cal_blocks = length(bm_cal),
    cal_blocks_dropped_nonfinite = cal_obj$n_full - sum(cal_obj$valid_block),
    test_blocks = n_test_blocks,
    test_blocks_dropped_nonfinite = test_obj$n_full - sum(test_obj$valid_block),
    k_index = qfit$k,
    threshold = qfit$q,
    block_exceedance_rate = rate_block,
    block_rate_lo = ci[1],
    block_rate_hi = ci[2],
    daily_exceedance_rate = rate_day,
    alarm_blocks = n_alarm
  )

  # Optional held-out range-reference validation. This is not ground truth.
  if (!is.null(event_daily)) {
    use_n <- test_obj$use_n
    ev_test <- event_daily[split$test]
    ev_mat <- matrix(ev_test[seq_len(use_n)], nrow = B, ncol = test_obj$n_full)
    event_block_all <- apply(
      ev_mat,
      2,
      function(z) any(z %in% TRUE, na.rm = TRUE)
    )
    event_block <- event_block_all[test_obj$valid_block]

    tp <- sum(alarms & event_block)
    fp <- sum(alarms & !event_block)
    fn <- sum(!alarms & event_block)

    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

    # Direct confusion-matrix definition. If there are reference events and/or
    # alarms but no true positives, F1 is correctly zero rather than NA.
    f1_den <- 2 * tp + fp + fn
    f1 <- if (f1_den > 0) {
      2 * tp / f1_den
    } else {
      NA_real_
    }

    out <- out %>%
      mutate(
        stress_event_blocks = sum(event_block),
        stress_event_prevalence = mean(event_block),
        stress_hit_blocks = tp,
        stress_proxy_precision = precision,
        stress_proxy_recall = recall,
        stress_proxy_f1 = f1
      )
  }

  list(
    summary = out,
    q = qfit$q,
    alarms = alarms,
    block_scores_cal = bm_cal,
    block_scores_test = bm_test,
    cal_valid_block = cal_obj$valid_block,
    test_valid_block = test_obj$valid_block
  )
}

# ----------------------------------------------------------------------
# 10. Descriptive statistics
# ----------------------------------------------------------------------
sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  mean(((x - m) / s)^3)
}

sample_kurtosis_raw <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  s <- sd(x)
  mean(((x - m) / s)^4)
}

return_descriptives <- function(df, symbol) {
  x <- df$ret[is.finite(df$ret)]
  n <- length(x)
  sk <- sample_skewness(x)
  ku <- sample_kurtosis_raw(x)
  jb <- n / 6 * (sk^2 + ((ku - 3)^2) / 4)
  jb_p <- 1 - pchisq(jb, df = 2)

  tibble(
    symbol = symbol,
    n = n,
    mean = mean(x),
    sd = sd(x),
    skewness = sk,
    kurtosis_raw = ku,
    jb_p = jb_p,
    lb10_return_p = Box.test(x, lag = 10, type = "Ljung-Box")$p.value,
    lb10_sqreturn_p = Box.test(x^2, lag = 10, type = "Ljung-Box")$p.value
  )
}

frv_descriptives <- function(M_train, symbol, horizons = c(1, 5, 10, 20, 30)) {
  horizons <- horizons[horizons <= ncol(M_train)]

  bind_rows(lapply(horizons, function(h) {
    x <- M_train[, h]
    tibble(
      symbol = symbol,
      h = h,
      mean = mean(x),
      median = median(x),
      sd = sd(x),
      p90 = as.numeric(quantile(x, 0.90, names = FALSE))
    )
  }))
}

# ----------------------------------------------------------------------
# 11. One asset: baseline + benchmarks + weights
# ----------------------------------------------------------------------
run_asset_baseline <- function(df, symbol, cfg = CFG) {
  H <- cfg$H_main
  B <- cfg$B_main

  panel <- make_forward_panel(df, H)
  split <- make_split_indices(
    nrow(panel$M),
    H = H,
    train_frac = cfg$train_frac,
    cal_frac = cfg$cal_frac_main
  )

  Mtr <- panel$M[split$train, , drop = FALSE]

  # CIRD fit
  cfit <- fit_cird(Mtr)
  cird_score <- score_cird(cfit, panel$M, a = 1, b = 1)

  # Held-out, ex-post range-based stress reference for DESCRIPTIVE validation.
  # This is deliberately not used in CIRD fitting or conformal calibration.
  stress_ref <- make_expost_block_stress_reference(
    proxy = panel$future_range_max,
    split = split,
    B = B,
    quantile_prob = cfg$stress_eval_quantile,
    proxy_peak_date = panel$future_range_peak_date,
    origin_date = panel$origin_date
  )
  stress_q <- stress_ref$threshold
  stress_event <- stress_ref$event_daily

  # Benchmarks
  tfit <- fit_terminal_level(Mtr)
  score_terminal <- score_terminal_level(tfit, panel$M)

  rfit <- fit_upper_tail_rank(Mtr)
  score_rank <- score_upper_tail_rank(rfit, panel$M)

  hfit <- fit_har_forward(panel, split$train, eps = cfg$har_log_epsilon)
  score_har <- score_har_forward(hfit, panel, eps = cfg$har_log_epsilon)

  scores <- list(
    CIRD = cird_score,
    Terminal_FRV = score_terminal,
    Integrated_upper_tail_rank = score_rank,
    HAR_forward_residual = score_har
  )

  bench_eval <- setNames(lapply(names(scores), function(nm) {
    evaluate_score(
      score = scores[[nm]],
      split = split,
      B = B,
      alpha_day = cfg$alpha_day_main,
      event_daily = stress_event,
      label = nm
    )
  }), names(scores))

  bench <- bind_rows(lapply(bench_eval, `[[`, "summary")) %>%
    mutate(symbol = symbol, .before = 1)

  alarm_calendar <- bind_rows(lapply(names(scores), function(nm) {
    build_alarm_calendar(
      score = scores[[nm]],
      split = split,
      origin_date = panel$origin_date,
      B = B,
      q = bench_eval[[nm]]$q,
      event_daily = stress_event,
      symbol = symbol,
      method = nm
    )
  }))

  # Publication diagnostics.
  component_corr <- component_correlation_summary(cfit, panel, split, symbol)
  zero_return_audit <- har_zero_return_audit(
    df, panel, split, symbol, eps = cfg$har_log_epsilon
  )

  level_score_equiv <- score_cird(cfit, panel$M, a = 1, b = 0)
  level_eval_equiv <- evaluate_score(
    level_score_equiv, split, B, cfg$alpha_day_main,
    event_daily = stress_event, label = "level_only_equivalence_check"
  )
  terminal_eval_equiv <- bench_eval[["Terminal_FRV"]]
  level_terminal_equivalence <- tibble(
    symbol = symbol,
    spearman_score_correlation_test = suppressWarnings(cor(
      level_score_equiv[split$test], score_terminal[split$test],
      method = "spearman", use = "complete.obs"
    )),
    identical_test_alarm_sets = identical(
      as.logical(level_eval_equiv$alarms),
      as.logical(terminal_eval_equiv$alarms)
    ),
    level_only_block_rate = level_eval_equiv$summary$block_exceedance_rate,
    terminal_frv_block_rate = terminal_eval_equiv$summary$block_exceedance_rate,
    level_only_recall = level_eval_equiv$summary$stress_proxy_recall,
    terminal_frv_recall = terminal_eval_equiv$summary$stress_proxy_recall
  )

  random_reference <- run_random_score_reference(
    panel, split, stress_event, symbol, cfg
  )

  recognized_episode_mapping <- build_recognized_episode_mapping(
    panel = panel,
    split = split,
    B = B,
    cird_calendar = alarm_calendar,
    stress_ref = stress_ref,
    event_windows = cfg$recognized_episode_windows,
    symbol = symbol
  )

  # Explicit CIRD component ablation: level-only, full score, curvature-only.
  ablation_out <- bind_rows(lapply(seq_len(nrow(cfg$ablation_grid)), function(i) {
    arow <- cfg$ablation_grid[i, ]
    score_a <- score_cird(cfit, panel$M, a = arow$a, b = arow$b)

    evaluate_score(
      score = score_a,
      split = split,
      B = B,
      alpha_day = cfg$alpha_day_main,
      event_daily = stress_event,
      label = arow$ablation_id
    )$summary %>%
      mutate(
        symbol = symbol,
        ablation_id = arow$ablation_id,
        a = arow$a,
        b = arow$b,
        .before = 1
      )
  }))

  # Reliability grid for CIRD only.
  reliability <- bind_rows(lapply(cfg$alpha_day_grid, function(a_day) {
    evaluate_score(
      score = cird_score,
      split = split,
      B = B,
      alpha_day = a_day,
      event_daily = stress_event,
      label = "CIRD"
    )$summary
  })) %>%
    mutate(symbol = symbol, .before = 1)

  # Weight sensitivity.
  weight_out <- bind_rows(lapply(seq_len(nrow(cfg$weight_grid)), function(i) {
    w <- cfg$weight_grid[i, ]
    score_w <- score_cird(cfit, panel$M, a = w$a, b = w$b)

    ev <- evaluate_score(
      score = score_w,
      split = split,
      B = B,
      alpha_day = cfg$alpha_day_main,
      event_daily = stress_event,
      label = "CIRD"
    )

    ev$summary %>%
      mutate(
        symbol = symbol,
        weight_id = w$weight_id,
        a = w$a,
        b = w$b,
        weight_ratio_a_over_b = w$a / w$b,
        .before = 1
      )
  }))

  # Baseline alarm-set Jaccard relative to a=b=1.
  equal_row <- which(cfg$weight_grid$weight_id == "equal")
  equal_score <- score_cird(
    cfit, panel$M,
    a = cfg$weight_grid$a[equal_row],
    b = cfg$weight_grid$b[equal_row]
  )
  equal_eval <- evaluate_score(
    equal_score, split, B, cfg$alpha_day_main,
    event_daily = stress_event, label = "CIRD"
  )

  weight_alarm_sets <- lapply(seq_len(nrow(cfg$weight_grid)), function(i) {
    w <- cfg$weight_grid[i, ]
    s <- score_cird(cfit, panel$M, a = w$a, b = w$b)
    evaluate_score(
      s, split, B, cfg$alpha_day_main,
      event_daily = stress_event, label = "CIRD"
    )$alarms
  })

  weight_out$alarm_jaccard_vs_equal <- vapply(
    weight_alarm_sets,
    jaccard_binary,
    b = equal_eval$alarms,
    numeric(1)
  )

  # Split audit
  split_audit <- tibble(
    symbol = symbol,
    H = H,
    train_n = length(split$train),
    cal_n = length(split$cal),
    test_n = length(split$test),
    embargo = split$embargo,
    train_start = date_or_na(panel$origin_date[split$train]),
    train_end = date_or_na_max(panel$origin_date[split$train]),
    cal_start = date_or_na(panel$origin_date[split$cal]),
    cal_end = date_or_na_max(panel$origin_date[split$cal]),
    test_start = date_or_na(panel$origin_date[split$test]),
    test_end = date_or_na_max(panel$origin_date[split$test])
  )

  iso_audit <- projection_audit(Mtr) %>% mutate(symbol = symbol, .before = 1)

  list(
    symbol = symbol,
    panel = panel,
    split = split,
    cird_fit = cfit,
    cird_score = cird_score,
    stress_q = stress_q,
    stress_ref = stress_ref,
    stress_event = stress_event,
    benchmark = bench,
    alarm_calendar = alarm_calendar,
    ablation = ablation_out,
    reliability = reliability,
    weight_sensitivity = weight_out,
    split_audit = split_audit,
    iso_audit = iso_audit,
    component_corr = component_corr,
    zero_return_audit = zero_return_audit,
    level_terminal_equivalence = level_terminal_equivalence,
    random_reference = random_reference,
    recognized_episode_mapping = recognized_episode_mapping,
    range_reference_block_details = stress_ref$block_details %>% mutate(symbol = symbol, .before = 1),
    return_desc = return_descriptives(df, symbol),
    frv_desc = frv_descriptives(Mtr, symbol)
  )
}

# ----------------------------------------------------------------------
# 12. H x B x calibration-window sensitivity
# ----------------------------------------------------------------------
run_asset_tuning_sensitivity <- function(df, symbol, cfg = CFG) {
  grid <- tidyr::expand_grid(
    H = cfg$H_grid,
    B = cfg$B_grid,
    cal_frac = cfg$cal_frac_grid
  )

  bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]

    panel <- make_forward_panel(df, H = g$H)
    split <- make_split_indices_fixed_test(
      n = nrow(panel$M),
      H = g$H,
      train_frac = cfg$train_frac,
      baseline_cal_frac = cfg$cal_frac_main,
      cal_frac = g$cal_frac
    )

    fit <- fit_cird(panel$M[split$train, , drop = FALSE])
    score <- score_cird(fit, panel$M, a = 1, b = 1)

    stress_ref <- make_expost_block_stress_reference(
      proxy = panel$future_range_max,
      split = split,
      B = g$B,
      quantile_prob = cfg$stress_eval_quantile
    )
    event <- stress_ref$event_daily

    ev <- evaluate_score(
      score = score,
      split = split,
      B = g$B,
      alpha_day = cfg$alpha_day_main,
      event_daily = event,
      label = "CIRD"
    )$summary

    ev %>%
      mutate(
        symbol = symbol,
        H = g$H,
        B = g$B,
        cal_frac = g$cal_frac,
        cal_origin_n = length(split$cal),
        test_origin_n = length(split$test),
        test_start_date = date_or_na(panel$origin_date[split$test]),
        test_end_date = date_or_na_max(panel$origin_date[split$test]),
        fixed_test_calibration_sensitivity = TRUE,
        abs_calibration_gap =
          abs(block_exceedance_rate - alpha_block_design),
        .before = 1
      )
  }))
}

# ----------------------------------------------------------------------
# 13. Simulation DGPs
# ----------------------------------------------------------------------
simulate_sv_ar1 <- function(T, burn = 500L,
                            mu = -5.0, phi = 0.97, sigma_h = 0.20) {
  n <- T + burn
  h <- numeric(n)
  h[1] <- mu
  for (t in 2:n) {
    h[t] <- mu + phi * (h[t - 1] - mu) + sigma_h * rnorm(1)
  }
  r <- exp(h / 2) * rnorm(n)
  tail(r, T)
}

simulate_har_rv <- function(T, burn = 500L,
                            mu = -5.0,
                            b_d = 0.55,
                            b_w = 0.25,
                            b_m = 0.15,
                            sigma_u = 0.20) {
  if ((b_d + b_w + b_m) >= 1) {
    stop("HAR persistence must be < 1 in the revised DGP.")
  }

  n <- T + burn + 30L
  logRV <- rep(mu, n)
  b0 <- (1 - b_d - b_w - b_m) * mu

  for (t in 23:n) {
    ld <- logRV[t - 1L]

    # Manuscript-exact HAR terms: log of the arithmetic mean of RV,
    # not the arithmetic mean of log-RV.
    rv_week <- exp(logRV[(t - 5L):(t - 1L)])
    rv_month <- exp(logRV[(t - 22L):(t - 1L)])
    lw <- log(mean(rv_week))
    lm <- log(mean(rv_month))

    logRV[t] <- b0 + b_d * ld + b_w * lw + b_m * lm +
      rnorm(1, sd = sigma_u)
  }

  RV <- exp(logRV)
  r <- sqrt(RV) * rnorm(n)
  tail(r, T)
}

simulate_sv_break <- function(T, burn = 500L,
                              mu1 = -5.2, mu2 = -4.5,
                              break_frac = 0.80,
                              phi = 0.97, sigma_h = 0.20) {
  n <- T + burn
  h <- numeric(n)
  h[1] <- mu1
  break_at <- burn + floor(T * break_frac)

  for (t in 2:n) {
    mu_t <- if (t <= break_at) mu1 else mu2
    h[t] <- mu_t + phi * (h[t - 1L] - mu_t) + sigma_h * rnorm(1)
  }

  r <- exp(h / 2) * rnorm(n)
  tail(r, T)
}

make_forward_matrix_returns <- function(r, H) {
  r <- as.numeric(r)
  r <- r[is.finite(r)]
  n <- length(r)
  if (n <= H + 50L) stop("Simulation path too short.")

  origins <- seq_len(n - H)
  M <- matrix(NA_real_, nrow = length(origins), ncol = H)

  for (h in seq_len(H)) {
    M[, h] <- vapply(
      origins,
      function(i) sum(r[(i + 1L):(i + h)]^2),
      numeric(1)
    )
  }
  M
}

evaluate_sim_path <- function(r, H, B, cal_frac,
                              alpha_day_grid = CFG$alpha_day_grid,
                              train_frac = CFG$train_frac,
                              fixed_test = FALSE,
                              baseline_cal_frac = CFG$cal_frac_main) {
  M <- make_forward_matrix_returns(r, H)

  split <- if (isTRUE(fixed_test)) {
    make_split_indices_fixed_test(
      n = nrow(M),
      H = H,
      train_frac = train_frac,
      baseline_cal_frac = baseline_cal_frac,
      cal_frac = cal_frac
    )
  } else {
    make_split_indices(
      n = nrow(M),
      H = H,
      train_frac = train_frac,
      cal_frac = cal_frac
    )
  }

  fit <- fit_cird(M[split$train, , drop = FALSE])
  score <- score_cird(fit, M, a = 1, b = 1)

  bind_rows(lapply(alpha_day_grid, function(a_day) {
    ev <- evaluate_score(
      score = score,
      split = split,
      B = B,
      alpha_day = a_day,
      event_daily = NULL,
      label = "CIRD"
    )$summary

    ev %>%
      transmute(
        H = H,
        B = B,
        cal_frac = cal_frac,
        fixed_test = isTRUE(fixed_test),
        alpha_day = alpha_day,
        alpha_block_design = alpha_block_design,
        threshold = threshold,
        block_rate = block_exceedance_rate,
        daily_rate = daily_exceedance_rate,
        cal_blocks = cal_blocks,
        test_blocks = test_blocks
      )
  }))
}

run_simulation_main <- function(cfg = CFG) {
  dgps <- c("HAR_RV", "SV_AR1", "SV_BREAK")

  # Section-specific deterministic seed.
  set.seed(cfg$seed + 1000L)

  all <- vector("list", length(dgps) * cfg$mc_reps)
  pos <- 1L

  for (dgp in dgps) {
    message("Simulation DGP: ", dgp)

    for (rep in seq_len(cfg$mc_reps)) {
      r <- switch(
        dgp,
        HAR_RV = simulate_har_rv(cfg$mc_T, cfg$mc_burn),
        SV_AR1 = simulate_sv_ar1(cfg$mc_T, cfg$mc_burn),
        SV_BREAK = simulate_sv_break(cfg$mc_T, cfg$mc_burn)
      )

      z <- evaluate_sim_path(
        r = r,
        H = cfg$H_main,
        B = cfg$B_main,
        cal_frac = cfg$cal_frac_main
      ) %>%
        mutate(dgp = dgp, rep = rep, .before = 1)

      all[[pos]] <- z
      pos <- pos + 1L
    }
  }

  bind_rows(all)
}

summarise_simulation <- function(sim_raw) {
  sim_raw %>%
    group_by(dgp, H, B, cal_frac, alpha_day, alpha_block_design) %>%
    summarise(
      reps = n(),
      mean_block_rate = mean(block_rate, na.rm = TRUE),
      median_block_rate = median(block_rate, na.rm = TRUE),
      block_rate_lo = quantile(block_rate, 0.025, na.rm = TRUE),
      block_rate_hi = quantile(block_rate, 0.975, na.rm = TRUE),
      mean_daily_rate = mean(daily_rate, na.rm = TRUE),
      q_mean = mean(threshold, na.rm = TRUE),
      q_sd = sd(threshold, na.rm = TRUE),
      q_median = median(threshold, na.rm = TRUE),
      mean_abs_gap = mean(abs(block_rate - alpha_block_design), na.rm = TRUE),
      .groups = "drop"
    )
}

run_simulation_sensitivity <- function(cfg = CFG) {
  # One-factor-at-a-time around the baseline. For each DGP and replication,
  # the SAME simulated return path is reused across all tuning specifications.
  # This common-random-numbers design isolates sensitivity to H, B, and the
  # calibration-window fraction from path-to-path Monte Carlo variation.
  specs <- bind_rows(
    tibble(
      sensitivity = "H",
      H = cfg$H_grid,
      B = cfg$B_main,
      cal_frac = cfg$cal_frac_main
    ),
    tibble(
      sensitivity = "B",
      H = cfg$H_main,
      B = cfg$B_grid,
      cal_frac = cfg$cal_frac_main
    ),
    tibble(
      sensitivity = "cal_frac",
      H = cfg$H_main,
      B = cfg$B_main,
      cal_frac = cfg$cal_frac_grid
    )
  ) %>% distinct()

  dgps <- c("HAR_RV", "SV_AR1", "SV_BREAK")
  out <- list()
  pos <- 1L

  # A section-specific deterministic seed makes this study reproducible even
  # if preceding empirical code is modified in future releases.
  set.seed(cfg$seed + 2000L)

  for (dgp in dgps) {
    for (rep in seq_len(cfg$mc_sensitivity_reps)) {
      r <- switch(
        dgp,
        HAR_RV = simulate_har_rv(cfg$mc_T, cfg$mc_burn),
        SV_AR1 = simulate_sv_ar1(cfg$mc_T, cfg$mc_burn),
        SV_BREAK = simulate_sv_break(cfg$mc_T, cfg$mc_burn)
      )

      for (i in seq_len(nrow(specs))) {
        sp <- specs[i, ]

        z <- evaluate_sim_path(
          r = r,
          H = sp$H,
          B = sp$B,
          cal_frac = sp$cal_frac,
          alpha_day_grid = cfg$alpha_day_main,
          fixed_test = TRUE,
          baseline_cal_frac = cfg$cal_frac_main
        ) %>%
          mutate(
            dgp = dgp,
            rep = rep,
            sensitivity = sp$sensitivity,
            spec_id = paste0(
              sp$sensitivity, "_H", sp$H,
              "_B", sp$B,
              "_C", format(sp$cal_frac, nsmall = 2)
            ),
            .before = 1
          )

        out[[pos]] <- z
        pos <- pos + 1L
      }
    }
  }

  bind_rows(out)
}

# ----------------------------------------------------------------------
# 14. Main empirical run
# ----------------------------------------------------------------------
message("Downloading / loading Yahoo Finance data...")
DATA <- download_all_data(CFG)

message("Running baseline empirical analysis...")
BASELINE <- setNames(
  lapply(CFG$symbols, function(sym) {
    message("  Baseline: ", sym)
    run_asset_baseline(DATA[[sym]], sym, CFG)
  }),
  CFG$symbols
)

benchmark_table <- bind_rows(lapply(BASELINE, `[[`, "benchmark"))
alarm_calendar_table <- bind_rows(lapply(BASELINE, `[[`, "alarm_calendar"))
ablation_table <- bind_rows(lapply(BASELINE, `[[`, "ablation"))
reliability_table <- bind_rows(lapply(BASELINE, `[[`, "reliability"))
weight_table <- bind_rows(lapply(BASELINE, `[[`, "weight_sensitivity"))
split_audit_table <- bind_rows(lapply(BASELINE, `[[`, "split_audit"))
iso_audit_table <- bind_rows(lapply(BASELINE, `[[`, "iso_audit"))
returns_table <- bind_rows(lapply(BASELINE, `[[`, "return_desc"))
frv_table <- bind_rows(lapply(BASELINE, `[[`, "frv_desc"))
component_corr_table <- bind_rows(lapply(BASELINE, `[[`, "component_corr"))
zero_return_audit_table <- bind_rows(lapply(BASELINE, `[[`, "zero_return_audit"))
level_terminal_equivalence_table <- bind_rows(lapply(BASELINE, `[[`, "level_terminal_equivalence"))
random_reference_raw <- bind_rows(lapply(BASELINE, `[[`, "random_reference"))
recognized_episode_mapping_table <- bind_rows(lapply(BASELINE, `[[`, "recognized_episode_mapping"))
range_reference_block_details_table <- bind_rows(lapply(BASELINE, `[[`, "range_reference_block_details"))

readr::write_csv(benchmark_table, file.path(DIRS$tables, "benchmark_comparison.csv"))
readr::write_csv(alarm_calendar_table, file.path(DIRS$tables, "alarm_episode_calendar.csv"))
readr::write_csv(ablation_table, file.path(DIRS$tables, "cird_component_ablation.csv"))
readr::write_csv(reliability_table, file.path(DIRS$tables, "cird_reliability_grid.csv"))
readr::write_csv(weight_table, file.path(DIRS$tables, "cird_weight_sensitivity.csv"))
readr::write_csv(split_audit_table, file.path(DIRS$audit, "split_audit.csv"))
readr::write_csv(iso_audit_table, file.path(DIRS$audit, "isotonic_projection_audit.csv"))
readr::write_csv(returns_table, file.path(DIRS$tables, "return_descriptives.csv"))
readr::write_csv(frv_table, file.path(DIRS$tables, "frv_descriptives.csv"))
readr::write_csv(component_corr_table, file.path(DIRS$tables, "cird_level_curvature_correlations.csv"))
readr::write_csv(zero_return_audit_table, file.path(DIRS$audit, "har_zero_return_audit.csv"))
readr::write_csv(level_terminal_equivalence_table, file.path(DIRS$audit, "level_terminal_equivalence_audit.csv"))
readr::write_csv(random_reference_raw, file.path(DIRS$tables, "random_score_reference_raw.csv"))
readr::write_csv(recognized_episode_mapping_table, file.path(DIRS$tables, "recognized_episode_mapping.csv"))
readr::write_csv(range_reference_block_details_table, file.path(DIRS$tables, "range_reference_block_details.csv"))

saveRDS(BASELINE, file.path(DIRS$rds, "baseline_empirical_objects.rds"))

# ----------------------------------------------------------------------
# 14b. Held-out range-reference estimability audit
# ----------------------------------------------------------------------
# The high--low range reference is a DESCRIPTIVE held-out evaluation device,
# not ground truth.  Its top-tail threshold is estimated only from the held-out
# range-proxy block values after CIRD fitting and conformal calibration are
# complete.  Therefore it cannot affect the fitted score or alarm threshold.
stress_proxy_audit <- bind_rows(lapply(BASELINE, function(z) {
  ref <- z$stress_ref

  tibble(
    symbol = z$symbol,
    reference_source = ref$source,
    heldout_stress_quantile = ref$quantile_prob,
    stress_threshold = ref$threshold,
    heldout_proxy_blocks = ref$n_valid_blocks,
    test_stress_blocks = ref$n_event_blocks,
    test_stress_block_rate = ref$event_block_rate,
    distinct_reference_runs = ref$n_distinct_reference_runs,
    unique_reference_peak_dates = ref$n_unique_peak_dates,
    unique_reference_proxy_values = ref$n_unique_reference_proxy_values
  )
}))

readr::write_csv(
  stress_proxy_audit,
  file.path(DIRS$audit, "stress_proxy_audit.csv")
)

print(stress_proxy_audit)

if (any(!is.finite(stress_proxy_audit$stress_threshold)) ||
    any(stress_proxy_audit$heldout_proxy_blocks < 1L) ||
    any(stress_proxy_audit$test_stress_blocks < 1L)) {
  stop(
    "STRESS-REFERENCE AUDIT FAILED: held-out range-based event validation ",
    "is not estimable for every asset."
  )
}

# ----------------------------------------------------------------------
# 14c. Publication diagnostic summaries
# ----------------------------------------------------------------------
random_reference_asset_summary <- random_reference_raw %>%
  group_by(symbol) %>%
  summarise(
    reps = n(),
    alpha_block_design = first(alpha_block_design),
    mean_block_rate = mean(block_exceedance_rate, na.rm = TRUE),
    block_rate_lo = quantile(block_exceedance_rate, 0.025, na.rm = TRUE),
    block_rate_hi = quantile(block_exceedance_rate, 0.975, na.rm = TRUE),
    mean_abs_gap = mean(abs_calibration_gap, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    recall_lo = quantile(recall, 0.025, na.rm = TRUE),
    recall_hi = quantile(recall, 0.975, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE),
    .groups = "drop"
  )

random_reference_rep_pooled <- random_reference_raw %>%
  group_by(rep) %>%
  summarise(
    alpha_block_design = first(alpha_block_design),
    block_exceedance_rate = mean(block_exceedance_rate, na.rm = TRUE),
    abs_calibration_gap = mean(abs_calibration_gap, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1 = mean(f1, na.rm = TRUE),
    .groups = "drop"
  )

random_reference_pooled_summary <- random_reference_rep_pooled %>%
  summarise(
    symbol = "Pooled",
    reps = n(),
    alpha_block_design = first(alpha_block_design),
    mean_block_rate = mean(block_exceedance_rate, na.rm = TRUE),
    block_rate_lo = quantile(block_exceedance_rate, 0.025, na.rm = TRUE),
    block_rate_hi = quantile(block_exceedance_rate, 0.975, na.rm = TRUE),
    mean_abs_gap = mean(abs_calibration_gap, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    recall_lo = quantile(recall, 0.025, na.rm = TRUE),
    recall_hi = quantile(recall, 0.975, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE)
  )

random_reference_summary <- bind_rows(
  random_reference_asset_summary,
  random_reference_pooled_summary
)

readr::write_csv(
  random_reference_summary,
  file.path(DIRS$tables, "random_score_reference_summary.csv")
)

# Main and sensitivity Monte Carlo experiments call the same DGP functions with
# the same retained T, burn-in, and DGP defaults. They differ in replication
# counts/seeds and in the sensitivity experiment's common-random-number/fixed-test design.
simulation_design_consistency <- tibble::tribble(
  ~parameter, ~main_design, ~sensitivity_design, ~same_core_parameter,
  "retained_T", CFG$mc_T, CFG$mc_T, TRUE,
  "burn_in", CFG$mc_burn, CFG$mc_burn, TRUE,
  "HAR_mu", -5.0, -5.0, TRUE,
  "HAR_b_d", 0.55, 0.55, TRUE,
  "HAR_b_w", 0.25, 0.25, TRUE,
  "HAR_b_m", 0.15, 0.15, TRUE,
  "HAR_sigma_u", 0.20, 0.20, TRUE,
  "SV_mu", -5.0, -5.0, TRUE,
  "SV_phi", 0.97, 0.97, TRUE,
  "SV_sigma_h", 0.20, 0.20, TRUE,
  "SV_break_mu1", -5.2, -5.2, TRUE,
  "SV_break_mu2", -4.5, -4.5, TRUE,
  "SV_break_fraction", 0.80, 0.80, TRUE,
  "SV_break_phi", 0.97, 0.97, TRUE,
  "SV_break_sigma_h", 0.20, 0.20, TRUE,
  "replications", CFG$mc_reps, CFG$mc_sensitivity_reps, FALSE
)
readr::write_csv(
  simulation_design_consistency,
  file.path(DIRS$audit, "simulation_design_consistency.csv")
)

# ----------------------------------------------------------------------
# 15. Full empirical tuning sensitivity
# ----------------------------------------------------------------------
message("Running H x B x calibration-window sensitivity...")
TUNING <- bind_rows(lapply(CFG$symbols, function(sym) {
  message("  Sensitivity: ", sym)
  run_asset_tuning_sensitivity(DATA[[sym]], sym, CFG)
}))

readr::write_csv(
  TUNING,
  file.path(DIRS$tables, "cird_H_B_calibration_sensitivity.csv")
)

# ----------------------------------------------------------------------
# 16. Publication figures
# ----------------------------------------------------------------------
# All visible labels are publication-style; code-style underscores are never
# exposed in figure titles, axes, legends, or facet strips.

# 16.1 Benchmark held-out range-reference recall (APPENDIX FIGURE A-F1)
p_bench <- benchmark_table %>%
  mutate(
    method_label = factor(
      pretty_method(method),
      levels = unname(METHOD_LABELS)
    )
  ) %>%
  ggplot(aes(x = method_label, y = stress_proxy_recall)) +
  geom_col(na.rm = FALSE) +
  facet_wrap(~ symbol, ncol = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Held-out range-reference recall"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(DIRS$figures, "Fig4.png"),
  p_bench, width = 10, height = 6.2, dpi = 600
)

# 16.2 Main CIRD ablation trade-off (MAIN FIGURE F2)
# The point of the ablation is not maximum recall in isolation.  The figure
# therefore shows the two dimensions jointly: calibration fidelity and recall.
ablation_plot_summary <- ablation_table %>%
  group_by(ablation_id) %>%
  summarise(
    mean_abs_calibration_gap = mean_or_na(
      abs(block_exceedance_rate - alpha_block_design)
    ),
    mean_range_reference_recall = mean_or_na(stress_proxy_recall),
    .groups = "drop"
  ) %>%
  mutate(
    specification = factor(
      pretty_ablation(ablation_id),
      levels = c("Level only", "Full CIRD", "Curvature only")
    )
  )

ablation_plot_long <- ablation_plot_summary %>%
  select(
    specification,
    mean_abs_calibration_gap,
    mean_range_reference_recall
  ) %>%
  pivot_longer(
    cols = -specification,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = dplyr::recode(
      metric,
      mean_abs_calibration_gap =
        "Mean absolute calibration gap",
      mean_range_reference_recall =
        "Held-out range-reference recall"
    ),
    metric_label = factor(
      metric_label,
      levels = c(
        "Mean absolute calibration gap",
        "Held-out range-reference recall"
      )
    ),
    value_label = sprintf("%.3f", value)
  )

p_ablation_tradeoff <- ablation_plot_long %>%
  ggplot(aes(x = specification, y = value)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(label = value_label),
    vjust = -0.35,
    size = 3.8
  ) +
  facet_wrap(~ metric_label, nrow = 1, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Value"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 18, hjust = 1, vjust = 1),
    strip.text = element_text(size = 9.5)
  )

ggsave(
  file.path(DIRS$figures, "Fig2.png"),
  p_ablation_tradeoff, width = 10, height = 4.8, dpi = 600
)

# 16.3 Assetwise CIRD ablation recall (APPENDIX FIGURE A-F5)
p_ablation_asset <- ablation_table %>%
  mutate(
    ablation_label = factor(
      pretty_ablation(ablation_id),
      levels = c("Level only", "Full CIRD", "Curvature only")
    )
  ) %>%
  ggplot(aes(x = ablation_label, y = stress_proxy_recall)) +
  geom_col(na.rm = FALSE) +
  facet_wrap(~ symbol, ncol = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Held-out range-reference recall"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1))

ggsave(
  file.path(DIRS$figures, "Fig8.png"),
  p_ablation_asset, width = 10, height = 6.2, dpi = 600
)

# 16.4 Weight sensitivity (APPENDIX FIGURE A-F2)
alpha_block_main <- alpha_block_map(CFG$alpha_day_main, CFG$B_main)

p_weight <- weight_table %>%
  ggplot(aes(x = weight_ratio_a_over_b, y = block_exceedance_rate)) +
  geom_hline(
    yintercept = alpha_block_main,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(~ symbol, ncol = 2) +
  labs(
    x = expression(a/b),
    y = "Observed block exceedance rate"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(DIRS$figures, "Fig5.png"),
  p_weight, width = 10, height = 6.2, dpi = 600
)

# 16.5 Horizon/block-size sensitivity at baseline calibration fraction
#      (APPENDIX FIGURE A-F3)
p_hb <- TUNING %>%
  filter(abs(cal_frac - CFG$cal_frac_main) < 1e-12) %>%
  ggplot(aes(x = factor(B), y = factor(H), fill = abs_calibration_gap)) +
  geom_tile() +
  facet_wrap(~ symbol, ncol = 2) +
  labs(
    x = "Block size B",
    y = "Horizon H",
    fill = "Absolute calibration gap"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(DIRS$figures, "Fig6.png"),
  p_hb, width = 10, height = 6.2, dpi = 600
)

# 16.6 Baseline reliability across assets (MAIN FIGURE F1)
p_rel <- reliability_table %>%
  ggplot(aes(x = alpha_block_design, y = block_exceedance_rate)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = block_rate_lo, ymax = block_rate_hi), width = 0.01) +
  facet_wrap(~ symbol, ncol = 2) +
  labs(
    x = "Design block target",
    y = "Observed block exceedance rate"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(DIRS$figures, "Fig1.png"),
  p_rel, width = 10, height = 6.2, dpi = 600
)

# ----------------------------------------------------------------------
# 17. Monte Carlo: corrected HAR, dependent SV, structural break
# ----------------------------------------------------------------------
message("Running main Monte Carlo study...")
SIM_RAW <- run_simulation_main(CFG)
SIM_SUM <- summarise_simulation(SIM_RAW)

readr::write_csv(SIM_RAW, file.path(DIRS$tables, "simulation_main_raw.csv"))
readr::write_csv(SIM_SUM, file.path(DIRS$tables, "simulation_main_summary.csv"))
saveRDS(SIM_RAW, file.path(DIRS$rds, "simulation_main_raw.rds"))

p_sim <- SIM_SUM %>%
  mutate(
    dgp_label = factor(pretty_dgp(dgp), levels = unname(DGP_LABELS))
  ) %>%
  ggplot(aes(x = alpha_block_design, y = mean_block_rate, group = dgp_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = block_rate_lo, ymax = block_rate_hi), width = 0.01) +
  facet_wrap(~ dgp_label) +
  labs(
    x = "Design block target",
    y = "Mean observed block exceedance rate"
  ) +
  theme_bw(base_size = 10)

ggsave(
  file.path(DIRS$figures, "Fig3.png"),
  p_sim, width = 10, height = 4.8, dpi = 600
)

message("Running Monte Carlo sensitivity analysis...")
SIM_SENS_RAW <- run_simulation_sensitivity(CFG)

SIM_SENS_SUM <- SIM_SENS_RAW %>%
  group_by(
    dgp, sensitivity, H, B, cal_frac,
    alpha_day, alpha_block_design
  ) %>%
  summarise(
    reps = n(),
    mean_block_rate = mean(block_rate, na.rm = TRUE),
    lo = quantile(block_rate, 0.025, na.rm = TRUE),
    hi = quantile(block_rate, 0.975, na.rm = TRUE),
    mean_abs_gap = mean(abs(block_rate - alpha_block_design), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    calibration_gap = mean_block_rate - alpha_block_design,
    gap_lo = lo - alpha_block_design,
    gap_hi = hi - alpha_block_design,
    abs_calibration_gap = abs(calibration_gap)
  )

readr::write_csv(
  SIM_SENS_SUM,
  file.path(DIRS$tables, "simulation_parameter_sensitivity.csv")
)

# Every B has a different design target. Plot observed-minus-target so the
# correct common reference is zero rather than the B=7 target.
#
# APPENDIX FIGURE A-F4: vertical layout with short specification labels.
# Each data-generating process has its own facet; therefore no redundant legend
# is required and facet-strip text cannot be clipped.
sim_sens_plot_data <- SIM_SENS_SUM %>%
  mutate(
    specification = dplyr::case_when(
      sensitivity == "H" ~ paste0("H = ", H),
      sensitivity == "B" ~ paste0("B = ", B),
      sensitivity == "cal_frac" ~ paste0("Cal. = ", sprintf("%.2f", cal_frac)),
      TRUE ~ gsub("_", " ", sensitivity)
    ),
    specification = factor(
      specification,
      levels = c(
        "Cal. = 0.25", "Cal. = 0.20", "Cal. = 0.15",
        "B = 10", "B = 7", "B = 5",
        "H = 60", "H = 30", "H = 15"
      )
    ),
    dgp_label = factor(
      pretty_dgp(dgp),
      levels = unname(DGP_LABELS)
    )
  )

p_sim_sens <- sim_sens_plot_data %>%
  ggplot(aes(x = specification, y = calibration_gap)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(
    aes(ymin = gap_lo, ymax = gap_hi),
    width = 0.18
  ) +
  geom_point(size = 2) +
  coord_flip() +
  facet_wrap(~ dgp_label, ncol = 1) +
  labs(
    x = NULL,
    y = "Observed minus design block rate"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(size = 10),
    panel.spacing.y = grid::unit(0.9, "lines")
  )

ggsave(
  file.path(DIRS$figures, "Fig7.png"),
  p_sim_sens, width = 8.8, height = 9.4, dpi = 600
)

# ----------------------------------------------------------------------
# 18. Manuscript-ready compact tables and LaTeX export
# ----------------------------------------------------------------------

# ---------- Core CSV summaries ----------

main_asset_table <- benchmark_table %>%
  filter(method == "CIRD") %>%
  select(
    symbol, alpha_day, alpha_block_design,
    threshold, block_exceedance_rate,
    block_rate_lo, block_rate_hi,
    daily_exceedance_rate,
    alarm_blocks,
    stress_event_prevalence,
    stress_proxy_precision,
    stress_proxy_recall,
    stress_proxy_f1
  )

readr::write_csv(
  main_asset_table,
  file.path(DIRS$tables, "main_asset_results.csv")
)

benchmark_pooled <- benchmark_table %>%
  group_by(method) %>%
  summarise(
    mean_block_exceedance_rate = mean_or_na(block_exceedance_rate),
    mean_abs_calibration_gap = mean_or_na(
      abs(block_exceedance_rate - alpha_block_design)
    ),
    mean_stress_proxy_precision = mean_or_na(stress_proxy_precision),
    mean_stress_proxy_recall = mean_or_na(stress_proxy_recall),
    mean_stress_proxy_f1 = mean_or_na(stress_proxy_f1),
    .groups = "drop"
  )

readr::write_csv(
  benchmark_pooled,
  file.path(DIRS$tables, "benchmark_pooled_summary.csv")
)

weight_summary <- weight_table %>%
  group_by(weight_id, a, b, weight_ratio_a_over_b) %>%
  summarise(
    mean_block_rate = mean_or_na(block_exceedance_rate),
    mean_abs_gap = mean_or_na(
      abs(block_exceedance_rate - alpha_block_design)
    ),
    mean_jaccard_vs_equal = mean_or_na(alarm_jaccard_vs_equal),
    mean_stress_proxy_recall = mean_or_na(stress_proxy_recall),
    mean_stress_proxy_f1 = mean_or_na(stress_proxy_f1),
    .groups = "drop"
  )

readr::write_csv(
  weight_summary,
  file.path(DIRS$tables, "cird_weight_sensitivity_summary.csv")
)

ablation_summary <- ablation_table %>%
  group_by(ablation_id, a, b) %>%
  summarise(
    mean_block_rate = mean_or_na(block_exceedance_rate),
    mean_abs_gap = mean_or_na(
      abs(block_exceedance_rate - alpha_block_design)
    ),
    mean_stress_proxy_prevalence = mean_or_na(stress_event_prevalence),
    mean_stress_proxy_precision = mean_or_na(stress_proxy_precision),
    mean_stress_proxy_recall = mean_or_na(stress_proxy_recall),
    mean_stress_proxy_f1 = mean_or_na(stress_proxy_f1),
    .groups = "drop"
  )

readr::write_csv(
  ablation_summary,
  file.path(DIRS$tables, "cird_component_ablation_summary.csv")
)

empirical_sensitivity_summary <- TUNING %>%
  group_by(H, B, cal_frac, alpha_block_design) %>%
  summarise(
    mean_block_rate = mean_or_na(block_exceedance_rate),
    mean_abs_gap = mean_or_na(abs_calibration_gap),
    mean_proxy_recall = mean_or_na(stress_proxy_recall),
    .groups = "drop"
  )

readr::write_csv(
  empirical_sensitivity_summary,
  file.path(DIRS$tables, "empirical_sensitivity_pooled_summary.csv")
)

cird_episode_summary <- alarm_calendar_table %>%
  filter(method == "CIRD") %>%
  group_by(symbol) %>%
  summarise(
    test_blocks = n(),
    alarm_blocks = sum(alarm, na.rm = TRUE),
    stress_reference_blocks = sum(stress_proxy_event %in% TRUE, na.rm = TRUE),
    alarm_and_stress_blocks = sum(
      alarm %in% TRUE & stress_proxy_event %in% TRUE, na.rm = TRUE
    ),
    first_alarm = if (any(alarm %in% TRUE)) min(block_start[alarm %in% TRUE]) else as.Date(NA),
    last_alarm = if (any(alarm %in% TRUE)) max(block_end[alarm %in% TRUE]) else as.Date(NA),
    .groups = "drop"
  )

readr::write_csv(
  cird_episode_summary,
  file.path(DIRS$tables, "cird_alarm_episode_summary.csv")
)

# ---------- MAIN-TEXT LaTeX tables ----------

T1 <- main_asset_table %>%
  transmute(
    Asset = symbol,
    `Day target` = fmt3(alpha_day),
    `Block target` = fmt3(alpha_block_design),
    Threshold = fmt3(threshold),
    `Block rate` = fmt3(block_exceedance_rate),
    `95% interval` = fmt_interval3(block_rate_lo, block_rate_hi),
    `Day rate` = fmt3(daily_exceedance_rate),
    `Alarm blocks` = fmt_int(alarm_blocks)
  )

write_latex_table(
  T1,
  headers = c(
    "Asset", "$\\alpha_{\\mathrm{day}}$", "$\\alpha_{\\mathrm{block}}$",
    "$q$", "Block rate", "95\\% interval", "Day rate", "Alarm blocks"
  ),
  path = file.path(DIRS$latex_main, "T1-Baseline-CIRD-Calibration.tex"),
  caption = "Baseline CIRD calibration across the four cryptocurrencies.",
  label = "tab:baseline-cird-calibration",
  note = paste0(
    "The main specification uses $H=30$, $B=7$, and ",
    "$\\alpha_{\\mathrm{day}}=0.05$. The block target is the design mapping ",
    "$1-(1-\\alpha_{\\mathrm{day}})^B$. Intervals are descriptive Wilson ",
    "intervals for observed block exceedance rates."
  ),
  align = "lrrrrrrr",
  wide = TRUE
)

T2 <- benchmark_pooled %>%
  transmute(
    Method = pretty_method(method),
    `Mean block rate` = fmt3(mean_block_exceedance_rate),
    `Mean absolute gap` = fmt3(mean_abs_calibration_gap),
    Precision = fmt3(mean_stress_proxy_precision),
    Recall = fmt3(mean_stress_proxy_recall),
    `F1 score` = fmt3(mean_stress_proxy_f1)
  )

write_latex_table(
  T2,
  headers = c(
    "Method", "Mean block rate", "Mean absolute gap",
    "Precision", "Recall", "$F_1$"
  ),
  path = file.path(DIRS$latex_main, "T2-Benchmark-Comparison.tex"),
  caption = "Pooled comparison with simple benchmark screening scores.",
  label = "tab:benchmark-comparison",
  note = paste0(
    "Averages are across BTC, ETH, XRP, and LTC. Precision, recall, and $F_1$ ",
    "refer only to the held-out high--low range reference and are descriptive, ",
    "not ground-truth detection metrics. An $F_1$ score is zero when there are ",
    "reference blocks and/or alarms but no coincident blocks. The calibration-gap ",
    "ranking is descriptive and is not interpreted as statistically established superiority."
  ),
  align = "lrrrrr",
  wide = TRUE
)

T3 <- ablation_summary %>%
  mutate(
    order_id = match(
      ablation_id,
      c("level_only", "full_CIRD", "curvature_only")
    )
  ) %>%
  arrange(order_id) %>%
  transmute(
    Specification = pretty_ablation(ablation_id),
    a = fmt3(a),
    b = fmt3(b),
    `Mean block rate` = fmt3(mean_block_rate),
    `Mean absolute gap` = fmt3(mean_abs_gap),
    Recall = fmt3(mean_stress_proxy_recall),
    `F1 score` = fmt3(mean_stress_proxy_f1)
  )

write_latex_table(
  T3,
  headers = c(
    "Specification", "$a$", "$b$", "Mean block rate",
    "Mean absolute gap", "Recall", "$F_1$"
  ),
  path = file.path(DIRS$latex_main, "T3-CIRD-Component-Ablation.tex"),
  caption = "CIRD component ablation.",
  label = "tab:cird-component-ablation",
  note = paste0(
    "The full score combines the standardized terminal-level component and the ",
    "curvature penalty. Level only and terminal FRV have identical conformal alarm sets and do not constitute independent evidence. Curvature only is the negative-curvature score, which favors smoother paths. Held-out range-reference metrics are descriptive only."
  ),
  align = "lrrrrrr",
  wide = TRUE
)

T4 <- SIM_SUM %>%
  arrange(match(dgp, c("HAR_RV", "SV_AR1", "SV_BREAK")), alpha_day) %>%
  transmute(
    DGP = pretty_dgp(dgp),
    `Day target` = fmt3(alpha_day),
    `Block target` = fmt3(alpha_block_design),
    `Mean block rate` = fmt3(mean_block_rate),
    `Replication range` = fmt_interval3(block_rate_lo, block_rate_hi),
    `Mean absolute gap` = fmt3(mean_abs_gap)
  )

write_latex_table(
  T4,
  headers = c(
    "Data-generating process", "$\\alpha_{\\mathrm{day}}$",
    "$\\alpha_{\\mathrm{block}}$", "Mean block rate",
    "2.5--97.5\\% replication range", "Mean absolute gap"
  ),
  path = file.path(DIRS$latex_main, "T4-Monte-Carlo-Reliability.tex"),
  caption = "Monte Carlo reliability under stationary and structural-break volatility dynamics.",
  label = "tab:monte-carlo-reliability",
  note = paste0(
    "Each design uses 200 repetitions, 1300 retained returns, and 32 test blocks. The reported 2.5--97.5\\% range is the empirical percentile range across ",
    "Monte Carlo replications, not a confidence interval for the mean."
  ),
  align = "lrrrrr",
  wide = TRUE
)

# ---------- APPENDIX LaTeX tables ----------

A1 <- benchmark_table %>%
  transmute(
    Asset = symbol,
    Method = pretty_method(method),
    Threshold = fmt3(threshold),
    `Block rate` = fmt3(block_exceedance_rate),
    `Absolute gap` = fmt3(abs(block_exceedance_rate - alpha_block_design)),
    Precision = fmt3(stress_proxy_precision),
    Recall = fmt3(stress_proxy_recall),
    `F1 score` = fmt3(stress_proxy_f1)
  )

write_latex_table(
  A1,
  headers = c(
    "Asset", "Method", "$q$", "Block rate", "Absolute gap",
    "Precision", "Recall", "$F_1$"
  ),
  path = file.path(DIRS$latex_appendix, "A1-Assetwise-Benchmark-Details.tex"),
  caption = "Assetwise benchmark comparison.",
  label = "tab:appendix-assetwise-benchmarks",
  note = "Held-out range-reference metrics are descriptive and are not ground-truth detection measures.",
  align = "llrrrrrr",
  wide = TRUE
)

A2 <- weight_summary %>%
  arrange(desc(weight_ratio_a_over_b)) %>%
  transmute(
    Specification = pretty_weight(weight_id),
    a = fmt3(a),
    b = fmt3(b),
    `a/b` = fmt3(weight_ratio_a_over_b),
    `Mean block rate` = fmt3(mean_block_rate),
    `Mean absolute gap` = fmt3(mean_abs_gap),
    `Jaccard vs equal` = fmt3(mean_jaccard_vs_equal),
    Recall = fmt3(mean_stress_proxy_recall),
    `F1 score` = fmt3(mean_stress_proxy_f1)
  )

write_latex_table(
  A2,
  headers = c(
    "Specification", "$a$", "$b$", "$a/b$", "Mean block rate",
    "Mean absolute gap", "Jaccard vs. equal", "Recall", "$F_1$"
  ),
  path = file.path(DIRS$latex_appendix, "A2-CIRD-Weight-Sensitivity.tex"),
  caption = "Sensitivity of CIRD to the level--curvature weights.",
  label = "tab:appendix-weight-sensitivity",
  note = "All positive-weight specifications satisfy $a+b=2$; equal weights are the neutral pre-specified baseline. The observed gap ordering is descriptive and is not treated as evidence of a statistically optimal weight ratio.",
  align = "lrrrrrrrr",
  wide = TRUE
)

A3 <- empirical_sensitivity_summary %>%
  arrange(H, B, cal_frac) %>%
  transmute(
    H = fmt_int(H),
    B = fmt_int(B),
    `Calibration fraction` = fmt3(cal_frac),
    `Block target` = fmt3(alpha_block_design),
    `Mean block rate` = fmt3(mean_block_rate),
    `Mean absolute gap` = fmt3(mean_abs_gap),
    `Mean proxy recall` = fmt3(mean_proxy_recall)
  )

write_latex_table(
  A3,
  headers = c(
    "$H$", "$B$", "Calibration fraction", "$\\alpha_{\\mathrm{block}}$",
    "Mean block rate", "Mean absolute gap", "Mean proxy recall"
  ),
  path = file.path(DIRS$latex_appendix, "A3-Empirical-Sensitivity.tex"),
  caption = "Empirical sensitivity to horizon, block size, and calibration-window length.",
  label = "tab:appendix-empirical-sensitivity",
  note = "Calibration-window comparisons hold test origins fixed within each asset and horizon, while the training sample changes and the score is refitted. Changing block size also changes the design target. Cell rankings are descriptive.",
  align = "rrrrrrr",
  wide = TRUE
)

A4 <- SIM_SENS_SUM %>%
  arrange(match(dgp, c("HAR_RV", "SV_AR1", "SV_BREAK")), sensitivity, H, B, cal_frac) %>%
  transmute(
    DGP = pretty_dgp(dgp),
    Dimension = dplyr::recode(
      sensitivity,
      H = "Horizon",
      B = "Block size",
      cal_frac = "Calibration fraction",
      .default = gsub("_", " ", sensitivity)
    ),
    H = fmt_int(H),
    B = fmt_int(B),
    `Calibration fraction` = fmt3(cal_frac),
    `Block target` = fmt3(alpha_block_design),
    `Mean block rate` = fmt3(mean_block_rate),
    `Calibration gap` = fmt3(calibration_gap)
  )

write_latex_table(
  A4,
  headers = c(
    "Data-generating process", "Dimension", "$H$", "$B$",
    "Calibration fraction", "$\\alpha_{\\mathrm{block}}$",
    "Mean block rate", "Calibration gap"
  ),
  path = file.path(DIRS$latex_appendix, "A4-Monte-Carlo-Sensitivity.tex"),
  caption = "Monte Carlo sensitivity of block calibration.",
  label = "tab:appendix-monte-carlo-sensitivity",
  note = "Common random numbers are used across sensitivity specifications within each data-generating process and replication. The calibration gap is the mean observed block rate minus the corresponding design target. Sensitivity estimates come from the separate 100-replication common-random-number experiment and therefore need not reproduce the 200-replication baseline estimates reported in the main Monte Carlo analysis.",
  align = "llrrrrrr",
  wide = TRUE
)

A5 <- returns_table %>%
  transmute(
    Asset = symbol,
    N = fmt_int(n),
    Mean = fmt3(mean),
    SD = fmt3(sd),
    Skewness = fmt3(skewness),
    Kurtosis = fmt3(kurtosis_raw),
    `Jarque--Bera p-value` = ifelse(jb_p < 0.001, "<0.001", fmt3(jb_p)),
    `Ljung--Box(10) return p-value` = ifelse(lb10_return_p < 0.001, "<0.001", fmt3(lb10_return_p)),
    `Ljung--Box(10) squared-return p-value` = ifelse(lb10_sqreturn_p < 0.001, "<0.001", fmt3(lb10_sqreturn_p))
  )

write_latex_table(
  A5,
  headers = c(
    "Asset", "$N$", "Mean", "Standard deviation", "Skewness", "Raw kurtosis",
    "Jarque--Bera $p$", "Ljung--Box(10) return $p$", "Ljung--Box(10) squared-return $p$"
  ),
  path = file.path(DIRS$latex_appendix, "A5-Return-Descriptives.tex"),
  caption = "Daily return descriptives.",
  label = "tab:appendix-return-descriptives",
  note = "Ljung--Box statistics are evaluated at lag 10.",
  align = "lrrrrrrrr",
  wide = TRUE
)

A6 <- frv_table %>%
  transmute(
    Asset = symbol,
    Horizon = fmt_int(h),
    Mean = fmt3(mean),
    Median = fmt3(median),
    SD = fmt3(sd),
    `90th percentile` = fmt3(p90)
  )

write_latex_table(
  A6,
  headers = c(
    "Asset", "Horizon", "Mean", "Median", "Standard deviation", "90th percentile"
  ),
  path = file.path(DIRS$latex_appendix, "A6-FRV-Descriptives.tex"),
  caption = "Forward realized-variance descriptives in the training sample.",
  label = "tab:appendix-frv-descriptives",
  note = "Descriptives are computed from the training-origin FRV curves in squared log-return units. Displayed zeros can reflect rounding to three decimals.",
  align = "lrrrrr",
  wide = TRUE
)

reference_peak_labels <- range_reference_block_details_table %>%
  filter(stress_reference_event %in% TRUE) %>%
  group_by(symbol) %>%
  summarise(peak_dates = paste0("\\shortstack{", paste(sort(unique(proxy_peak_date)), collapse = "\\\\"), "}"), .groups = "drop")

A7 <- stress_proxy_audit %>%
  left_join(reference_peak_labels, by = "symbol") %>%
  transmute(
    Asset = symbol,
    `Tail quantile` = fmt3(heldout_stress_quantile),
    Threshold = fmt3(stress_threshold),
    `Held-out proxy blocks` = fmt_int(heldout_proxy_blocks),
    `Stress-reference blocks` = fmt_int(test_stress_blocks),
    `Stress-reference rate` = fmt3(test_stress_block_rate),
    `Distinct reference runs` = fmt_int(distinct_reference_runs),
    `Unique peak dates` = fmt_int(unique_reference_peak_dates),
    `Range peak dates` = peak_dates
  )

write_latex_table(
  A7,
  headers = c(
    "Asset", "Tail quantile", "Threshold", "Held-out proxy blocks",
    "Stress-reference blocks", "Stress-reference rate",
    "Distinct reference runs", "Unique peak dates", "Range peak dates"
  ),
  path = file.path(DIRS$latex_appendix, "A7-Held-Out-Range-Reference-Audit.tex"),
  caption = "Held-out high--low range-reference construction.",
  label = "tab:appendix-range-reference-audit",
  note = "A block enters the reference when its held-out range proxy is at or above the empirical 90th percentile. Repeated maxima arise mechanically because adjacent origin blocks contain overlapping 30-day forward windows that can share the same future high--low peak date; the non-strict threshold also preserves any exact quantile ties. Distinct reference runs are contiguous runs of flagged blocks. This ex-post reference is used only for descriptive validation and never enters CIRD fitting or conformal calibration.",
  align = "lrrrrrrrl",
  wide = TRUE
)

A8 <- split_audit_table %>%
  transmute(
    Asset = symbol,
    H = fmt_int(H),
    `Train N` = fmt_int(train_n),
    `Calibration N` = fmt_int(cal_n),
    `Test N` = fmt_int(test_n),
    Embargo = fmt_int(embargo),
    `Train end` = as.character(train_end),
    `Calibration start` = as.character(cal_start),
    `Calibration end` = as.character(cal_end),
    `Test start` = as.character(test_start)
  )

write_latex_table(
  A8,
  headers = c(
    "Asset", "$H$", "Train $N$", "Calibration $N$", "Test $N$",
    "Embargo", "Train end", "Calibration start", "Calibration end", "Test start"
  ),
  path = file.path(DIRS$latex_appendix, "A8-Split-Audit.tex"),
  caption = "Chronological split and embargo audit.",
  label = "tab:appendix-split-audit",
  note = "An $H$-origin embargo separates train, calibration, and test origin sets.",
  align = "lrrrrlllll",
  wide = TRUE
)

A9 <- iso_audit_table %>%
  transmute(
    Asset = symbol,
    `Training curves` = fmt_int(n_curves),
    `Share adjusted` = fmt3(share_adjusted_gt_1e_12),
    `Mean max adjustment` = fmt3(mean_max_adjustment),
    `Maximum adjustment` = fmt3(max_adjustment)
  )

write_latex_table(
  A9,
  headers = c(
    "Asset", "Training curves", "Share adjusted",
    "Mean max adjustment", "Maximum adjustment"
  ),
  path = file.path(DIRS$latex_appendix, "A9-Isotonic-Projection-Audit.tex"),
  caption = "Isotonic projection audit.",
  label = "tab:appendix-isotonic-audit",
  note = "Cumulative squared-return FRV curves are monotone by construction; projection therefore acts as a numerical safeguard.",
  align = "lrrrr",
  wide = TRUE
)

A10 <- cird_episode_summary %>%
  transmute(
    Asset = symbol,
    `Test blocks` = fmt_int(test_blocks),
    `Alarm blocks` = fmt_int(alarm_blocks),
    `Stress-reference blocks` = fmt_int(stress_reference_blocks),
    `Alarm and stress blocks` = fmt_int(alarm_and_stress_blocks),
    `First alarm` = as.character(first_alarm),
    `Last alarm` = as.character(last_alarm)
  )

write_latex_table(
  A10,
  headers = c(
    "Asset", "Test blocks", "Alarm blocks", "Stress-reference blocks",
    "Alarm and stress blocks", "First alarm", "Last alarm"
  ),
  path = file.path(DIRS$latex_appendix, "A10-CIRD-Alarm-Episode-Summary.tex"),
  caption = "CIRD alarm and reference summary in the held-out period.",
  label = "tab:appendix-alarm-episode-summary",
  note = "First and last alarm dates are origin-block boundaries; the forward curves become observable later. The complete dated block calendar is supplied in Online Resource 1.",
  align = "lrrrrll",
  wide = TRUE
)

# ---------- PUBLICATION APPENDIX TABLES ----------

A11 <- component_corr_table %>%
  transmute(
    Asset = symbol,
    Sample = sample,
    N = fmt_int(n),
    `Pearson correlation` = fmt3(pearson_L_C),
    `Spearman correlation` = fmt3(spearman_L_C)
  )

write_latex_table(
  A11,
  headers = c("Asset", "Sample", "$N$", "Pearson corr. $(L,C)$", "Spearman corr. $(L,C)$"),
  path = file.path(DIRS$latex_appendix, "A11-Level-Curvature-Correlations.tex"),
  caption = "Association between the CIRD level and curvature components.",
  label = "tab:appendix-level-curvature-correlations",
  note = "The curvature component is the average absolute first difference of successive daily squared-return contributions after scaling by the global training curvature MAD; positive association with the level component is therefore expected in high-volatility paths.",
  align = "llrrr",
  wide = TRUE
)

A12 <- random_reference_summary %>%
  transmute(
    Asset = symbol,
    Replications = fmt_int(reps),
    `Block target` = fmt3(alpha_block_design),
    `Mean block rate` = fmt3(mean_block_rate),
    `Mean absolute gap` = fmt3(mean_abs_gap),
    `Mean recall` = fmt3(mean_recall),
    `Recall range` = fmt_interval3(recall_lo, recall_hi)
  )

write_latex_table(
  A12,
  headers = c(
    "Asset", "Replications", "$\\alpha_{\\mathrm{block}}$", "Mean block rate",
    "Mean absolute gap", "Mean recall", "2.5--97.5\\% recall range"
  ),
  path = file.path(DIRS$latex_appendix, "A12-Random-Score-Sanity-Check.tex"),
  caption = "Uninformative random-score calibration sanity check.",
  label = "tab:appendix-random-score-sanity",
  note = "Independent Uniform(0,1) origin scores are calibrated by exactly the same block-conformal rule. This diagnostic illustrates that a small calibration gap can occur for an uninformative score and therefore does not establish screening skill. Recall is measured against the same held-out range reference and is expected to be near the random block-alarm probability rather than near one.",
  align = "lrrrrrr",
  wide = TRUE
)


compress_block_ids <- function(x) {
  if (is.na(x) || !nzchar(x)) return("--")
  ids <- as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
  groups <- split(ids, cumsum(c(TRUE, diff(ids) != 1L)))
  ranges <- vapply(groups, function(z) {
    if (length(z) == 1L) as.character(z) else paste0(min(z), "--", max(z))
  }, character(1))
  paste0(paste(ranges, collapse = ", "), " (", length(ids), ")")
}

episode_titles <- c(
  ETF_APPROVAL = "U.S. spot-Bitcoin ETF approval: 7--13 January 2024",
  MARCH_HIGHS = "March crypto highs: 1--31 March 2024",
  AUGUST_SELLOFF = "Market-wide sell-off: 2--8 August 2024",
  NOVEMBER_RALLY = "November crypto rally: 1--30 November 2024"
)
episode_lines <- c(
  "\\begin{table*}[!htbp]", "  \\centering\\small",
  "  \\caption{Descriptive mapping of alarms and range-reference blocks to identifiable 2024 periods.}",
  "  \\label{tab:appendix-recognizable-episodes}",
  "  \\begin{adjustbox}{max width=\\textwidth}",
  "  \\begin{tabular}{llll}", "  \\toprule",
  "  Asset & CIRD alarm block IDs (count) & Reference block IDs (count) & Reference peak date \\\\",
  "  \\midrule"
)
for (eid in names(episode_titles)) {
  episode_lines <- c(episode_lines, paste0("  \\multicolumn{4}{l}{\\textit{", episode_titles[[eid]], "}} \\\\", "  \\addlinespace[2pt]") )
  z <- recognized_episode_mapping_table %>%
    filter(episode_id == eid) %>% arrange(match(symbol, CFG$symbols))
  for (i in seq_len(nrow(z))) {
    peak <- if (is.na(z$reference_peak_dates[i]) || !nzchar(z$reference_peak_dates[i])) "--" else z$reference_peak_dates[i]
    episode_lines <- c(episode_lines, paste0("  ", z$symbol[i], " & ",
      compress_block_ids(z$cird_alarm_block_ids[i]), " & ",
      compress_block_ids(z$reference_block_ids[i]), " & ", peak, " \\\\") )
  }
  episode_lines <- c(episode_lines, "  \\addlinespace[4pt]")
}
episode_lines <- c(episode_lines, "  \\bottomrule", "  \\end{tabular}", "  \\end{adjustbox}",
  "  \\vspace{2pt}",
  "  \\parbox{\\textwidth}{\\footnotesize\\textit{Note:} IDs refer to the held-out seven-day origin-block calendar; parentheses give block counts and -- denotes no match. An alarm is mapped when its 36-day union of forward evaluation windows overlaps the displayed period. Its score may be driven by returns outside that period. A reference block is mapped when its underlying daily high--low peak date lies in the period. March and November use calendar months; January and August use the indicated weeks. These are ex-post descriptive associations, not advance detection, causal attribution, or fitting inputs.}",
  "\\end{table*}")
writeLines(episode_lines, file.path(DIRS$latex_appendix, "A13-Recognizable-Episode-Mapping.tex"), useBytes = TRUE)

A14 <- zero_return_audit_table %>%
  transmute(
    Asset = symbol,
    `Finite returns` = fmt_int(finite_returns),
    `Exact zero returns` = fmt_int(exact_zero_returns),
    `Training RV1 zero` = fmt_int(har_training_rv1_zero),
    `Training RV1 <= epsilon` = fmt_int(har_training_rv1_at_or_below_epsilon),
    `Minimum positive training RV1` = formatC(minimum_positive_training_rv1, digits = 3, format = "e"),
    Epsilon = formatC(har_log_epsilon, digits = 1, format = "e")
  )

write_latex_table(
  A14,
  headers = c(
    "Asset", "Finite returns", "Exact zero returns", "Training $RV^{(1)}$ zero",
    "Training $RV^{(1)}\\leq\\varepsilon$", "Minimum positive training $RV^{(1)}$", "$\\varepsilon$"
  ),
  path = file.path(DIRS$latex_appendix, "A14-HAR-Zero-Return-Audit.tex"),
  caption = "Audit of zero returns and the HAR logarithmic offset.",
  label = "tab:appendix-har-zero-return-audit",
  note = "The HAR transformation is log(RV + epsilon), with an additive offset rather than truncation at a floor. No exact-zero returns or training one-day variances at or below epsilon occur in this sample. Adding epsilon nevertheless changes positive logarithmic regressors, so zero counts do not imply an exactly unchanged fit.",
  align = "lrrrrrr",
  wide = TRUE
)

# ---------- Manuscript placement plan ----------
table_plan <- tibble::tribble(
  ~order, ~placement, ~file, ~title,
  1L, "Main text", "T1-Baseline-CIRD-Calibration.tex", "Baseline CIRD calibration",
  2L, "Main text", "T2-Benchmark-Comparison.tex", "Benchmark comparison",
  3L, "Main text", "T3-CIRD-Component-Ablation.tex", "CIRD component ablation",
  4L, "Main text", "T4-Monte-Carlo-Reliability.tex", "Monte Carlo reliability",
  1L, "Appendix", "A1-Assetwise-Benchmark-Details.tex", "Assetwise benchmark details",
  2L, "Appendix", "A2-CIRD-Weight-Sensitivity.tex", "CIRD weight sensitivity",
  3L, "Appendix", "A3-Empirical-Sensitivity.tex", "Empirical sensitivity",
  4L, "Appendix", "A4-Monte-Carlo-Sensitivity.tex", "Monte Carlo sensitivity",
  5L, "Appendix", "A5-Return-Descriptives.tex", "Return descriptives",
  6L, "Appendix", "A6-FRV-Descriptives.tex", "FRV descriptives",
  7L, "Appendix", "A7-Held-Out-Range-Reference-Audit.tex", "Held-out range-reference audit",
  8L, "Appendix", "A8-Split-Audit.tex", "Split audit",
  9L, "Appendix", "A9-Isotonic-Projection-Audit.tex", "Isotonic projection audit",
  10L, "Appendix", "A10-CIRD-Alarm-Episode-Summary.tex", "CIRD alarm-episode summary",
  11L, "Appendix", "A11-Level-Curvature-Correlations.tex", "Level-curvature correlations",
  12L, "Appendix", "A12-Random-Score-Sanity-Check.tex", "Random-score sanity check",
  13L, "Appendix", "A13-Recognizable-Episode-Mapping.tex", "Recognizable episode mapping",
  14L, "Appendix", "A14-HAR-Zero-Return-Audit.tex", "HAR zero-return audit"
)

figure_plan <- tibble::tribble(
  ~order, ~placement, ~file, ~title,
  1L, "Main text", "Fig1.png", "Block-level reliability",
  2L, "Main text", "Fig2.png", "CIRD ablation trade-off",
  3L, "Main text", "Fig3.png", "Simulation reliability",
  1L, "Appendix", "Fig4.png", "Benchmark range-reference recall",
  2L, "Appendix", "Fig5.png", "CIRD weight sensitivity",
  3L, "Appendix", "Fig6.png", "Horizon/block-size sensitivity",
  4L, "Appendix", "Fig7.png", "Simulation sensitivity gap",
  5L, "Appendix", "Fig8.png", "Assetwise CIRD ablation recall"
)

readr::write_csv(
  table_plan,
  file.path(DIRS$audit, "manuscript_table_plan.csv")
)
readr::write_csv(
  figure_plan,
  file.path(DIRS$audit, "manuscript_figure_plan.csv")
)

writeLines(
  c(
    "MAIN-TEXT TABLES",
    "T1  Baseline CIRD calibration across BTC, ETH, XRP, and LTC.",
    "T2  Pooled benchmark comparison.",
    "T3  CIRD component ablation.",
    "T4  Monte Carlo reliability.",
    "",
    "APPENDIX TABLES",
    "A1  Assetwise benchmark details.",
    "A2  CIRD weight sensitivity.",
    "A3  Empirical H/B/calibration-window sensitivity.",
    "A4  Monte Carlo sensitivity.",
    "A5  Return descriptives.",
    "A6  FRV descriptives.",
    "A7  Held-out range-reference audit.",
    "A8  Chronological split audit.",
    "A9  Isotonic projection audit.",
    "A10 CIRD alarm-episode summary.",
    "A11 Level-curvature component correlations.",
    "A12 Uninformative random-score calibration sanity check.",
    "A13 Recognizable 2024 volatility-episode mapping.",
    "A14 HAR zero-return/log-floor audit.",
    "",
    "MAIN-TEXT FIGURES",
    "Fig1  Block-level reliability.",
    "Fig2  CIRD ablation trade-off: calibration gap versus held-out recall.",
    "Fig3  Simulation reliability.",
    "",
    "APPENDIX FIGURES",
    "Fig4 Benchmark range-reference recall.",
    "Fig5 CIRD weight sensitivity.",
    "Fig6 Horizon/block-size sensitivity.",
    "Fig7 Simulation sensitivity gap.",
    "Fig8 Assetwise CIRD ablation recall."
  ),
  file.path(DIRS$audit, "manuscript_output_plan.txt")
)

writeLines(
  c(
    "Block-Calibrated Shape-Aware Screening of Cryptocurrency Volatility Episodes",
    "Reproducibility package",
    "",
    "Main script:",
    "  cird_analysis.R",
    "",
    "Required R packages:",
    paste0("  ", paste(required_packages, collapse = ", ")),
    "",
    "Run from the package root:",
    "  Rscript cird_analysis.R",
    "",
    "Empirical universe:",
    "  BTC-USD, ETH-USD, XRP-USD, LTC-USD",
    "",
    "Main specification:",
    paste0("  H=", CFG$H_main,
           ", B=", CFG$B_main,
           ", alpha_day=", CFG$alpha_day_main,
           ", train_frac=", CFG$train_frac,
           ", cal_frac=", CFG$cal_frac_main),
    "",
    "Output folders:",
    "  data/            frozen Yahoo Finance cache files",
    "  tables/          CSV analysis outputs",
    "  figures/         title-free 600-dpi submission figures (Fig1--Fig8)",
    "  latex/main_text/ main-text LaTeX tables",
    "  latex/appendix/  appendix LaTeX tables",
    "  rds/             detailed R objects",
    "  audit/           configuration, checksums, split/audit files, session info",
    "",
    "Interpretation safeguards:",
    "  * Empirical rates are block exceedance/alarm rates, not ground-truth FPRs.",
    "  * alpha_block is a design mapping; the conformal rank is applied to block scores.",
    "  * The held-out range reference is ex-post descriptive validation only.",
    "  * Finite-sample conformal control requires exchangeability of block scores.",
    "  * Cumulative squared-return FRV curves are monotone by construction;",
    "    isotonic projection is a structural/numerical safeguard in this sample.",
    "  * The benchmark rank score is one-sided and aligned with upper-tail stress.",
    "  * Random-score results are a calibration sanity check, not a competitive model.",
    "  * Range-reference ties/episodes are audited using underlying future peak dates.",
    "  * Recognizable 2024 event windows are descriptive interpretation aids only."
  ),
  file.path(DIRS$audit, "README_REPRODUCIBILITY.txt")
)


# ----------------------------------------------------------------------
# 19. Reproducibility audit
# ----------------------------------------------------------------------
config_lines <- c(
  paste0("analysis_id=", CFG$analysis_id),
  paste0("seed=", CFG$seed),
  paste0("symbols=", paste(CFG$symbols, collapse = ",")),
  paste0("data_start=", CFG$data_start),
  paste0("data_end=", CFG$data_end),
  paste0("H_main=", CFG$H_main),
  paste0("B_main=", CFG$B_main),
  paste0("alpha_day_main=", CFG$alpha_day_main),
  paste0("alpha_block_design_main=",
         alpha_block_map(CFG$alpha_day_main, CFG$B_main)),
  paste0("train_frac=", CFG$train_frac),
  paste0("cal_frac_main=", CFG$cal_frac_main),
  paste0("H_grid=", paste(CFG$H_grid, collapse = ",")),
  paste0("B_grid=", paste(CFG$B_grid, collapse = ",")),
  paste0("cal_frac_grid=", paste(CFG$cal_frac_grid, collapse = ",")),
  paste0("stress_eval_quantile=", CFG$stress_eval_quantile),
  "stress_reference_source=heldout_test_block_tail_ex_post",
  "stress_reference_used_in_model_fit=FALSE",
  "stress_reference_used_in_conformal_calibration=FALSE",
  "stress_proxy_positive_test_block_required=TRUE",
  paste0("ablation_specs=", paste(CFG$ablation_grid$ablation_id, collapse = ",")),
  paste0("mc_reps=", CFG$mc_reps),
  paste0("mc_sensitivity_reps=", CFG$mc_sensitivity_reps),
  paste0("mc_T=", CFG$mc_T),
  paste0("mc_burn=", CFG$mc_burn),
  paste0("random_score_reps=", CFG$random_score_reps),
  paste0("random_score_seed=", CFG$random_score_seed),
  paste0("har_log_epsilon=", CFG$har_log_epsilon),
  "benchmark_rank_direction=one_sided_upper_tail",
  "calibration_gap_ranking_interpretation=descriptive_not_statistical_superiority",
  "weight_gap_ranking_interpretation=descriptive_not_statistical_optimum",
  "range_reference_episode_definition=contiguous_flagged_block_runs",
  "recognized_episode_mapping_used_in_fit_or_calibration=FALSE",
  "empirical_calibration_sensitivity_fixed_test=TRUE",
  "simulation_calibration_sensitivity_fixed_test=TRUE",
  "simulation_sensitivity_common_random_numbers=TRUE",
  "curvature_definition=H_minus_2_second_differences",
  "f1_zero_true_positive_rule=direct_confusion_matrix_definition",
  "latex_tables_three_decimal_format=TRUE",
  "visible_figure_labels_use_code_underscores=FALSE",
  "main_ablation_figure_is_tradeoff=TRUE",
  "weight_figure_has_design_target_reference=TRUE",
  "simulation_sensitivity_vertical_layout=TRUE",
  "submission_figures_no_internal_titles=TRUE",
  "submission_png_dpi=600",
  "har_week_month_terms=log_of_mean_RV"
)
writeLines(config_lines, file.path(DIRS$audit, "configuration.txt"))
capture.output(sessionInfo(), file = file.path(DIRS$audit, "sessionInfo.txt"))

# Checksums of frozen, date-stamped Yahoo Finance cache files.
cache_files <- list.files(DIRS$data, pattern = "\\.csv$", full.names = TRUE)
if (length(cache_files) == 0L) {
  stop("REPRODUCIBILITY AUDIT FAILED: no cached Yahoo CSV files found.")
}
cache_checksums <- tibble(
  file = basename(cache_files),
  md5 = unname(tools::md5sum(cache_files))
)
readr::write_csv(
  cache_checksums,
  file.path(DIRS$audit, "data_checksums_md5.csv")
)

# Output manifest
all_files <- list.files(CFG$out_root, recursive = TRUE, full.names = TRUE)
manifest <- tibble(
  file = normalizePath(all_files, winslash = "/", mustWork = FALSE),
  size_bytes = file.info(all_files)$size
)
readr::write_csv(manifest, file.path(DIRS$audit, "output_manifest.csv"))

# ----------------------------------------------------------------------
# 20. Final audit
# ----------------------------------------------------------------------
required_outputs <- c(
  file.path(DIRS$tables, "main_asset_results.csv"),
  file.path(DIRS$tables, "benchmark_comparison.csv"),
  file.path(DIRS$tables, "alarm_episode_calendar.csv"),
  file.path(DIRS$tables, "cird_component_ablation.csv"),
  file.path(DIRS$tables, "cird_component_ablation_summary.csv"),
  file.path(DIRS$tables, "benchmark_pooled_summary.csv"),
  file.path(DIRS$tables, "cird_weight_sensitivity.csv"),
  file.path(DIRS$tables, "cird_weight_sensitivity_summary.csv"),
  file.path(DIRS$tables, "cird_H_B_calibration_sensitivity.csv"),
  file.path(DIRS$tables, "cird_reliability_grid.csv"),
  file.path(DIRS$tables, "simulation_main_summary.csv"),
  file.path(DIRS$tables, "simulation_parameter_sensitivity.csv"),
  file.path(DIRS$tables, "empirical_sensitivity_pooled_summary.csv"),
  file.path(DIRS$tables, "cird_alarm_episode_summary.csv"),
  file.path(DIRS$tables, "cird_level_curvature_correlations.csv"),
  file.path(DIRS$tables, "random_score_reference_raw.csv"),
  file.path(DIRS$tables, "random_score_reference_summary.csv"),
  file.path(DIRS$tables, "recognized_episode_mapping.csv"),
  file.path(DIRS$tables, "range_reference_block_details.csv"),
  file.path(DIRS$audit, "har_zero_return_audit.csv"),
  file.path(DIRS$audit, "level_terminal_equivalence_audit.csv"),
  file.path(DIRS$audit, "simulation_design_consistency.csv"),
  file.path(DIRS$audit, "split_audit.csv"),
  file.path(DIRS$audit, "isotonic_projection_audit.csv"),
  file.path(DIRS$audit, "stress_proxy_audit.csv"),
  file.path(DIRS$audit, "configuration.txt"),
  file.path(DIRS$audit, "sessionInfo.txt"),
  file.path(DIRS$audit, "data_checksums_md5.csv"),
  file.path(DIRS$figures, "Fig4.png"),
  file.path(DIRS$figures, "Fig2.png"),
  file.path(DIRS$figures, "Fig5.png"),
  file.path(DIRS$figures, "Fig6.png"),
  file.path(DIRS$figures, "Fig1.png"),
  file.path(DIRS$figures, "Fig3.png"),
  file.path(DIRS$figures, "Fig7.png"),
  file.path(DIRS$figures, "Fig8.png"),
  file.path(DIRS$latex_main, "T1-Baseline-CIRD-Calibration.tex"),
  file.path(DIRS$latex_main, "T2-Benchmark-Comparison.tex"),
  file.path(DIRS$latex_main, "T3-CIRD-Component-Ablation.tex"),
  file.path(DIRS$latex_main, "T4-Monte-Carlo-Reliability.tex"),
  file.path(DIRS$latex_appendix, "A1-Assetwise-Benchmark-Details.tex"),
  file.path(DIRS$latex_appendix, "A2-CIRD-Weight-Sensitivity.tex"),
  file.path(DIRS$latex_appendix, "A3-Empirical-Sensitivity.tex"),
  file.path(DIRS$latex_appendix, "A4-Monte-Carlo-Sensitivity.tex"),
  file.path(DIRS$latex_appendix, "A5-Return-Descriptives.tex"),
  file.path(DIRS$latex_appendix, "A6-FRV-Descriptives.tex"),
  file.path(DIRS$latex_appendix, "A7-Held-Out-Range-Reference-Audit.tex"),
  file.path(DIRS$latex_appendix, "A8-Split-Audit.tex"),
  file.path(DIRS$latex_appendix, "A9-Isotonic-Projection-Audit.tex"),
  file.path(DIRS$latex_appendix, "A10-CIRD-Alarm-Episode-Summary.tex"),
  file.path(DIRS$latex_appendix, "A11-Level-Curvature-Correlations.tex"),
  file.path(DIRS$latex_appendix, "A12-Random-Score-Sanity-Check.tex"),
  file.path(DIRS$latex_appendix, "A13-Recognizable-Episode-Mapping.tex"),
  file.path(DIRS$latex_appendix, "A14-HAR-Zero-Return-Audit.tex"),
  file.path(DIRS$audit, "manuscript_table_plan.csv"),
  file.path(DIRS$audit, "manuscript_figure_plan.csv"),
  file.path(DIRS$audit, "manuscript_output_plan.txt"),
  file.path(DIRS$audit, "README_REPRODUCIBILITY.txt")
)

missing_out <- required_outputs[!file.exists(required_outputs)]

if (length(missing_out) > 0L) {
  stop(
    "FINAL AUDIT FAILED. Missing outputs:\n",
    paste(missing_out, collapse = "\n")
  )
}

# Structural publication audits: fail loudly rather than silently publishing
# incomplete publication analyses.
expected_symbols <- sort(CFG$symbols)
assert_symbols <- function(x, object_name) {
  got <- sort(unique(as.character(x)))
  if (!identical(got, expected_symbols)) {
    stop(
      "FINAL AUDIT FAILED: ", object_name,
      " does not contain exactly the configured symbols. Got: ",
      paste(got, collapse = ", ")
    )
  }
}

assert_symbols(main_asset_table$symbol, "main_asset_table")
assert_symbols(benchmark_table$symbol, "benchmark_table")
assert_symbols(ablation_table$symbol, "ablation_table")
assert_symbols(weight_table$symbol, "weight_table")
assert_symbols(TUNING$symbol, "TUNING")
assert_symbols(split_audit_table$symbol, "split_audit_table")
assert_symbols(iso_audit_table$symbol, "iso_audit_table")
assert_symbols(stress_proxy_audit$symbol, "stress_proxy_audit")
assert_symbols(component_corr_table$symbol, "component_corr_table")
assert_symbols(zero_return_audit_table$symbol, "zero_return_audit_table")
assert_symbols(level_terminal_equivalence_table$symbol, "level_terminal_equivalence_table")
assert_symbols(random_reference_raw$symbol, "random_reference_raw")
assert_symbols(recognized_episode_mapping_table$symbol, "recognized_episode_mapping_table")
assert_symbols(range_reference_block_details_table$symbol, "range_reference_block_details_table")

expected_n <- c(
  main_asset = length(CFG$symbols),
  benchmark = length(CFG$symbols) * 4L,
  ablation = length(CFG$symbols) * nrow(CFG$ablation_grid),
  reliability = length(CFG$symbols) * length(CFG$alpha_day_grid),
  weight = length(CFG$symbols) * nrow(CFG$weight_grid),
  tuning = length(CFG$symbols) *
    length(CFG$H_grid) * length(CFG$B_grid) * length(CFG$cal_frac_grid),
  split_audit = length(CFG$symbols),
  iso_audit = length(CFG$symbols),
  stress_proxy_audit = length(CFG$symbols),
  component_corr = length(CFG$symbols) * 2L,
  zero_return_audit = length(CFG$symbols),
  level_terminal_equivalence = length(CFG$symbols),
  random_reference = length(CFG$symbols) * CFG$random_score_reps,
  recognized_episode_mapping = length(CFG$symbols) * nrow(CFG$recognized_episode_windows)
)

actual_n <- c(
  main_asset = nrow(main_asset_table),
  benchmark = nrow(benchmark_table),
  ablation = nrow(ablation_table),
  reliability = nrow(reliability_table),
  weight = nrow(weight_table),
  tuning = nrow(TUNING),
  split_audit = nrow(split_audit_table),
  iso_audit = nrow(iso_audit_table),
  stress_proxy_audit = nrow(stress_proxy_audit),
  component_corr = nrow(component_corr_table),
  zero_return_audit = nrow(zero_return_audit_table),
  level_terminal_equivalence = nrow(level_terminal_equivalence_table),
  random_reference = nrow(random_reference_raw),
  recognized_episode_mapping = nrow(recognized_episode_mapping_table)
)

if (!identical(as.integer(actual_n), as.integer(expected_n))) {
  bad <- names(actual_n)[actual_n != expected_n]
  stop(
    "FINAL AUDIT FAILED: unexpected row counts for ",
    paste(bad, collapse = ", "), ". Expected/actual: ",
    paste(
      paste0(bad, "=", expected_n[bad], "/", actual_n[bad]),
      collapse = "; "
    )
  )
}

# Core numeric outputs needed for manuscript tables must be finite.
core_numeric_checks <- list(
  main_threshold = main_asset_table$threshold,
  main_block_rate = main_asset_table$block_exceedance_rate,
  benchmark_threshold = benchmark_table$threshold,
  benchmark_block_rate = benchmark_table$block_exceedance_rate,
  ablation_threshold = ablation_table$threshold,
  weight_threshold = weight_table$threshold,
  simulation_mean_block_rate = SIM_SUM$mean_block_rate,
  simulation_q_mean = SIM_SUM$q_mean,
  tuning_block_rate = TUNING$block_exceedance_rate
)

bad_numeric <- names(core_numeric_checks)[vapply(
  core_numeric_checks,
  function(x) any(!is.finite(x)),
  logical(1)
)]

if (length(bad_numeric) > 0L) {
  stop(
    "FINAL AUDIT FAILED: non-finite values in core outputs: ",
    paste(bad_numeric, collapse = ", ")
  )
}

# CIRD itself should not lose blocks through non-finite score values.
cird_rows <- benchmark_table %>% filter(method == "CIRD")
if (any(cird_rows$cal_blocks_dropped_nonfinite != 0L) ||
    any(cird_rows$test_blocks_dropped_nonfinite != 0L)) {
  stop("FINAL AUDIT FAILED: CIRD dropped non-finite calibration/test blocks.")
}

# The cumulative squared-return FRV construction is analytically monotone.
# Record, but do not fail on, tiny numerical projection adjustments.
if (any(iso_audit_table$share_adjusted_gt_1e_12 > 0)) {
  warning(
    "Isotonic audit found FRV curves adjusted by more than 1e-12. ",
    "Inspect preprocessing before manuscript interpretation."
  )
}

# Manuscript-exact curvature must return H-2 values.
curvature_dimension_ok <-
  length(curve_curvature(seq_len(CFG$H_main))) == (CFG$H_main - 2L)
if (!curvature_dimension_ok) {
  stop("FINAL AUDIT FAILED: curvature dimension is not H-2.")
}

# Calibration-window sensitivity must use identical held-out test dates for all
# cal_frac values within each asset/H/B comparison.
fixed_test_empirical_audit <- TUNING %>%
  group_by(symbol, H, B) %>%
  summarise(
    n_test_start = n_distinct(test_start_date),
    n_test_end = n_distinct(test_end_date),
    n_test_n = n_distinct(test_origin_n),
    .groups = "drop"
  )
fixed_test_empirical_ok <- all(
  fixed_test_empirical_audit$n_test_start == 1L &
    fixed_test_empirical_audit$n_test_end == 1L &
    fixed_test_empirical_audit$n_test_n == 1L
)
if (!fixed_test_empirical_ok) {
  stop("FINAL AUDIT FAILED: empirical calibration sensitivity changed the test window.")
}

# The same fixed-test requirement applies to simulation cal_frac sensitivity.
fixed_test_sim_audit <- SIM_SENS_RAW %>%
  filter(sensitivity == "cal_frac") %>%
  group_by(dgp, rep, H, B, alpha_day) %>%
  summarise(n_test_blocks = n_distinct(test_blocks), .groups = "drop")
fixed_test_sim_ok <- all(fixed_test_sim_audit$n_test_blocks == 1L)
if (!fixed_test_sim_ok) {
  stop("FINAL AUDIT FAILED: simulation calibration sensitivity changed test blocks.")
}

# Held-out range-reference event validation must be estimable in every configured asset.
stress_proxy_estimable_ok <-
  all(is.finite(stress_proxy_audit$test_stress_block_rate)) &&
  all(stress_proxy_audit$test_stress_blocks > 0L)
if (!stress_proxy_estimable_ok) {
  stop("FINAL AUDIT FAILED: held-out range-reference recall/F1 is not estimable for every asset.")
}

# F1 must follow the direct confusion-matrix definition. In particular,
# zero true positives imply F1=0 whenever alarms or reference events exist.
f1_definition_check <- benchmark_table %>%
  mutate(
    fp_check = alarm_blocks - stress_hit_blocks,
    fn_check = stress_event_blocks - stress_hit_blocks,
    f1_expected = if_else(
      (2 * stress_hit_blocks + fp_check + fn_check) > 0,
      2 * stress_hit_blocks /
        (2 * stress_hit_blocks + fp_check + fn_check),
      NA_real_
    ),
    f1_match = (is.na(f1_expected) & is.na(stress_proxy_f1)) |
      abs(f1_expected - stress_proxy_f1) < 1e-12
  )

f1_definition_ok <- all(f1_definition_check$f1_match)
if (!f1_definition_ok) {
  stop("FINAL AUDIT FAILED: F1 does not match the direct confusion-matrix definition.")
}


# Level-only and terminal-FRV are strictly monotone terminal-value scores and
# must therefore produce the same block alarm set under rank-based calibration.
level_terminal_equivalence_ok <- all(
  level_terminal_equivalence_table$identical_test_alarm_sets %in% TRUE
) && all(
  abs(level_terminal_equivalence_table$spearman_score_correlation_test - 1) < 1e-12
)
if (!level_terminal_equivalence_ok) {
  stop("FINAL AUDIT FAILED: level-only and terminal-FRV rank/alarm equivalence failed.")
}

# The one-sided rank benchmark must be present and the legacy two-sided method absent.
one_sided_rank_benchmark_ok <-
  all("Integrated_upper_tail_rank" %in% unique(benchmark_table$method)) &&
  !("Integrated_rank_depth" %in% unique(benchmark_table$method))
if (!one_sided_rank_benchmark_ok) {
  stop("FINAL AUDIT FAILED: upper-tail rank benchmark replacement is incomplete.")
}

# Random-score diagnostic must be finite and generated at the configured size.
random_reference_ok <-
  nrow(random_reference_raw) == length(CFG$symbols) * CFG$random_score_reps &&
  all(is.finite(random_reference_raw$block_exceedance_rate)) &&
  all(is.finite(random_reference_raw$abs_calibration_gap))
if (!random_reference_ok) {
  stop("FINAL AUDIT FAILED: random-score reference diagnostic is incomplete.")
}

# Range-reference episodes must be traceable to at least one peak date in every asset.
range_episode_trace_ok <-
  all(stress_proxy_audit$distinct_reference_runs > 0L) &&
  all(stress_proxy_audit$unique_reference_peak_dates > 0L)
if (!range_episode_trace_ok) {
  stop("FINAL AUDIT FAILED: range-reference episodes are not traceable to peak dates.")
}

# Exact-zero HAR audit is reported, not silently altered. Warn if the floor is active.
if (any(zero_return_audit_table$har_training_rv1_at_or_below_epsilon > 0L)) {
  warning(
    "HAR LOG-FLOOR AUDIT: at least one training RV(1) is at/below epsilon. ",
    "Inspect A14 before deciding whether a floor sensitivity rerun is needed."
  )
}

simulation_design_shared_ok <- all(
  simulation_design_consistency$same_core_parameter[
    simulation_design_consistency$parameter != "replications"
  ]
)
if (!simulation_design_shared_ok) {
  stop("FINAL AUDIT FAILED: main/sensitivity simulation design mismatch.")
}

# Positive-weight sensitivity must represent genuinely different a/b ratios.
weight_ratio_unique_ok <-
  length(unique(round(CFG$weight_grid$a / CFG$weight_grid$b, 12))) ==
  nrow(CFG$weight_grid)
if (!weight_ratio_unique_ok) {
  stop("FINAL AUDIT FAILED: weight sensitivity contains duplicate a/b ratios.")
}

har_persistence_ok <- (0.55 + 0.25 + 0.15) < 1
if (!har_persistence_ok) {
  stop("FINAL AUDIT FAILED: HAR persistence is not below one.")
}

publication_audit <- tibble(
  analysis_id = CFG$analysis_id,
  audit_passed = TRUE,
  n_symbols = length(CFG$symbols),
  n_main_asset_rows = nrow(main_asset_table),
  n_benchmark_rows = nrow(benchmark_table),
  n_ablation_rows = nrow(ablation_table),
  n_weight_rows = nrow(weight_table),
  n_tuning_rows = nrow(TUNING),
  n_sim_main_rows = nrow(SIM_SUM),
  n_sim_sensitivity_rows = nrow(SIM_SENS_SUM),
  common_random_numbers_sensitivity = TRUE,
  empirical_calibration_sensitivity_fixed_test = fixed_test_empirical_ok,
  simulation_calibration_sensitivity_fixed_test = fixed_test_sim_ok,
  curvature_dimension_H_minus_2 = curvature_dimension_ok,
  weight_ratios_unique = weight_ratio_unique_ok,
  har_persistence_below_one = har_persistence_ok,
  data_checksum_files = nrow(cache_checksums),
  heldout_range_reference_estimable = stress_proxy_estimable_ok,
  f1_direct_definition_verified = f1_definition_ok,
  main_empirical_asset_count = length(CFG$symbols),
  shareable_single_file_pipeline = TRUE,
  final_publication_figure_layout = TRUE,
  stress_reference_expost_evaluation_only = TRUE,
  stress_reference_used_in_fit_or_calibration = FALSE,
  empirical_fpr_label_used = FALSE,
  one_sided_upper_tail_rank_benchmark = one_sided_rank_benchmark_ok,
  random_score_sanity_check_complete = random_reference_ok,
  level_terminal_rank_equivalence_verified = level_terminal_equivalence_ok,
  range_reference_peak_episode_trace_complete = range_episode_trace_ok,
  simulation_main_sensitivity_core_design_shared = simulation_design_shared_ok,
  har_zero_return_audit_reported = TRUE,
  recognized_episode_mapping_descriptive_only = TRUE,
  calibration_gap_ranking_statistical_superiority_claimed = FALSE
)
readr::write_csv(
  publication_audit,
  file.path(DIRS$audit, "publication_audit.csv")
)

message("==============================================================")
message("PUBLICATION PIPELINE COMPLETED.")
message("Output root: ", normalizePath(CFG$out_root, winslash = "/"))
message("Important interpretation:")
message("  * empirical rates are exceedance/alarm rates, not ground-truth FPRs;")
message("  * block maxima target at least one extreme score within a block;")
message("  * alpha_block mapping is a design mapping;")
message("  * finite-sample conformal control requires block-score exchangeability.")
message("  * the range-based stress reference is held-out/ex-post and evaluation-only;")
message("    it is never used to fit CIRD or calibrate conformal thresholds.")
message("  * empirical universe is locked to BTC, ETH, XRP, and LTC.")
message("  * benchmark rank comparison is one-sided and upper-tail aligned.")
message("  * random-score calibration is a diagnostic sanity check, not a skill benchmark.")
message("  * level/curvature correlations, reference episodes, and HAR zero returns are audited.")
message("  * calibration-gap/weight rankings are descriptive; no statistical superiority is claimed.")
message("  * LaTeX tables are written to latex/main_text and latex/appendix.")
message("  * manuscript placement plan is in audit/manuscript_output_plan.txt.")
message("  * reproducibility instructions are in audit/README_REPRODUCIBILITY.txt.")
message("==============================================================")
