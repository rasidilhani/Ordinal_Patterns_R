# =============================================================================
# Ordinal patterns under transformations of a seismic amplitude series
# -----------------------------------------------------------------------------
# Message for the audience:
#   * Ordinal patterns depend only on the RANK ORDER of values in each window.
#   * Any strictly increasing transformation (x^2, sqrt(x), log(x) for x > 0)
#     keeps the order -> identical patterns, identical permutation entropy.
#   * Variance (an amplitude-based statistic) changes drastically.
#   * A non-monotone transformation, e.g. (x - mean(x))^2, DOES change patterns.
#
# Data: a few days of real RSAM read from the local rsync cache (daily_cache).
# Nothing is downloaded; only files already in the cache are used.
# Base R only, no packages needed.
# =============================================================================

## ---- CONFIG ------------------------------------------------------------------
STATION    <- "WIZ"
CACHE_DIR  <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/daily_cache"
OUT_PDF    <- "C:/Users/UserA1/Documents/GitHub/Seismic_Amplitude_Timeseries_Analysis/results/ordinal_pattern_demo.pdf"

DEMO_START <- as.Date("2019-01-01")  # first day to use (falls back if not cached)
N_DAYS     <- 30                       # how many daily files to read
N_POINTS   <- 30*6*24                     # keep only the first N_POINTS values

D   <- 4   # embedding dimension (pattern length); 3 keeps the bar plots readable
tau <- 1   # embedding delay

## ---- 1. Data from the cache -------------------------------------------------------
# rsync --files-from keeps the server layout, so cached files sit at
#   CACHE_DIR/<year>/<year>.<DOY>/<STATION>.NZ/<year>.<DOY>.<STATION>.timeseries.csv
cached <- list.files(CACHE_DIR,
                     pattern = paste0("\\.", STATION, "\\.timeseries\\.csv$"),
                     recursive = TRUE, full.names = TRUE)
if (length(cached) == 0) {
  stop("No cached ", STATION, " files found in ", CACHE_DIR,
       ". Run the rsync pipeline for a few days first.")
}

# Date of each file from its name "YYYY.DDD.STATION.timeseries.csv"
yd         <- substr(basename(cached), 1, 8)
file_dates <- as.Date(paste0(substr(yd, 1, 4), "-01-01")) +
  as.integer(substr(yd, 6, 8)) - 1
o          <- order(file_dates)
cached     <- cached[o]
file_dates <- file_dates[o]

want <- DEMO_START + 0:(N_DAYS - 1)
use  <- file_dates %in% want
if (!any(use)) {
  message("Days from ", DEMO_START, " are not in the cache; ",
          "using the first ", N_DAYS, " cached days instead.")
  use <- seq_along(cached) <= N_DAYS
}

read_day <- function(f) {
  d <- read.csv(f, check.names = FALSE)   # 1st column = unnamed timestamp
  data.frame(time = d[[1]], rsam = d[["RSAM"]])
}
raw <- do.call(rbind, lapply(cached[use], read_day))

ok <- is.finite(raw$rsam) & raw$rsam > 0   # log and sqrt need x > 0
if (any(!ok)) message("Dropped ", sum(!ok), " missing or non-positive RSAM values.")
raw <- raw[ok, ]
raw <- head(raw, N_POINTS)
x   <- raw$rsam

# Timestamps: parse as UTC first, then display in NZ local time
tt <- as.POSIXct(raw$time, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
if (anyNA(tt)) stop("Could not parse some timestamps, e.g. '", raw$time[which(is.na(tt))[1]], "'")
attr(tt, "tzone") <- "Pacific/Auckland"

data_label <- sprintf("%s RSAM, %s to %s (n = %d)",
                      STATION, raw$time[1], raw$time[nrow(raw)], length(x))
cat("Using", data_label, "\n")
cat("Rows per daily file:",
    paste(sapply(cached[use], function(f) nrow(read.csv(f))), collapse = ", "), "\n")

stopifnot(length(x) > 50, all(x > 0))

## ---- 2. Transformations ------------------------------------------------------
series <- list(
  "Original x"                = x,
  "Square x^2"                = x^2,
  "Square root sqrt(x)"       = sqrt(x),
  "Log log(x)"                = log(x),
  "Non-monotone (x - mean)^2" = (x - mean(x))^2   # contrast case
)

## ---- 3. Ordinal pattern functions ------------------------------------------------
# Pattern of each window = permutation that sorts it, written as a label "123" etc.
# order(..., method = "radix") is stable: tied values are ranked by time order
# (relevant for quantised RSAM, which can contain ties).
ordinal_patterns <- function(x, D = 3, tau = 1) {
  idx <- seq_len(length(x) - (D - 1) * tau)
  vapply(idx, function(i) {
    w <- x[i + (0:(D - 1)) * tau]
    paste(order(w, method = "radix"), collapse = "")
  }, character(1))
}

all_patterns <- function(D) {
  perms <- function(v) if (length(v) == 1) list(v) else
    do.call(c, lapply(seq_along(v), function(i)
      lapply(perms(v[-i]), function(p) c(v[i], p))))
  vapply(perms(1:D), paste, character(1), collapse = "")
}

pattern_probs <- function(pat, D) {
  tab <- table(factor(pat, levels = all_patterns(D)))
  as.numeric(tab) / sum(tab)
}

# Normalised permutation (Shannon) entropy, in [0, 1]
perm_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p)) / log(factorial(D))
}

## ---- 4. Compute everything ---------------------------------------------------------
pats  <- lapply(series, ordinal_patterns, D = D, tau = tau)
probs <- lapply(pats, pattern_probs, D = D)

summary_tab <- data.frame(
  Series            = names(series),
  Mean              = sapply(series, mean),
  Variance          = sapply(series, var),
  CV                = sapply(series, function(s) sd(s) / mean(s)),
  PermEntropy       = sapply(probs, perm_entropy),
  Same_patterns_pct = sapply(pats, function(p) 100 * mean(p == pats[[1]])),
  row.names = NULL
)
print(format(summary_tab, digits = 4, scientific = 3), row.names = FALSE)

cat("\nFirst 6 ordinal patterns (D = ", D, "):\n", sep = "")
print(sapply(pats, head, 6))

## ---- 5. Figure -------------------------------------------------------------------------
cols <- c("#1b6ca8", "#d1495b", "#edae49", "#00798c", "grey30")

# y-axis labels: RSAM, with the transformation applied
ylabs <- list(expression(RSAM),
              expression(RSAM^2),
              expression(sqrt(RSAM)),
              expression(log(RSAM)),
              expression((RSAM - bar(RSAM))^2))

# Date ticks for the time axis: daily for short spans, weekly beyond 10 days
span_days <- as.numeric(difftime(max(tt), min(tt), units = "days"))
day_ticks <- seq(trunc(min(tt), "days"), max(tt),
                 by = if (span_days > 10) "week" else "day")

# Drawing code lives in a function so it can go to the Plots pane AND the PDF.
draw_figure <- function() {
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) })
  
  nr <- length(series)
  # Last row is a bit taller to make room for the shared x-axis labels
  layout(matrix(1:(2 * nr), ncol = 2, byrow = TRUE),
         widths = c(2.3, 1), heights = c(rep(1, nr - 1), 1.45))
  par(oma = c(0, 0, 0, 0), cex.main = 1.05, family = "serif")
  ymax <- max(unlist(probs)) * 1.15
  pat_names <- all_patterns(D)
  
  for (k in seq_len(nr)) {
    s      <- series[[k]]
    last   <- k == nr                 # only the bottom row gets x-axis labels
    bottom <- if (last) 4.5 else 1
    
    # Time series panel
    par(mar = c(bottom, 4.5, 2, 1), mgp = c(2.2, 0.7, 0))
    plot(tt, s, type = "l", col = cols[k], lwd = 1, xaxt = "n",
         xlab = "", ylab = ylabs[[k]],
         main = sprintf("%s    var = %s", names(series)[k],
                        format(var(s), digits = 3, scientific = TRUE)))
    axis.POSIXct(1, at = day_ticks, format = "%d %b %Y", labels = last)
    if (last) mtext("Time", side = 1, line = 2.5)
    
    # Pattern distribution panel: title shows only H; D goes on the bottom axis
    par(mar = c(bottom, 4.5, 2, 1), mgp = c(3, 0.7, 0))
    barplot(probs[[k]], col = cols[k], border = NA, ylim = c(0, ymax),
            names.arg = if (last) pat_names else rep("", length(pat_names)),
            cex.names = 0.7, las = 2,
            ylab = if (k == ceiling(nr / 2)) "Relative frequency" else "",
            main = sprintf("H = %.3f", perm_entropy(probs[[k]])))
    if (last) mtext(sprintf("Ordinal pattern (D = %d)", D),
                    side = 1, line = 3.3, cex = 0.8)
  }
}

# 1) Show on screen. Try the RStudio Plots pane first; if the pane is too
#    small for 10 panels ("figure margins too large"), open a separate,
#    full-size plot window instead.
shown <- tryCatch({ draw_figure(); TRUE }, error = function(e) FALSE)
if (!shown) {
  message("Plots pane too small; opening a separate plot window.")
  dev.new(width = 11, height = 12, noRStudioGD = TRUE)
  draw_figure()
}

# 2) Save the same figure to PDF.
dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 11, height = 12)
draw_figure()
dev.off()
cat("\nFigure saved to", OUT_PDF, "\n")

#################################################################
#Using synthetic data
# =============================================================================
# Ordinal patterns under transformations of a SYNTHETIC amplitude series
# -----------------------------------------------------------------------------
# Message for the audience:
#   * Ordinal patterns depend only on the RANK ORDER of values in each window.
#   * Any strictly increasing transformation (x^2, sqrt(x), log(x) for x > 0)
#     keeps the order -> identical patterns, identical permutation entropy.
#   * Variance (an amplitude-based statistic) changes drastically.
#   * A non-monotone transformation, e.g. (x - mean(x))^2, DOES change patterns.
#
# Data: simulated with a random number generator (no real-world data).
# Base R only, no packages needed.
# =============================================================================

## ---- CONFIG ------------------------------------------------------------------
SEED       <- 123456                       # change for a different random series
START_TIME <- as.POSIXct("2019-12-05 00:00:00", tz = "Pacific/Auckland")
STEP_MIN   <- 10                         # sampling interval (minutes)
N_DAYS     <- 30                         # length of the series in days
N_POINTS   <- N_DAYS * 24 * 60 / STEP_MIN   # 144 points per day at 10 min

OUT_PDF    <- "ordinal_pattern_demo_synthetic.pdf"

D   <- 4   # embedding dimension (pattern length)
tau <- 1   # embedding delay

## ---- 1. Synthetic data ----------------------------------------------------------
# An amplitude-like series: positive, persistent background (log-AR(1))
# plus a few transient bursts that decay exponentially.
set.seed(SEED)
n <- N_POINTS

z <- as.numeric(arima.sim(list(ar = 0.8), n = n, sd = 0.25))

bursts <- rep(0, n)
burst_at <- round(n * c(0.25, 0.60, 0.82))        # burst positions
bursts[burst_at] <- c(1.5, 2.2, 1.2)              # burst sizes (log scale)
bursts <- as.numeric(stats::filter(bursts, 0.85, method = "recursive"))

x  <- 200 * exp(z + bursts)                       # always > 0
tt <- START_TIME + (0:(n - 1)) * STEP_MIN * 60    # synthetic time stamps

stopifnot(all(x > 0))                             # needed for log and sqrt

cat(sprintf("Synthetic series: %s to %s (n = %d)\n",
            format(min(tt), "%d %b %Y %H:%M"),
            format(max(tt), "%d %b %Y %H:%M"), n))

## ---- 2. Transformations ------------------------------------------------------
series <- list(
  "Original x"                = x,
  "Square x^2"                = x^2,
  "Square root sqrt(x)"       = sqrt(x),
  "Log log(x)"                = log(x),
  "Non-monotone (x - mean)^2" = (x - mean(x))^2   # contrast case
)

## ---- 3. Ordinal pattern functions ------------------------------------------------
# Pattern of each window = permutation that sorts it, written as a label "1234" etc.
# order(..., method = "radix") is stable: ties are ranked by time order.
ordinal_patterns <- function(x, D = 3, tau = 1) {
  idx <- seq_len(length(x) - (D - 1) * tau)
  vapply(idx, function(i) {
    w <- x[i + (0:(D - 1)) * tau]
    paste(order(w, method = "radix"), collapse = "")
  }, character(1))
}

all_patterns <- function(D) {
  perms <- function(v) if (length(v) == 1) list(v) else
    do.call(c, lapply(seq_along(v), function(i)
      lapply(perms(v[-i]), function(p) c(v[i], p))))
  vapply(perms(1:D), paste, character(1), collapse = "")
}

pattern_probs <- function(pat, D) {
  tab <- table(factor(pat, levels = all_patterns(D)))
  as.numeric(tab) / sum(tab)
}

# Normalised permutation (Shannon) entropy, in [0, 1]
perm_entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p)) / log(factorial(D))
}

## ---- 4. Compute everything ---------------------------------------------------------
pats  <- lapply(series, ordinal_patterns, D = D, tau = tau)
probs <- lapply(pats, pattern_probs, D = D)

summary_tab <- data.frame(
  Series            = names(series),
  Mean              = sapply(series, mean),
  Variance          = sapply(series, var),
  CV                = sapply(series, function(s) sd(s) / mean(s)),
  PermEntropy       = sapply(probs, perm_entropy),
  Same_patterns_pct = sapply(pats, function(p) 100 * mean(p == pats[[1]])),
  row.names = NULL
)
print(format(summary_tab, digits = 4, scientific = 3), row.names = FALSE)

cat("\nFirst 6 ordinal patterns (D = ", D, "):\n", sep = "")
print(sapply(pats, head, 6))

## ---- 5. Figure -------------------------------------------------------------------------
cols <- c("#1b6ca8", "#d1495b", "#edae49", "#00798c", "grey30")

# y-axis labels, with the transformation applied
ylabs <- list(expression(x),
              expression(x^2),
              expression(sqrt(x)),
              expression(log(x)),
              expression((x - bar(x))^2))

# Date ticks: daily for short spans, weekly beyond 10 days
span_days <- as.numeric(difftime(max(tt), min(tt), units = "days"))
day_ticks <- seq(trunc(min(tt), "days"), max(tt),
                 by = if (span_days > 10) "week" else "day")
date_fmt  <- if (span_days > 10) "%d %b" else "%d %b %Y"

draw_figure <- function() {
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); layout(1) })
  
  nr <- length(series)
  # Last row is a bit taller to make room for the shared x-axis labels
  layout(matrix(1:(2 * nr), ncol = 2, byrow = TRUE),
         widths = c(2.3, 1), heights = c(rep(1, nr - 1), 1.45))
  par(oma = c(0, 0, 0, 0), cex.main = 1.05, family = "serif")
  ymax <- max(unlist(probs)) * 1.15
  pat_names <- all_patterns(D)
  
  for (k in seq_len(nr)) {
    s      <- series[[k]]
    last   <- k == nr                 # only the bottom row gets x-axis labels
    bottom <- if (last) 4.5 else 1
    
    # Time series panel
    par(mar = c(bottom, 4.5, 2, 1), mgp = c(2.2, 0.7, 0))
    plot(tt, s, type = "l", col = cols[k], lwd = 0.6, xaxt = "n",
         xlab = "", ylab = ylabs[[k]],
         main = sprintf("%s    var = %s", names(series)[k],
                        format(var(s), digits = 3, scientific = TRUE)))
    axis.POSIXct(1, at = day_ticks, format = date_fmt, labels = last)
    if (last) mtext("Time", side = 1, line = 2.5)
    
    # Pattern distribution panel: title shows only H; D goes on the bottom axis
    par(mar = c(bottom, 4.5, 2, 1), mgp = c(3, 0.7, 0))
    barplot(probs[[k]], col = cols[k], border = NA, ylim = c(0, ymax),
            names.arg = if (last) pat_names else rep("", length(pat_names)),
            cex.names = 0.7, las = 2,
            ylab = if (k == ceiling(nr / 2)) "Relative frequency" else "",
            main = sprintf("H = %.3f", perm_entropy(probs[[k]])))
    if (last) mtext(sprintf("Ordinal pattern (D = %d)", D),
                    side = 1, line = 3.3, cex = 0.8)
  }
}

# 1) Show on screen. Try the RStudio Plots pane first; if it is too small
#    ("figure margins too large"), open a separate full-size window.
shown <- tryCatch({ draw_figure(); TRUE }, error = function(e) FALSE)
if (!shown) {
  message("Plots pane too small; opening a separate plot window.")
  dev.new(width = 11, height = 12, noRStudioGD = TRUE)
  draw_figure()
}

# 2) Save the same figure to PDF (in the working directory).
pdf(OUT_PDF, width = 11, height = 12)
draw_figure()
dev.off()
cat("\nFigure saved to", normalizePath(OUT_PDF), "\n")