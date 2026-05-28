# =========================================================
# Beyond Volatility: Leakage-Safe Residual-Stress Signal
# =========================================================

options(stringsAsFactors = FALSE)
options(xts.message.period.apply.mean = FALSE)
graphics.off()
rm(list = ls())

# =========================================================
# 0. Packages
# =========================================================

cat("STEP 0: Loading packages...\n")

pkg_needed <- c(
  "quantmod", "xts", "zoo", "PerformanceAnalytics", "pROC", "PRROC", "scales"
)

pkg_install_if_missing <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss) > 0) install.packages(miss, dependencies = TRUE)
}

pkg_install_if_missing(pkg_needed)

suppressPackageStartupMessages({
  library(quantmod)
  library(xts)
  library(zoo)
  library(PerformanceAnalytics)
  library(pROC)
  library(PRROC)
  library(scales)
})

cat("Packages loaded successfully.\n\n")

# =========================================================
# 1. CONFIG
# =========================================================

CFG <- list(
  start = as.Date("2015-01-01"),
  end   = as.Date("2025-12-31"),
  
  subperiods = list(
    Main_2015_2025 = c(as.Date("2015-01-01"), as.Date("2025-12-31")),
    PrePandemic_2015_2019 = c(as.Date("2015-01-01"), as.Date("2019-12-31")),
    StressHeavy_2020_2025 = c(as.Date("2020-01-01"), as.Date("2025-12-31"))
  ),
  
  ticker_spy = "SPY",
  tickers_sectors = c("XLC", "XLY", "XLP", "XLE", "XLF", "XLV", "XLI", "XLK", "XLB", "XLU", "XLRE"),
  
  K = 2,
  pca_type = "expanding",
  pca_burn_in = 252,
  pca_rolling_windows = c(252, 504),
  
  threshold_window = 252,
  stress_q = 0.95,
  vol_q = 0.95,
  corr_q = 0.95,
  
  drawdown_threshold = -0.10,
  warning_horizon = 21,
  
  vol_window = 21,
  corr_window = 63,
  
  event_alarm_lookback = 63,
  event_study_window = 30,
  
  overlay_lambda = 0.30,
  baseline_cost_bps = 5,
  cost_scenarios_bps = c(5, 10, 25),
  
  show_table_views = FALSE,
  show_diagnostic_plots = TRUE
)

UNIV <- CFG$tickers_sectors
TICKERS_ALL <- c(CFG$ticker_spy, UNIV)

COL <- list(
  stress = "#1B9E77",
  threshold = "#B2182B",
  drawdown = "#333333",
  onset = "#D95F02",
  vol = "#7570B3",
  corr = "#E6AB02",
  spy = "#000000",
  overlay_stress = "#1F78B4",
  overlay_vol = "#A6761D",
  grid = "#D0D7E2"
)

# =========================================================
# 2. Utility functions
# =========================================================

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

round_df_numeric <- function(df, digits = 6) {
  out <- df
  is_num <- sapply(out, is.numeric)
  out[is_num] <- lapply(out[is_num], round, digits = digits)
  out
}

safe_last <- function(x, default = NA_real_) {
  if (is.null(x) || NROW(x) == 0) return(default)
  v <- safe_numeric(last(x))
  if (length(v) == 0 || !is.finite(v)) return(default)
  v
}

finite_range <- function(x, pad_frac = 0.08, default = c(0, 1)) {
  z <- safe_numeric(x)
  z <- z[is.finite(z)]
  if (length(z) == 0) return(default)
  r <- range(z, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) <= 0) {
    pad <- max(abs(r[1]) * 0.05, 1e-4)
    return(c(r[1] - pad, r[2] + pad))
  }
  pad <- diff(r) * pad_frac
  c(r[1] - pad, r[2] + pad)
}

print_table <- function(df, table_name) {
  cat("\n====================================================\n")
  cat(table_name, "\n")
  cat("====================================================\n")
  print(round_df_numeric(df, 8), row.names = FALSE)
}

# Alias: prints only, does not save a file.
save_table <- function(df, filename) {
  table_name <- gsub("_", " ", filename)
  table_name <- gsub("\\.csv$", "", table_name)
  print_table(df, table_name)
}

base_par <- function(mar = c(4.2, 4.6, 3.0, 0.8) + 0.1,
                     cex_axis = 0.90,
                     cex_lab = 1.00,
                     cex_main = 1.00) {
  par(
    mfrow = c(1, 1),
    mar = mar,
    oma = c(0, 0, 0, 0),
    las = 1,
    bty = "l",
    tck = -0.015,
    mgp = c(2.4, 0.75, 0),
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main
  )
}

add_grid <- function() grid(col = COL$grid, lty = 3)

nice_date_axis <- function(dates) {
  axis.Date(
    1,
    at = seq(min(dates), max(dates), by = "2 years"),
    format = "%Y"
  )
}

plot_message <- function(title, message) {
  base_par(mar = c(3.8, 3.8, 3.0, 0.8) + 0.1)
  plot.new()
  title(main = title)
  wrapped <- paste(strwrap(message, width = 85), collapse = "\n")
  text(0.5, 0.5, wrapped, cex = 0.90)
}

# Robust NAV: first NA return is treated as zero so a curve starts at 1.
calc_nav <- function(ret_xts) {
  r <- safe_numeric(ret_xts)
  r[!is.finite(r)] <- 0
  nav <- cumprod(1 + r)
  if (length(nav) > 0 && is.finite(nav[1]) && abs(nav[1]) > 1e-12) {
    nav <- nav / nav[1]
  }
  out <- xts(nav, order.by = index(ret_xts))
  colnames(out) <- paste0(colnames(ret_xts)[1], "_NAV")
  out
}

calc_drawdown <- function(ret_xts) {
  nav <- calc_nav(ret_xts)
  x <- safe_numeric(nav)
  dd <- x / cummax(x) - 1
  out <- xts(dd, order.by = index(nav))
  colnames(out) <- "Drawdown"
  out
}

ann_stats <- function(ret_xts) {
  if (is.null(ret_xts) || NROW(ret_xts) < 2) {
    return(c(Ann_Return = NA_real_, Ann_Vol = NA_real_, Sharpe = NA_real_, MaxDD = NA_real_))
  }
  r <- safe_numeric(ret_xts)
  r <- r[is.finite(r)]
  if (length(r) < 2) {
    return(c(Ann_Return = NA_real_, Ann_Vol = NA_real_, Sharpe = NA_real_, MaxDD = NA_real_))
  }
  ann_ret <- prod(1 + r, na.rm = TRUE)^(252 / length(r)) - 1
  ann_vol <- sd(r, na.rm = TRUE) * sqrt(252)
  sharpe <- ifelse(is.finite(ann_vol) && ann_vol > 1e-12, ann_ret / ann_vol, NA_real_)
  nav <- cumprod(1 + r)
  maxdd <- min(nav / cummax(nav) - 1, na.rm = TRUE)
  c(Ann_Return = ann_ret, Ann_Vol = ann_vol, Sharpe = sharpe, MaxDD = maxdd)
}

render_figure <- function(filename_stub, plot_fun) {
  cat("\nDisplaying figure:", gsub("_", " ", filename_stub), "\n")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(try(par(oldpar), silent = TRUE), add = TRUE)
  
  tryCatch(
    plot_fun(),
    error = function(e) {
      plot_message(
        paste("Plot failed:", gsub("_", " ", filename_stub)),
        paste("Error:", e$message)
      )
    }
  )
  
  invisible(NULL)
}

# =========================================================
# 3. Data download and construction
# =========================================================

cat("STEP 1: Downloading adjusted prices from Yahoo Finance...\n")

for (tk in TICKERS_ALL) {
  if (exists(tk, inherits = TRUE)) rm(list = tk, inherits = TRUE)
}

invisible(getSymbols(
  Symbols = TICKERS_ALL,
  src = "yahoo",
  from = CFG$start,
  to = CFG$end + 1,
  auto.assign = TRUE,
  warnings = FALSE
))

adj_prices_raw <- do.call(merge, lapply(TICKERS_ALL, function(tk) Ad(get(tk))))
colnames(adj_prices_raw) <- TICKERS_ALL

adj_prices <- adj_prices_raw[paste(CFG$start, CFG$end, sep = "/")]
adj_prices <- adj_prices[complete.cases(adj_prices), ]

ret_simple <- na.omit(ROC(adj_prices, type = "discrete"))
colnames(ret_simple) <- colnames(adj_prices)

spy_ret <- ret_simple[, CFG$ticker_spy, drop = FALSE]
sector_ret <- ret_simple[, UNIV]
sector_excess <- sector_ret - as.numeric(spy_ret)
colnames(sector_excess) <- UNIV

cat("Data range:", as.character(first(index(ret_simple))), "to", as.character(last(index(ret_simple))), "\n")
cat("Price dimensions:", paste(dim(adj_prices), collapse = " x "), "\n\n")

# =========================================================
# 4. Leakage-safe PCA residual-stress construction
# =========================================================

compute_oos_residual_stress <- function(X_xts, K = 2, burn_in = 252,
                                        pca_type = "expanding", rolling_window = NULL,
                                        scale_residual = FALSE, idio_window = 252) {
  X <- as.matrix(X_xts)
  dates <- index(X_xts)
  n <- nrow(X)
  p <- ncol(X)
  
  stress <- rep(NA_real_, n)
  scaled_stress <- rep(NA_real_, n)
  residual_mat <- matrix(NA_real_, nrow = n, ncol = p)
  colnames(residual_mat) <- colnames(X_xts)
  
  for (i in seq_len(n)) {
    if (i <= burn_in) next
    
    if (pca_type == "rolling") {
      if (is.null(rolling_window)) stop("rolling_window must be supplied for rolling PCA.")
      start_i <- max(1, i - rolling_window)
      train_idx <- start_i:(i - 1)
    } else {
      train_idx <- 1:(i - 1)
    }
    
    X_train <- X[train_idx, , drop = FALSE]
    X_train <- X_train[complete.cases(X_train), , drop = FALSE]
    x_t <- X[i, ]
    
    if (nrow(X_train) < max(30, K + 5) || any(!is.finite(x_t))) next
    
    mu <- colMeans(X_train, na.rm = TRUE)
    X_centered <- sweep(X_train, 2, mu, "-")
    
    pca_fit <- tryCatch(prcomp(X_centered, center = FALSE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pca_fit)) next
    
    K_use <- min(K, ncol(pca_fit$rotation))
    L <- pca_fit$rotation[, 1:K_use, drop = FALSE]
    
    xt_centered <- x_t - mu
    fitted <- L %*% (t(L) %*% xt_centered)
    e_t <- as.numeric(xt_centered - fitted)
    
    residual_mat[i, ] <- e_t
    stress[i] <- sqrt(mean(e_t^2, na.rm = TRUE))
    
    if (isTRUE(scale_residual)) {
      prev_resid <- residual_mat[max(1, i - idio_window):(i - 1), , drop = FALSE]
      idio_sd <- apply(prev_resid, 2, sd, na.rm = TRUE)
      idio_sd[!is.finite(idio_sd) | idio_sd <= 1e-12] <- NA_real_
      z <- e_t / idio_sd
      scaled_stress[i] <- sqrt(mean(z^2, na.rm = TRUE))
    }
  }
  
  stress_xts <- xts(stress, order.by = dates)
  colnames(stress_xts) <- "ResidualStress"
  
  scaled_xts <- xts(scaled_stress, order.by = dates)
  colnames(scaled_xts) <- "ScaledResidualStress"
  
  residual_xts <- xts(residual_mat, order.by = dates)
  colnames(residual_xts) <- colnames(X_xts)
  
  list(stress = stress_xts, scaled_stress = scaled_xts, residuals = residual_xts)
}

rolling_train_quantile <- function(x_xts, window, q) {
  q_raw <- rollapplyr(
    x_xts,
    width = window,
    FUN = function(z) quantile(as.numeric(z), probs = q, na.rm = TRUE),
    fill = NA,
    by.column = TRUE
  )
  out <- lag(q_raw, 1)
  colnames(out) <- paste0(colnames(x_xts), "_Threshold")
  out
}

compute_realized_vol <- function(spy_ret_xts, window = 21) {
  vol <- rollapplyr(spy_ret_xts, width = window, FUN = sd, fill = NA) * sqrt(252)
  colnames(vol) <- "RealizedVol"
  vol
}

compute_pairwise_corr <- function(excess_xts, window = 63) {
  X <- as.matrix(excess_xts)
  out <- rep(NA_real_, nrow(X))
  for (i in seq_len(nrow(X))) {
    if (i < window) next
    Z <- X[(i - window + 1):i, , drop = FALSE]
    if (any(!is.finite(Z))) next
    C <- cor(Z)
    out[i] <- mean(C[upper.tri(C)], na.rm = TRUE)
  }
  ans <- xts(out, order.by = index(excess_xts))
  colnames(ans) <- "PairwiseCorr"
  ans
}

cat("STEP 2: Constructing leakage-safe residual stress...\n")

base_sig <- compute_oos_residual_stress(
  sector_excess,
  K = CFG$K,
  burn_in = CFG$pca_burn_in,
  pca_type = CFG$pca_type
)

stress <- base_sig$stress
residuals <- base_sig$residuals
stress_thr <- rolling_train_quantile(stress, CFG$threshold_window, CFG$stress_q)

vol <- compute_realized_vol(spy_ret, CFG$vol_window)
vol_thr <- rolling_train_quantile(vol, CFG$threshold_window, CFG$vol_q)

corr_score <- compute_pairwise_corr(sector_excess, CFG$corr_window)
corr_thr <- rolling_train_quantile(corr_score, CFG$threshold_window, CFG$corr_q)

# =========================================================
# 5. Drawdown-onset and H-day early-warning label
# =========================================================

spy_nav <- calc_nav(spy_ret)
colnames(spy_nav) <- "SPY_NAV"

spy_dd <- spy_nav / cummax(spy_nav) - 1
colnames(spy_dd) <- "SPY_Drawdown"

dd_region <- spy_dd <= CFG$drawdown_threshold
dd_region_lag <- lag(dd_region, 1)
dd_region_lag[is.na(dd_region_lag)] <- FALSE

onset <- dd_region & !dd_region_lag
colnames(onset) <- "DrawdownOnset"

make_horizon_label <- function(onset_xts, H) {
  y <- rep(0, NROW(onset_xts))
  onset_vec <- as.logical(onset_xts)
  onset_vec[is.na(onset_vec)] <- FALSE
  n <- length(onset_vec)
  
  for (i in seq_len(n)) {
    j1 <- i + 1
    j2 <- min(n, i + H)
    if (j1 <= j2) y[i] <- as.integer(any(onset_vec[j1:j2]))
  }
  
  out <- xts(y, order.by = index(onset_xts))
  colnames(out) <- paste0("Y_H", H)
  out
}

y_H <- make_horizon_label(onset, CFG$warning_horizon)

# =========================================================
# 6. Metrics and diagnostics
# =========================================================

compute_auc_metrics <- function(y_xts, score_xts, threshold_xts = NULL, name = "Signal") {
  z <- na.omit(merge(y_xts, score_xts, all = FALSE))
  colnames(z) <- c("y", "score")
  
  y <- safe_numeric(z$y)
  s <- safe_numeric(z$score)
  keep <- is.finite(y) & is.finite(s)
  y <- y[keep]
  s <- s[keep]
  
  if (length(y) < 2 || length(unique(y)) < 2) {
    roc_auc <- NA_real_
    pr_auc <- NA_real_
  } else {
    roc_auc <- safe_numeric(pROC::auc(pROC::roc(y, s, quiet = TRUE, direction = "<")))
    fg <- s[y == 1]
    bg <- s[y == 0]
    pr_auc <- tryCatch(
      PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)$auc.integral,
      error = function(e) NA_real_
    )
  }
  
  if (!is.null(threshold_xts)) {
    zz <- na.omit(merge(y_xts, score_xts, threshold_xts, all = FALSE))
    yy <- safe_numeric(zz[, 1])
    ss <- safe_numeric(zz[, 2])
    tt <- safe_numeric(zz[, 3])
    
    pred <- as.integer(ss > tt)
    TP <- sum(pred == 1 & yy == 1, na.rm = TRUE)
    FP <- sum(pred == 1 & yy == 0, na.rm = TRUE)
    FN <- sum(pred == 0 & yy == 1, na.rm = TRUE)
    TN <- sum(pred == 0 & yy == 0, na.rm = TRUE)
    
    precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
    recall <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
    f1 <- ifelse(is.finite(precision + recall) && (precision + recall) > 0,
                 2 * precision * recall / (precision + recall), NA_real_)
  } else {
    TP <- FP <- FN <- TN <- precision <- recall <- f1 <- NA_real_
  }
  
  data.frame(
    Signal = name,
    N = length(y),
    Positive_Rate = mean(y == 1, na.rm = TRUE),
    ROC_AUC = roc_auc,
    PR_AUC = pr_auc,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    TP = TP, FP = FP, FN = FN, TN = TN,
    stringsAsFactors = FALSE
  )
}

stress_metrics <- compute_auc_metrics(y_H, stress, stress_thr, "Residual stress")
vol_metrics <- compute_auc_metrics(y_H, vol, vol_thr, "SPY realized volatility")
corr_metrics <- compute_auc_metrics(y_H, corr_score, corr_thr, "Pairwise sector correlation")

Table_1 <- rbind(stress_metrics, vol_metrics, corr_metrics)
save_table(Table_1, "Table_1_early_warning_metrics.csv")
save_table(Table_1[Table_1$Signal %in% c("Residual stress", "SPY realized volatility"), ],
           "Table_2A_standalone_metrics.csv")

conditional_regime_table <- function(y_xts, stress_xts, stress_threshold_xts, vol_xts, vol_threshold_xts) {
  z <- na.omit(merge(y_xts, stress_xts, stress_threshold_xts, vol_xts, vol_threshold_xts, all = FALSE))
  colnames(z) <- c("y", "stress", "stress_thr", "vol", "vol_thr")
  
  y <- safe_numeric(z$y)
  high_stress <- safe_numeric(z$stress) > safe_numeric(z$stress_thr)
  high_vol <- safe_numeric(z$vol) > safe_numeric(z$vol_thr)
  
  regimes <- data.frame(
    Stress_Regime = ifelse(high_stress, "High stress", "Low stress"),
    Vol_Regime = ifelse(high_vol, "High volatility", "Low volatility"),
    y = y
  )
  
  cells <- expand.grid(
    Stress_Regime = c("Low stress", "High stress"),
    Vol_Regime = c("Low volatility", "High volatility"),
    stringsAsFactors = FALSE
  )
  
  out <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    ss <- cells$Stress_Regime[i]
    vv <- cells$Vol_Regime[i]
    sub <- regimes[regimes$Stress_Regime == ss & regimes$Vol_Regime == vv, ]
    
    data.frame(
      Stress_Regime = ss,
      Vol_Regime = vv,
      N = nrow(sub),
      Onset_Probability_H = ifelse(nrow(sub) > 0, mean(sub$y == 1, na.rm = TRUE), NA_real_),
      stringsAsFactors = FALSE
    )
  }))
  
  q00 <- out$Onset_Probability_H[out$Stress_Regime == "Low stress" & out$Vol_Regime == "Low volatility"]
  q10 <- out$Onset_Probability_H[out$Stress_Regime == "High stress" & out$Vol_Regime == "Low volatility"]
  
  out$Q10_minus_Q00 <- NA_real_
  out$Q10_minus_Q00[1] <- q10 - q00
  out
}

Table_2B <- conditional_regime_table(y_H, stress, stress_thr, vol, vol_thr)
save_table(Table_2B, "Table_2B_conditional_regimes.csv")

event_alarm_diagnostics <- function(onset_xts, stress_xts, stress_thr_xts,
                                    vol_xts, vol_thr_xts, lookback = 63) {
  z <- na.omit(merge(onset_xts, stress_xts, stress_thr_xts, vol_xts, vol_thr_xts, all = FALSE))
  colnames(z) <- c("onset", "stress", "stress_thr", "vol", "vol_thr")
  
  dates <- index(z)
  onset_vec <- as.logical(z$onset)
  stress_alarm <- safe_numeric(z$stress) > safe_numeric(z$stress_thr)
  vol_alarm <- safe_numeric(z$vol) > safe_numeric(z$vol_thr)
  
  onset_idx <- which(onset_vec)
  
  if (length(onset_idx) == 0) {
    overlap <- data.frame(Category = c("Both", "Stress only", "Vol only", "Neither"),
                          Count = 0, Share = NA_real_)
    lead <- data.frame(Metric = c("Median stress lead", "Median volatility lead",
                                  "Mean stress lead", "Mean volatility lead",
                                  "Stress earlier share"),
                       Value = NA_real_)
    return(list(overlap = overlap, lead = lead, event_detail = data.frame()))
  }
  
  detail <- do.call(rbind, lapply(onset_idx, function(i) {
    j1 <- max(1, i - lookback)
    j2 <- i - 1
    
    if (j2 < j1) {
      sa_idx <- integer(0)
      va_idx <- integer(0)
    } else {
      sa_idx <- which(stress_alarm[j1:j2]) + j1 - 1
      va_idx <- which(vol_alarm[j1:j2]) + j1 - 1
    }
    
    stress_flag <- length(sa_idx) > 0
    vol_flag <- length(va_idx) > 0
    
    stress_lead <- ifelse(stress_flag, i - max(sa_idx), NA_real_)
    vol_lead <- ifelse(vol_flag, i - max(va_idx), NA_real_)
    
    data.frame(
      Onset_Date = as.character(dates[i]),
      Stress_Flag = stress_flag,
      Vol_Flag = vol_flag,
      Stress_Lead_Days = stress_lead,
      Vol_Lead_Days = vol_lead,
      stringsAsFactors = FALSE
    )
  }))
  
  both <- sum(detail$Stress_Flag & detail$Vol_Flag)
  stress_only <- sum(detail$Stress_Flag & !detail$Vol_Flag)
  vol_only <- sum(!detail$Stress_Flag & detail$Vol_Flag)
  neither <- sum(!detail$Stress_Flag & !detail$Vol_Flag)
  total <- nrow(detail)
  
  overlap <- data.frame(
    Category = c("Both", "Stress only", "Vol only", "Neither"),
    Count = c(both, stress_only, vol_only, neither),
    Share = c(both, stress_only, vol_only, neither) / total,
    stringsAsFactors = FALSE
  )
  
  paired <- detail[detail$Stress_Flag & detail$Vol_Flag, ]
  
  lead <- data.frame(
    Metric = c("Median stress lead", "Median volatility lead",
               "Mean stress lead", "Mean volatility lead",
               "Stress earlier share"),
    Value = c(
      median(paired$Stress_Lead_Days, na.rm = TRUE),
      median(paired$Vol_Lead_Days, na.rm = TRUE),
      mean(paired$Stress_Lead_Days, na.rm = TRUE),
      mean(paired$Vol_Lead_Days, na.rm = TRUE),
      mean(paired$Stress_Lead_Days > paired$Vol_Lead_Days, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  list(overlap = overlap, lead = lead, event_detail = detail)
}

event_diag <- event_alarm_diagnostics(onset, stress, stress_thr, vol, vol_thr, CFG$event_alarm_lookback)
save_table(event_diag$overlap, "Table_2C_event_overlap.csv")
save_table(event_diag$lead, "Table_2D_lead_time.csv")
save_table(event_diag$event_detail, "Event_Detail_overlap_lead_time.csv")

# =========================================================
# 7. Figure functions
# =========================================================

plot_figure_1A <- function() {
  z <- na.omit(merge(stress, stress_thr, all = FALSE))
  
  if (NROW(z) < 2) {
    plot_message("Figure 1A. Residual stress", "No aligned stress and threshold values.")
    return(invisible(NULL))
  }
  
  dates <- index(z)
  
  base_par(mar = c(4.0, 4.7, 3.0, 0.7) + 0.1)
  plot(
    dates,
    safe_numeric(z$ResidualStress),
    type = "l",
    lwd = 1.3,
    col = COL$stress,
    xaxt = "n",
    xlab = "Date",
    ylab = "Residual stress",
    main = "Figure 1A. Residual stress with train-only threshold",
    ylim = finite_range(c(safe_numeric(z$ResidualStress), safe_numeric(z[, 2])), 0.10)
  )
  nice_date_axis(dates)
  lines(dates, safe_numeric(z[, 2]), lwd = 1.4, lty = 2, col = COL$threshold)
  add_grid()
  
  legend(
    "topright",
    legend = c("Residual stress", "Train-only threshold"),
    col = c(COL$stress, COL$threshold),
    lty = c(1, 2),
    lwd = c(1.3, 1.4),
    bty = "o",
    bg = adjustcolor("white", 0.92),
    cex = 0.72
  )
}

plot_figure_1B <- function() {
  z <- na.omit(merge(spy_dd, onset, all = FALSE))
  
  if (NROW(z) < 2) {
    plot_message("Figure 1B. SPY drawdown", "No aligned drawdown and onset values.")
    return(invisible(NULL))
  }
  
  dates <- index(z)
  
  base_par(mar = c(4.0, 4.7, 3.0, 0.7) + 0.1)
  plot(
    dates,
    safe_numeric(z$SPY_Drawdown),
    type = "l",
    lwd = 1.35,
    col = COL$drawdown,
    xaxt = "n",
    xlab = "Date",
    ylab = "SPY drawdown",
    main = "Figure 1B. SPY drawdown and drawdown-onset events",
    ylim = finite_range(c(safe_numeric(z$SPY_Drawdown), CFG$drawdown_threshold), 0.10)
  )
  nice_date_axis(dates)
  abline(h = CFG$drawdown_threshold, lty = 2, col = COL$threshold, lwd = 1.25)
  
  onset_dates <- dates[as.logical(z$DrawdownOnset)]
  if (length(onset_dates) > 0) rug(onset_dates, col = COL$onset, lwd = 1.8, ticksize = 0.06)
  
  add_grid()
  
  legend(
    "bottomleft",
    legend = c("SPY drawdown", "Drawdown threshold", "Onset rug"),
    col = c(COL$drawdown, COL$threshold, COL$onset),
    lty = c(1, 2, NA),
    pch = c(NA, NA, "|"),
    lwd = c(1.35, 1.25, 1.8),
    bty = "o",
    bg = adjustcolor("white", 0.92),
    cex = 0.70
  )
}

plot_figure_2 <- function() {
  z <- na.omit(merge(stress, onset, all = FALSE))
  
  if (NROW(z) < 10 || sum(as.logical(z$DrawdownOnset), na.rm = TRUE) == 0) {
    plot_message("Figure 2. Event-study of residual stress", "No drawdown-onset events are available after alignment.")
    return(invisible(NULL))
  }
  
  s <- safe_numeric(z$ResidualStress)
  onset_idx <- which(as.logical(z$DrawdownOnset))
  W <- CFG$event_study_window
  grid_tau <- -W:W
  mat <- matrix(NA_real_, nrow = length(onset_idx), ncol = length(grid_tau))
  
  for (k in seq_along(onset_idx)) {
    i <- onset_idx[k]
    rng <- (i - W):(i + W)
    ok <- rng >= 1 & rng <= length(s)
    mat[k, ok] <- s[rng[ok]]
  }
  
  keep_cols <- colSums(is.finite(mat)) > 0
  if (!any(keep_cols)) {
    plot_message("Figure 2. Event-study of residual stress", "The event-study matrix has no finite values.")
    return(invisible(NULL))
  }
  
  grid_tau <- grid_tau[keep_cols]
  mat <- mat[, keep_cols, drop = FALSE]
  
  avg <- apply(mat, 2, mean, na.rm = TRUE)
  lo <- apply(mat, 2, function(x) quantile(x, 0.25, na.rm = TRUE))
  hi <- apply(mat, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  plot(
    grid_tau,
    avg,
    type = "n",
    xlab = "Event time (trading days)",
    ylab = "Residual stress",
    main = "Figure 2. Event-study of residual stress around drawdown onsets",
    ylim = finite_range(c(lo, hi, avg), 0.12)
  )
  
  polygon(c(grid_tau, rev(grid_tau)), c(lo, rev(hi)),
          border = NA, col = adjustcolor(COL$stress, 0.22))
  lines(grid_tau, avg, lwd = 2.1, col = COL$stress)
  points(grid_tau, avg, pch = 19, cex = 0.38, col = COL$stress)
  abline(v = 0, lty = 2, col = COL$onset, lwd = 1.35)
  add_grid()
  
  legend(
    "topright",
    legend = c("Cross-event mean", "25%-75% band", "Onset"),
    col = c(COL$stress, adjustcolor(COL$stress, 0.22), COL$onset),
    lty = c(1, NA, 2),
    lwd = c(2.1, NA, 1.35),
    pch = c(19, 15, NA),
    bty = "o",
    bg = adjustcolor("white", 0.92),
    cex = 0.70
  )
}

plot_figure_3A <- function() {
  z1 <- na.omit(merge(y_H, stress, vol, corr_score, all = FALSE))
  
  if (NROW(z1) < 10) {
    plot_message("Figure 3A. ROC curves", "Not enough aligned observations.")
    return(invisible(NULL))
  }
  
  y <- safe_numeric(z1[, 1])
  
  if (length(unique(y[is.finite(y)])) < 2) {
    plot_message("Figure 3A. ROC curves", "ROC is not defined because the label has only one class.")
    return(invisible(NULL))
  }
  
  scores <- list(
    "Residual stress" = safe_numeric(z1[, 2]),
    "Realized volatility" = safe_numeric(z1[, 3]),
    "Pairwise correlation" = safe_numeric(z1[, 4])
  )
  cols <- c(COL$stress, COL$vol, COL$corr)
  
  base_par(mar = c(4.2, 4.6, 3.0, 0.7) + 0.1)
  
  plot(
    0, 0,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "False positive rate",
    ylab = "True positive rate",
    main = "Figure 3A. ROC curves"
  )
  abline(0, 1, lty = 2, col = "gray50")
  add_grid()
  
  leg <- c()
  leg_cols <- c()
  
  for (i in seq_along(scores)) {
    sc <- scores[[i]]
    ok <- is.finite(y) & is.finite(sc)
    
    if (sum(ok) < 10 || length(unique(y[ok])) < 2) next
    
    roc_obj <- pROC::roc(y[ok], sc[ok], quiet = TRUE, direction = "<")
    
    lines(
      1 - roc_obj$specificities,
      roc_obj$sensitivities,
      col = cols[i],
      lwd = 2.1
    )
    
    leg <- c(leg, paste0(names(scores)[i], " (AUC=", sprintf("%.3f", pROC::auc(roc_obj)), ")"))
    leg_cols <- c(leg_cols, cols[i])
  }
  
  if (length(leg) > 0) {
    legend(
      "bottomright",
      legend = leg,
      col = leg_cols,
      lty = 1,
      lwd = 2.1,
      bty = "o",
      bg = adjustcolor("white", 0.92),
      cex = 0.68
    )
  }
}

plot_figure_3B <- function() {
  z1 <- na.omit(merge(y_H, stress, vol, corr_score, all = FALSE))
  
  if (NROW(z1) < 10) {
    plot_message("Figure 3B. Precision-recall curves", "Not enough aligned observations.")
    return(invisible(NULL))
  }
  
  y <- safe_numeric(z1[, 1])
  
  if (sum(y == 1, na.rm = TRUE) == 0 || sum(y == 0, na.rm = TRUE) == 0) {
    plot_message("Figure 3B. Precision-recall curves", "PR curve is not defined because one class is missing.")
    return(invisible(NULL))
  }
  
  scores <- list(
    "Residual stress" = safe_numeric(z1[, 2]),
    "Realized volatility" = safe_numeric(z1[, 3]),
    "Pairwise correlation" = safe_numeric(z1[, 4])
  )
  cols <- c(COL$stress, COL$vol, COL$corr)
  
  base_par(mar = c(4.2, 4.6, 3.0, 0.7) + 0.1)
  
  plot(
    0, 0,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "Recall",
    ylab = "Precision",
    main = "Figure 3B. Precision-recall curves"
  )
  
  base_rate <- mean(y == 1, na.rm = TRUE)
  abline(h = base_rate, lty = 2, col = "gray50")
  add_grid()
  
  leg <- c()
  leg_cols <- c()
  
  for (i in seq_along(scores)) {
    sc <- scores[[i]]
    ok <- is.finite(y) & is.finite(sc)
    
    if (sum(ok) < 10 || sum(y[ok] == 1) == 0 || sum(y[ok] == 0) == 0) next
    
    fg <- sc[ok][y[ok] == 1]
    bg <- sc[ok][y[ok] == 0]
    
    pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
    
    lines(pr$curve[, 1], pr$curve[, 2], col = cols[i], lwd = 2.1)
    
    leg <- c(leg, paste0(names(scores)[i], " (AUC=", sprintf("%.3f", pr$auc.integral), ")"))
    leg_cols <- c(leg_cols, cols[i])
  }
  
  if (length(leg) > 0) {
    legend(
      "topright",
      legend = c(leg, paste0("Baseline=", sprintf("%.3f", base_rate))),
      col = c(leg_cols, "gray50"),
      lty = c(rep(1, length(leg)), 2),
      lwd = c(rep(2.1, length(leg)), 1.3),
      bty = "o",
      bg = adjustcolor("white", 0.92),
      cex = 0.68
    )
  }
}

plot_event_overlap <- function() {
  df <- event_diag$overlap
  
  if (is.null(df) || nrow(df) == 0) {
    plot_message("Diagnostic plot. Event overlap", "No event-overlap table is available.")
    return(invisible(NULL))
  }
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  
  vals <- df$Count
  ylim_max <- max(vals, na.rm = TRUE) * 1.25 + 1
  
  bp <- barplot(
    vals,
    names.arg = df$Category,
    col = c(COL$stress, COL$threshold, COL$vol, "gray70"),
    border = "gray25",
    ylab = "Number of onset events",
    main = "Diagnostic plot. Event overlap by alarm type",
    ylim = c(0, ylim_max)
  )
  
  add_grid()
  
  text(
    bp,
    vals,
    labels = paste0(vals, " (", sprintf("%.1f%%", 100 * df$Share), ")"),
    pos = 3,
    cex = 0.75
  )
}

# =========================================================
# 8. Overlay figures and tables
# =========================================================

stress_regime <- stress > stress_thr
vol_regime <- vol > vol_thr
stress_regime[is.na(stress_regime)] <- FALSE
vol_regime[is.na(vol_regime)] <- FALSE

run_spy_overlay <- function(spy_ret_xts, trigger_xts, lambda = 0.30, cost_bps = 5) {
  z <- merge(spy_ret_xts, trigger_xts, all = FALSE)
  r <- safe_numeric(z[, 1])
  trig <- as.logical(z[, 2])
  trig[is.na(trig)] <- FALSE
  
  trig_lag <- c(FALSE, head(trig, -1))
  exposure <- ifelse(trig_lag, 1 - lambda, 1.0)
  exposure_lag <- c(1.0, head(exposure, -1))
  turnover <- abs(exposure - exposure_lag)
  cost <- turnover * (cost_bps / 10000)
  net_ret <- exposure * r - cost
  
  out <- xts(cbind(net_ret, exposure, turnover, cost), order.by = index(z))
  colnames(out) <- c("Net_Return", "Exposure", "Turnover", "Cost")
  out
}

spy_bh_ret <- spy_ret
colnames(spy_bh_ret) <- "SPY_BuyHold"

stress_overlay <- run_spy_overlay(spy_ret, stress_regime, CFG$overlay_lambda, CFG$baseline_cost_bps)
vol_overlay <- run_spy_overlay(spy_ret, vol_regime, CFG$overlay_lambda, CFG$baseline_cost_bps)

stress_overlay_ret <- stress_overlay[, "Net_Return"]
vol_overlay_ret <- vol_overlay[, "Net_Return"]

colnames(stress_overlay_ret) <- "StressOverlay"
colnames(vol_overlay_ret) <- "VolOverlay"

Table_6 <- data.frame(
  Strategy = c("SPY buy-and-hold", "SPY overlay, residual stress, 30%", "SPY overlay, volatility, 30%"),
  rbind(
    ann_stats(spy_bh_ret),
    ann_stats(stress_overlay_ret),
    ann_stats(vol_overlay_ret)
  ),
  stringsAsFactors = FALSE
)
save_table(Table_6, "Table_6_overlay_performance.csv")

cost_sens <- do.call(rbind, lapply(CFG$cost_scenarios_bps, function(cst) {
  rr <- run_spy_overlay(spy_ret, stress_regime, CFG$overlay_lambda, cst)[, "Net_Return"]
  stats <- ann_stats(rr)
  data.frame(Cost_bps = cst, t(stats), stringsAsFactors = FALSE)
}))
save_table(cost_sens, "Table_7_transaction_cost_sensitivity.csv")

plot_figure_4 <- function() {
  z <- na.omit(merge(calc_nav(spy_bh_ret), calc_nav(stress_overlay_ret), calc_nav(vol_overlay_ret), all = FALSE))
  
  if (NROW(z) < 2) {
    plot_message("Figure 4. Cumulative wealth", "No aligned overlay return data.")
    return(invisible(NULL))
  }
  
  colnames(z) <- c("SPY buy-and-hold", "Stress overlay", "Volatility overlay")
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  
  matplot(
    index(z),
    coredata(z),
    type = "l",
    lty = 1,
    lwd = c(2.0, 2.0, 2.0),
    col = c(COL$spy, COL$overlay_stress, COL$overlay_vol),
    xlab = "Date",
    ylab = "Cumulative wealth",
    main = "Figure 4. Cumulative wealth for SPY and stress-managed overlays"
  )
  
  add_grid()
  
  legend(
    "topleft",
    legend = colnames(z),
    col = c(COL$spy, COL$overlay_stress, COL$overlay_vol),
    lty = 1,
    lwd = 2.0,
    bty = "o",
    bg = adjustcolor("white", 0.90),
    cex = 0.82
  )
}

plot_figure_5 <- function() {
  z <- na.omit(merge(calc_drawdown(spy_bh_ret), calc_drawdown(stress_overlay_ret), calc_drawdown(vol_overlay_ret), all = FALSE))
  
  if (NROW(z) < 2) {
    plot_message("Figure 5. Drawdowns", "No aligned drawdown data.")
    return(invisible(NULL))
  }
  
  colnames(z) <- c("SPY buy-and-hold", "Stress overlay", "Volatility overlay")
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  
  matplot(
    index(z),
    coredata(z),
    type = "l",
    lty = 1,
    lwd = c(2.0, 2.0, 2.0),
    col = c(COL$spy, COL$overlay_stress, COL$overlay_vol),
    xlab = "Date",
    ylab = "Drawdown",
    main = "Figure 5. Drawdowns for SPY and stress-managed overlays"
  )
  
  add_grid()
  
  legend(
    "bottomleft",
    legend = colnames(z),
    col = c(COL$spy, COL$overlay_stress, COL$overlay_vol),
    lty = 1,
    lwd = 2.0,
    bty = "o",
    bg = adjustcolor("white", 0.90),
    cex = 0.82
  )
}

# =========================================================
# 9. Robustness tables
# =========================================================

run_period_metrics <- function(period_name, start_date, end_date) {
  win <- paste(start_date, end_date, sep = "/")
  
  y_p <- y_H[win]
  stress_p <- stress[win]
  stress_thr_p <- stress_thr[win]
  vol_p <- vol[win]
  vol_thr_p <- vol_thr[win]
  corr_p <- corr_score[win]
  corr_thr_p <- corr_thr[win]
  
  rbind(
    compute_auc_metrics(y_p, stress_p, stress_thr_p, paste0(period_name, ": Residual stress")),
    compute_auc_metrics(y_p, vol_p, vol_thr_p, paste0(period_name, ": Realized volatility")),
    compute_auc_metrics(y_p, corr_p, corr_thr_p, paste0(period_name, ": Pairwise correlation"))
  )
}

sub_tbl <- do.call(rbind, lapply(names(CFG$subperiods), function(nm) {
  pp <- CFG$subperiods[[nm]]
  run_period_metrics(nm, pp[1], pp[2])
}))
save_table(sub_tbl, "Robustness_subperiod_metrics.csv")

cat("STEP 3: Running rolling-PCA robustness checks...\n")

robust_pca_tbl <- data.frame()

for (rw in CFG$pca_rolling_windows) {
  sig_rw <- compute_oos_residual_stress(
    sector_excess,
    K = CFG$K,
    burn_in = CFG$pca_burn_in,
    pca_type = "rolling",
    rolling_window = rw
  )
  
  st_rw <- sig_rw$stress
  th_rw <- rolling_train_quantile(st_rw, CFG$threshold_window, CFG$stress_q)
  
  tmp <- compute_auc_metrics(y_H, st_rw, th_rw, paste0("Rolling PCA ", rw, "d"))
  robust_pca_tbl <- rbind(robust_pca_tbl, tmp)
}

sig_scaled <- compute_oos_residual_stress(
  sector_excess,
  K = CFG$K,
  burn_in = CFG$pca_burn_in,
  pca_type = CFG$pca_type,
  scale_residual = TRUE
)

st_scaled <- sig_scaled$scaled_stress
th_scaled <- rolling_train_quantile(st_scaled, CFG$threshold_window, CFG$stress_q)

robust_pca_tbl <- rbind(
  robust_pca_tbl,
  compute_auc_metrics(y_H, st_scaled, th_scaled, "Idiosyncratic-volatility-scaled stress")
)

save_table(robust_pca_tbl, "Robustness_PCA_window_and_scaled_stress.csv")

# =========================================================
# 10. Residual-ranked sector long-short application
# =========================================================

make_long_short_weights <- function(e_vec) {
  e <- safe_numeric(e_vec)
  names(e) <- names(e_vec)
  
  w_long <- rep(0, length(e))
  names(w_long) <- names(e)
  
  w_short <- rep(0, length(e))
  names(w_short) <- names(e)
  
  pos <- which(is.finite(e) & e > 0)
  neg <- which(is.finite(e) & e < 0)
  
  if (length(pos) > 0) w_long[pos] <- abs(e[pos]) / sum(abs(e[pos]))
  if (length(neg) > 0) w_short[neg] <- abs(e[neg]) / sum(abs(e[neg]))
  
  list(long = w_long, short_abs = w_short)
}

run_residual_ranked_ls <- function(sector_ret_xts, residual_xts, cost_bps = 5) {
  common <- index(na.omit(merge(sector_ret_xts, residual_xts, all = FALSE)))
  R <- sector_ret_xts[common]
  E <- residual_xts[common]
  
  n <- NROW(R)
  ret <- rep(NA_real_, n)
  turnover <- rep(NA_real_, n)
  
  prev_w <- rep(0, NCOL(R))
  names(prev_w) <- colnames(R)
  
  for (i in 2:n) {
    e_lag <- safe_numeric(E[i - 1, ])
    names(e_lag) <- colnames(E)
    
    ww <- make_long_short_weights(e_lag)
    
    w_net <- 0.5 * ww$long - 0.5 * ww$short_abs
    r_i <- safe_numeric(R[i, names(w_net)])
    
    raw_ret <- sum(w_net * r_i, na.rm = TRUE)
    to_i <- sum(abs(w_net - prev_w), na.rm = TRUE)
    
    ret[i] <- raw_ret - to_i * cost_bps / 10000
    turnover[i] <- to_i
    
    prev_w <- w_net
  }
  
  ret[!is.finite(ret)] <- 0
  turnover[!is.finite(turnover)] <- 0
  
  out <- xts(cbind(ret, turnover), order.by = index(R))
  colnames(out) <- c("LS_Return", "Turnover")
  out
}

ls_app <- run_residual_ranked_ls(sector_ret, residuals, CFG$baseline_cost_bps)

ls_stats <- data.frame(
  Strategy = "Residual-ranked sector long-short",
  t(ann_stats(ls_app[, "LS_Return"])),
  Avg_Turnover = mean(safe_numeric(ls_app[, "Turnover"]), na.rm = TRUE),
  stringsAsFactors = FALSE
)
save_table(ls_stats, "Residual_ranked_sector_long_short_application.csv")

plot_ls_application_nav <- function() {
  z <- calc_nav(ls_app[, "LS_Return"])
  z <- na.omit(z)
  
  if (NROW(z) < 2 || all(!is.finite(safe_numeric(z)))) {
    plot_message("Application plot. Residual-ranked sector long-short", "No finite long-short return series is available.")
    return(invisible(NULL))
  }
  
  colnames(z) <- "Residual-ranked long-short"
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  
  plot(
    index(z),
    safe_numeric(z),
    type = "l",
    lwd = 2.0,
    col = COL$stress,
    xlab = "Date",
    ylab = "Cumulative wealth",
    main = "Application plot. Residual-ranked sector long-short"
  )
  
  add_grid()
  
  legend(
    "topright",
    legend = colnames(z),
    col = COL$stress,
    lty = 1,
    lwd = 2.0,
    bty = "o",
    bg = adjustcolor("white", 0.92),
    cex = 0.78
  )
}

plot_ls_application_turnover <- function() {
  z <- na.omit(ls_app[, "Turnover"])
  
  if (NROW(z) < 2 || all(!is.finite(safe_numeric(z)))) {
    plot_message("Application plot. Residual-ranked long-short turnover", "No finite turnover series is available.")
    return(invisible(NULL))
  }
  
  base_par(mar = c(4.2, 4.7, 3.0, 0.7) + 0.1)
  
  plot(
    index(z),
    safe_numeric(z),
    type = "h",
    lwd = 1.0,
    col = adjustcolor(COL$stress, 0.70),
    xlab = "Date",
    ylab = "One-way turnover",
    main = "Application plot. Residual-ranked long-short turnover"
  )
  
  add_grid()
}

# =========================================================
# 11. Reference values and result object
# =========================================================

reference_values <- data.frame(
  Item = c(
    "Table 6 SPY buy-and-hold Ann_Return",
    "Table 6 SPY buy-and-hold Ann_Vol",
    "Table 6 SPY buy-and-hold Sharpe",
    "Table 6 SPY buy-and-hold MaxDD",
    "Table 6 Stress overlay Ann_Return",
    "Table 6 Stress overlay Ann_Vol",
    "Table 6 Stress overlay Sharpe",
    "Table 6 Stress overlay MaxDD",
    "Table 6 Vol overlay Ann_Return",
    "Table 6 Vol overlay Ann_Vol",
    "Table 6 Vol overlay Sharpe",
    "Table 6 Vol overlay MaxDD"
  ),
  Reference_Value = c(
    0.1470, 0.1961, 0.7494, -0.3372,
    0.1447, 0.1827, 0.7924, -0.3006,
    0.1499, 0.1747, 0.8579, -0.2744
  ),
  stringsAsFactors = FALSE
)
save_table(reference_values, "Reference_values_for_checking.csv")

RESULT_TABLES <- list(
  Table_1 = Table_1,
  Table_2A = Table_1[Table_1$Signal %in% c("Residual stress", "SPY realized volatility"), ],
  Table_2B = Table_2B,
  Table_2C = event_diag$overlap,
  Table_2D = event_diag$lead,
  Event_Detail = event_diag$event_detail,
  Table_6 = Table_6,
  Table_7 = cost_sens,
  Robustness_Subperiods = sub_tbl,
  Robustness_PCA = robust_pca_tbl,
  LongShort_Application = ls_stats,
  Reference_Values = reference_values
)

print_result_tables <- function() {
  print_table(Table_1, "Table 1. Early-warning metrics")
  print_table(Table_1[Table_1$Signal %in% c("Residual stress", "SPY realized volatility"), ],
              "Table 2A. Standalone metrics")
  print_table(Table_2B, "Table 2B. Conditional regimes")
  print_table(event_diag$overlap, "Table 2C. Event overlap")
  print_table(event_diag$lead, "Table 2D. Lead time")
  print_table(Table_6, "Table 6. Overlay performance")
  print_table(cost_sens, "Table 7. Transaction-cost sensitivity")
  print_table(sub_tbl, "Robustness. Subperiod metrics")
  print_table(robust_pca_tbl, "Robustness. PCA-window and scaled-stress checks")
  print_table(ls_stats, "Application. Residual-ranked sector long-short")
  invisible(NULL)
}

view_result_tables <- function() {
  if (!interactive()) return(invisible(NULL))
  if (!isTRUE(CFG$show_table_views)) return(invisible(NULL))
  
  for (nm in names(RESULT_TABLES)) {
    try(utils::View(round_df_numeric(RESULT_TABLES[[nm]], 6), title = nm), silent = TRUE)
  }
  
  invisible(NULL)
}

# =========================================================
# 12. Draw all figures directly in RStudio
# =========================================================

show_main_figures <- function() {
  render_figure("Figure_1A_residual_stress_threshold", plot_figure_1A)
  render_figure("Figure_1B_spy_drawdown_onsets", plot_figure_1B)
  render_figure("Figure_2_event_study_residual_stress", plot_figure_2)
  render_figure("Figure_3A_ROC_curves", plot_figure_3A)
  render_figure("Figure_3B_precision_recall_curves", plot_figure_3B)
  render_figure("Figure_4_overlay_wealth", plot_figure_4)
  render_figure("Figure_5_overlay_drawdown", plot_figure_5)
  invisible(NULL)
}

show_diagnostic_figures <- function() {
  render_figure("Diagnostic_event_overlap_bar", plot_event_overlap)
  render_figure("Application_residual_ranked_long_short_wealth", plot_ls_application_nav)
  render_figure("Application_residual_ranked_long_short_turnover", plot_ls_application_turnover)
  invisible(NULL)
}

show_all_figures <- function() {
  show_main_figures()
  if (isTRUE(CFG$show_diagnostic_plots)) show_diagnostic_figures()
  invisible(NULL)
}

# Draw now.
show_all_figures()

if (isTRUE(CFG$show_table_views)) view_result_tables()

cat("\n====================================================\n")
cat("REPLICATION COMPLETE\n")
cat("====================================================\n")
cat("Result tables were printed in the R console and stored in RESULT_TABLES.\n")
cat("Figures were drawn directly in the RStudio Plots pane.\n")
cat("No CSV, PNG, image, or data files were written.\n")
cat("Use the Plots pane left/right arrows to review earlier figures.\n")
cat("Run show_all_figures() to redraw all figures.\n")
cat("Run print_result_tables() to reprint all tables.\n")
cat("Done.\n")


# =========================================================
# 13. Article-complete tables and figures
# =========================================================
# This add-on extends the main replication script so that the GitHub version
# produces all named figures and result tables discussed in the article.
# It keeps the same leakage-safe timing logic. Small numerical differences can
# still occur when Yahoo Finance revises adjusted price histories.

cat("\nSTEP 4: Creating article-complete tables and figures...\n")

# ---------------------------------------------------------
# 13.1 Additional benchmark signals for Table 1
# ---------------------------------------------------------

compute_cross_section_dispersion <- function(excess_xts) {
  X <- as.matrix(excess_xts)
  out <- apply(X, 1, sd, na.rm = TRUE)
  ans <- xts(out, order.by = index(excess_xts))
  colnames(ans) <- "CrossSectionDispersion"
  ans
}

compute_downside_semivol <- function(spy_ret_xts, window = 21) {
  r <- safe_numeric(spy_ret_xts)
  out <- rep(NA_real_, length(r))
  for (i in seq_along(r)) {
    if (i < window) next
    z <- r[(i - window + 1):i]
    z <- z[is.finite(z)]
    z <- pmin(z, 0)
    if (length(z) > 1) out[i] <- sd(z, na.rm = TRUE) * sqrt(252)
  }
  ans <- xts(out, order.by = index(spy_ret_xts))
  colnames(ans) <- "DownsideSemivolatility"
  ans
}

compute_breadth <- function(excess_xts) {
  X <- as.matrix(excess_xts)
  out <- apply(X, 1, function(z) mean(z < 0, na.rm = TRUE))
  ans <- xts(out, order.by = index(excess_xts))
  colnames(ans) <- "Breadth"
  ans
}

download_vix_score <- function(start_date, end_date) {
  vix <- tryCatch({
    suppressWarnings(getSymbols("^VIX", src = "yahoo",
                                from = start_date,
                                to = end_date + 1,
                                auto.assign = FALSE))
  }, error = function(e) NULL)
  
  if (is.null(vix)) {
    ans <- xts(rep(NA_real_, NROW(spy_ret)), order.by = index(spy_ret))
    colnames(ans) <- "VIX"
    return(ans)
  }
  
  v <- Cl(vix)
  colnames(v) <- "VIX"
  v[paste(start_date, end_date, sep = "/")]
}

dispersion_score <- compute_cross_section_dispersion(sector_excess)
dispersion_thr <- rolling_train_quantile(dispersion_score, CFG$threshold_window, CFG$stress_q)

downside_semivol <- compute_downside_semivol(spy_ret, CFG$vol_window)
downside_semivol_thr <- rolling_train_quantile(downside_semivol, CFG$threshold_window, CFG$vol_q)

breadth_score <- compute_breadth(sector_excess)
breadth_thr <- rolling_train_quantile(breadth_score, CFG$threshold_window, CFG$stress_q)

vix_score <- download_vix_score(CFG$start, CFG$end)
vix_thr <- rolling_train_quantile(vix_score, CFG$threshold_window, CFG$vol_q)

dispersion_metrics <- compute_auc_metrics(y_H, dispersion_score, dispersion_thr, "Cross-sectional dispersion")
downside_metrics <- compute_auc_metrics(y_H, downside_semivol, downside_semivol_thr, "Downside semivolatility")
breadth_metrics <- compute_auc_metrics(y_H, breadth_score, breadth_thr, "Breadth")
vix_metrics <- compute_auc_metrics(y_H, vix_score, vix_thr, "VIX")

Table_1 <- rbind(
  stress_metrics,
  vol_metrics,
  corr_metrics,
  dispersion_metrics,
  downside_metrics,
  breadth_metrics,
  vix_metrics
)

Table_1_Full <- data.frame(
  Signal = Table_1$Signal,
  ROC_AUC = Table_1$ROC_AUC,
  PR_AUC = Table_1$PR_AUC,
  Precision = Table_1$Precision,
  Recall = Table_1$Recall,
  F1 = Table_1$F1,
  Event_Rate = Table_1$Positive_Rate,
  stringsAsFactors = FALSE
)

save_table(Table_1_Full, "Table_1_full_early_warning_metrics.csv")

# ---------------------------------------------------------
# 13.2 Combined Figure 1 and Figure 3
# ---------------------------------------------------------

plot_figure_1_combined <- function() {
  z <- na.omit(merge(stress, stress_thr, spy_dd, onset, all = FALSE))
  if (NROW(z) < 2) {
    plot_message("Figure 1. Residual stress and SPY drawdown",
                 "No aligned data are available.")
    return(invisible(NULL))
  }
  
  dates <- index(z)
  par(mfrow = c(2, 1), mar = c(2.9, 4.6, 2.4, 0.8) + 0.1,
      oma = c(2.2, 0, 0.2, 0), las = 1, bty = "l",
      tck = -0.015, mgp = c(2.2, 0.70, 0),
      cex.axis = 0.78, cex.lab = 0.90, cex.main = 0.90)
  
  plot(dates, safe_numeric(z$ResidualStress), type = "l",
       lwd = 1.0, col = COL$drawdown, xaxt = "n",
       xlab = "", ylab = "Residual Stress",
       main = "Residual Stress vs SPY Drawdown",
       ylim = finite_range(c(safe_numeric(z$ResidualStress), safe_numeric(z[, 2])), 0.10))
  lines(dates, safe_numeric(z[, 2]), lwd = 1.1, lty = 2, col = COL$corr)
  onset_dates <- dates[as.logical(z$DrawdownOnset)]
  if (length(onset_dates) > 0) rug(onset_dates, col = COL$onset, lwd = 1.2, ticksize = 0.06)
  add_grid()
  legend("topright",
         legend = c("Stress score", "OOS threshold", "Drawdown onset"),
         col = c(COL$drawdown, COL$corr, COL$onset),
         lty = c(1, 2, NA), pch = c(NA, NA, "|"),
         lwd = c(1.0, 1.1, 1.2), bty = "n", cex = 0.70)
  
  plot(dates, safe_numeric(z$SPY_Drawdown), type = "l",
       lwd = 1.0, col = COL$overlay_stress, xaxt = "n",
       xlab = "", ylab = "SPY Drawdown",
       ylim = finite_range(c(safe_numeric(z$SPY_Drawdown), CFG$drawdown_threshold), 0.10))
  abline(h = CFG$drawdown_threshold, lwd = 1.0, lty = 2, col = COL$corr)
  if (length(onset_dates) > 0) rug(onset_dates, col = COL$onset, lwd = 1.2, ticksize = 0.06)
  nice_date_axis(dates)
  add_grid()
  mtext("Date", side = 1, outer = TRUE, line = 0.7)
}

plot_figure_3_combined <- function() {
  z1 <- na.omit(merge(y_H, stress, all = FALSE))
  if (NROW(z1) < 10) {
    plot_message("Figure 3. Early-warning performance",
                 "Not enough aligned observations.")
    return(invisible(NULL))
  }
  
  y <- safe_numeric(z1[, 1])
  s <- safe_numeric(z1[, 2])
  ok <- is.finite(y) & is.finite(s)
  
  if (sum(ok) < 10 || length(unique(y[ok])) < 2) {
    plot_message("Figure 3. Early-warning performance",
                 "ROC and PR curves are not defined because the label has one class.")
    return(invisible(NULL))
  }
  
  roc_obj <- pROC::roc(y[ok], s[ok], quiet = TRUE, direction = "<")
  fg <- s[ok][y[ok] == 1]
  bg <- s[ok][y[ok] == 0]
  pr <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
  base_rate <- mean(y[ok] == 1, na.rm = TRUE)
  
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3.0, 0.8) + 0.1,
      oma = c(0, 0, 0, 0), las = 1, bty = "l",
      tck = -0.015, mgp = c(2.3, 0.75, 0),
      cex.axis = 0.82, cex.lab = 0.92, cex.main = 0.95)
  
  plot(1 - roc_obj$specificities, roc_obj$sensitivities, type = "l",
       lwd = 2.0, col = COL$corr, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = paste0("ROC Curve (H=", CFG$warning_horizon, ")"))
  abline(0, 1, lty = 3, col = "gray60")
  add_grid()
  
  plot(pr$curve[, 1], pr$curve[, 2], type = "l",
       lwd = 2.0, col = COL$stress, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Recall", ylab = "Precision",
       main = paste0("Precision--Recall Curve (H=", CFG$warning_horizon, ")"))
  abline(h = base_rate, lty = 3, col = "gray60")
  add_grid()
}

# ---------------------------------------------------------
# 13.3 Table 3: subperiod stability in article format
# ---------------------------------------------------------

count_onsets <- function(onset_xts) {
  sum(as.logical(onset_xts), na.rm = TRUE)
}

make_table3_row <- function(sample_name, start_date, end_date) {
  win <- paste(start_date, end_date, sep = "/")
  onsets_n <- count_onsets(onset[win])
  
  s_met <- compute_auc_metrics(y_H[win], stress[win], stress_thr[win], "Residual stress")
  v_met <- compute_auc_metrics(y_H[win], vol[win], vol_thr[win], "SPY realized volatility")
  
  data.frame(
    Sample_Window = sample_name,
    Onsets = onsets_n,
    ROC_AUC_S = ifelse(onsets_n <= 1, NA_real_, s_met$ROC_AUC),
    ROC_AUC_V = ifelse(onsets_n <= 1, NA_real_, v_met$ROC_AUC),
    PR_AUC_S = ifelse(onsets_n <= 1, NA_real_, s_met$PR_AUC),
    PR_AUC_V = ifelse(onsets_n <= 1, NA_real_, v_met$PR_AUC),
    Main_Interpretation = ifelse(onsets_n <= 1,
                                 "Too few onsets for reliable classification",
                                 "Vol stronger; stress complementary"),
    stringsAsFactors = FALSE
  )
}

Table_3 <- do.call(rbind, lapply(names(CFG$subperiods), function(nm) {
  pp <- CFG$subperiods[[nm]]
  label <- switch(nm,
                  Main_2015_2025 = "2015-2025",
                  PrePandemic_2015_2019 = "2015-2019",
                  StressHeavy_2020_2025 = "2020-2025",
                  nm)
  make_table3_row(label, pp[1], pp[2])
}))

# Put article order: 2015-2019, 2020-2025, 2015-2025
Table_3 <- Table_3[match(c("2015-2019", "2020-2025", "2015-2025"), Table_3$Sample_Window), ]
rownames(Table_3) <- NULL
save_table(Table_3, "Table_3_subperiod_stability.csv")

# ---------------------------------------------------------
# 13.4 Table 4: PCA-window, selected-K, and scaled-stress robustness
# ---------------------------------------------------------

q10_minus_q00 <- function(st_xts, th_xts) {
  tmp <- conditional_regime_table(y_H, st_xts, th_xts, vol, vol_thr)
  safe_numeric(tmp$Q10_minus_Q00[1])
}

make_table4_row <- function(spec_name, st_xts, th_xts) {
  met <- compute_auc_metrics(y_H, st_xts, th_xts, spec_name)
  data.frame(
    Specification = spec_name,
    ROC_AUC = met$ROC_AUC,
    PR_AUC = met$PR_AUC,
    Precision = met$Precision,
    Recall = met$Recall,
    Q10_minus_Q00 = q10_minus_q00(st_xts, th_xts),
    Interpretation = "Complementary",
    stringsAsFactors = FALSE
  )
}

select_k_bai_ng_simple <- function(X_centered, max_k = 8) {
  X_centered <- as.matrix(X_centered)
  n <- nrow(X_centered)
  p <- ncol(X_centered)
  max_k <- min(max_k, p - 1, n - 2)
  if (max_k < 1) return(1)
  
  sv <- tryCatch(svd(X_centered, nu = 0, nv = max_k), error = function(e) NULL)
  if (is.null(sv)) return(1)
  
  vals <- sv$d^2
  total_ss <- sum(X_centered^2, na.rm = TRUE)
  penalty <- ((n + p) / (n * p)) * log((n * p) / (n + p))
  
  ic <- rep(NA_real_, max_k)
  for (k in 1:max_k) {
    resid_ss <- max(total_ss - sum(vals[1:k], na.rm = TRUE), .Machine$double.eps)
    V_k <- resid_ss / (n * p)
    ic[k] <- log(V_k) + k * penalty
  }
  
  which.min(ic)
}

compute_oos_residual_stress_selectedK <- function(X_xts, burn_in = 252, max_k = 8) {
  X <- as.matrix(X_xts)
  dates <- index(X_xts)
  n <- nrow(X)
  p <- ncol(X)
  
  stress <- rep(NA_real_, n)
  selected_k <- rep(NA_integer_, n)
  
  for (i in seq_len(n)) {
    if (i <= burn_in) next
    
    train_idx <- 1:(i - 1)
    X_train <- X[train_idx, , drop = FALSE]
    X_train <- X_train[complete.cases(X_train), , drop = FALSE]
    x_t <- X[i, ]
    
    if (nrow(X_train) < max(30, max_k + 5) || any(!is.finite(x_t))) next
    
    mu <- colMeans(X_train, na.rm = TRUE)
    X_centered <- sweep(X_train, 2, mu, "-")
    K_use <- select_k_bai_ng_simple(X_centered, max_k = max_k)
    selected_k[i] <- K_use
    
    pca_fit <- tryCatch(prcomp(X_centered, center = FALSE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pca_fit)) next
    
    K_use <- min(K_use, ncol(pca_fit$rotation))
    L <- pca_fit$rotation[, 1:K_use, drop = FALSE]
    xt_centered <- x_t - mu
    fitted <- L %*% (t(L) %*% xt_centered)
    e_t <- as.numeric(xt_centered - fitted)
    stress[i] <- sqrt(mean(e_t^2, na.rm = TRUE))
  }
  
  stress_xts <- xts(stress, order.by = dates)
  colnames(stress_xts) <- "BaiNgSelectedKStress"
  
  k_xts <- xts(selected_k, order.by = dates)
  colnames(k_xts) <- "SelectedK"
  
  list(stress = stress_xts, selected_k = k_xts)
}

cat("Running selected-K robustness check for Table 4...\n")
selected_k_sig <- compute_oos_residual_stress_selectedK(
  sector_excess,
  burn_in = CFG$pca_burn_in,
  max_k = 8
)
stress_selected_k <- selected_k_sig$stress
stress_selected_k_thr <- rolling_train_quantile(stress_selected_k, CFG$threshold_window, CFG$stress_q)

# Existing rolling and scaled objects are recomputed here so that Table 4 is self-contained.
sig_roll_252 <- compute_oos_residual_stress(sector_excess, K = CFG$K, burn_in = CFG$pca_burn_in,
                                            pca_type = "rolling", rolling_window = 252)
st_roll_252 <- sig_roll_252$stress
th_roll_252 <- rolling_train_quantile(st_roll_252, CFG$threshold_window, CFG$stress_q)

sig_roll_504 <- compute_oos_residual_stress(sector_excess, K = CFG$K, burn_in = CFG$pca_burn_in,
                                            pca_type = "rolling", rolling_window = 504)
st_roll_504 <- sig_roll_504$stress
th_roll_504 <- rolling_train_quantile(st_roll_504, CFG$threshold_window, CFG$stress_q)

sig_scaled2 <- compute_oos_residual_stress(sector_excess, K = CFG$K, burn_in = CFG$pca_burn_in,
                                           pca_type = CFG$pca_type, scale_residual = TRUE)
st_scaled2 <- sig_scaled2$scaled_stress
th_scaled2 <- rolling_train_quantile(st_scaled2, CFG$threshold_window, CFG$stress_q)

Table_4 <- rbind(
  make_table4_row("Expanding, K = 2, unscaled RMS", stress, stress_thr),
  make_table4_row("Rolling 252, K = 2", st_roll_252, th_roll_252),
  make_table4_row("Rolling 504, K = 2", st_roll_504, th_roll_504),
  make_table4_row("Expanding, Bai-Ng selected K", stress_selected_k, stress_selected_k_thr),
  make_table4_row("Expanding, K = 2, idio-vol scaled", st_scaled2, th_scaled2)
)
save_table(Table_4, "Table_4_PCA_window_factor_selection_scaled_stress.csv")

Selected_K_Summary <- data.frame(
  Selected_K = as.integer(names(table(na.omit(as.numeric(selected_k_sig$selected_k))))),
  Count = as.integer(table(na.omit(as.numeric(selected_k_sig$selected_k)))),
  stringsAsFactors = FALSE
)
save_table(Selected_K_Summary, "Selected_K_summary.csv")

# ---------------------------------------------------------
# 13.5 Table 5: broader-sample proxy evidence
# ---------------------------------------------------------
# The long-history proxy uses sector ETFs with histories reaching back to 1999.
# XLC and XLRE are omitted because they do not have the same long history.
# If the article was produced with a different proxy file, replace
# proxy_sector_tickers below with that exact proxy universe.

run_proxy_pipeline <- function() {
  proxy_start <- as.Date("1999-01-01")
  proxy_end <- CFG$end
  proxy_spy <- "SPY"
  proxy_sector_tickers <- c("XLY", "XLP", "XLE", "XLF", "XLV", "XLI", "XLK", "XLB", "XLU")
  proxy_tickers <- c(proxy_spy, proxy_sector_tickers)
  
  px_raw <- tryCatch({
    invisible(getSymbols(Symbols = proxy_tickers, src = "yahoo",
                         from = proxy_start, to = proxy_end + 1,
                         auto.assign = TRUE, warnings = FALSE))
    do.call(merge, lapply(proxy_tickers, function(tk) Ad(get(tk))))
  }, error = function(e) NULL)
  
  if (is.null(px_raw)) return(NULL)
  
  colnames(px_raw) <- proxy_tickers
  px <- px_raw[paste(proxy_start, proxy_end, sep = "/")]
  px <- px[complete.cases(px), ]
  rr <- na.omit(ROC(px, type = "discrete"))
  colnames(rr) <- colnames(px)
  
  spy_r <- rr[, proxy_spy, drop = FALSE]
  sec_r <- rr[, proxy_sector_tickers]
  sec_x <- sec_r - as.numeric(spy_r)
  colnames(sec_x) <- proxy_sector_tickers
  
  sig <- compute_oos_residual_stress(sec_x, K = CFG$K, burn_in = CFG$pca_burn_in,
                                     pca_type = "expanding")
  st <- sig$stress
  th <- rolling_train_quantile(st, CFG$threshold_window, CFG$stress_q)
  
  vv <- compute_realized_vol(spy_r, CFG$vol_window)
  vth <- rolling_train_quantile(vv, CFG$threshold_window, CFG$vol_q)
  
  nav <- calc_nav(spy_r)
  dd <- nav / cummax(nav) - 1
  colnames(dd) <- "SPY_Drawdown"
  
  dd_reg <- dd <= CFG$drawdown_threshold
  dd_lag <- lag(dd_reg, 1)
  dd_lag[is.na(dd_lag)] <- FALSE
  ons <- dd_reg & !dd_lag
  colnames(ons) <- "DrawdownOnset"
  
  yy <- make_horizon_label(ons, CFG$warning_horizon)
  
  sm <- compute_auc_metrics(yy, st, th, "Proxy residual stress")
  vm <- compute_auc_metrics(yy, vv, vth, "Proxy SPY realized volatility")
  qtbl <- conditional_regime_table(yy, st, th, vv, vth)
  
  ed <- event_alarm_diagnostics(ons, st, th, vv, vth, CFG$event_alarm_lookback)
  s_only <- ed$overlap$Count[ed$overlap$Category == "Stress only"]
  paired <- ed$event_detail[ed$event_detail$Stress_Flag & ed$event_detail$Vol_Flag, ]
  p_earlier <- ifelse(nrow(paired) > 0,
                      mean(paired$Stress_Lead_Days > paired$Vol_Lead_Days, na.rm = TRUE),
                      NA_real_)
  
  list(
    onsets = sum(as.logical(ons), na.rm = TRUE),
    stress_metrics = sm,
    vol_metrics = vm,
    q10_minus_q00 = qtbl$Q10_minus_Q00[1],
    stress_only = s_only,
    p_stress_earlier = p_earlier,
    proxy_detail = list(stress = st, vol = vv, onset = ons, y = yy,
                        event_diag = ed, tickers = proxy_sector_tickers)
  )
}

proxy_results <- tryCatch(run_proxy_pipeline(), error = function(e) NULL)

if (!is.null(proxy_results)) {
  Table_5 <- data.frame(
    Sample = c("Modern ETF sample", "Long-history proxy sample"),
    Onsets = c(count_onsets(onset), proxy_results$onsets),
    ROC_AUC_S = c(stress_metrics$ROC_AUC, proxy_results$stress_metrics$ROC_AUC),
    ROC_AUC_V = c(vol_metrics$ROC_AUC, proxy_results$vol_metrics$ROC_AUC),
    PR_AUC_S = c(stress_metrics$PR_AUC, proxy_results$stress_metrics$PR_AUC),
    PR_AUC_V = c(vol_metrics$PR_AUC, proxy_results$vol_metrics$PR_AUC),
    Q10_minus_Q00 = c(Table_2B$Q10_minus_Q00[1], proxy_results$q10_minus_q00),
    Stress_Only = c(event_diag$overlap$Count[event_diag$overlap$Category == "Stress only"],
                    proxy_results$stress_only),
    Pr_Stress_Earlier = c(event_diag$lead$Value[event_diag$lead$Metric == "Stress earlier share"],
                          proxy_results$p_stress_earlier),
    stringsAsFactors = FALSE
  )
} else {
  Table_5 <- data.frame(
    Sample = "Long-history proxy sample",
    Note = "Proxy sample could not be downloaded or computed in this run.",
    stringsAsFactors = FALSE
  )
}

save_table(Table_5, "Table_5_baseline_and_broader_sample_evidence.csv")

# ---------------------------------------------------------
# 13.6 Transaction-cost table with state-dependent cost row
# ---------------------------------------------------------

run_spy_overlay_state_dependent_cost <- function(spy_ret_xts, trigger_xts,
                                                 lambda = 0.30,
                                                 normal_cost_bps = 5,
                                                 stress_cost_bps = 25) {
  z <- merge(spy_ret_xts, trigger_xts, all = FALSE)
  r <- safe_numeric(z[, 1])
  trig <- as.logical(z[, 2])
  trig[is.na(trig)] <- FALSE
  
  trig_lag <- c(FALSE, head(trig, -1))
  exposure <- ifelse(trig_lag, 1 - lambda, 1.0)
  exposure_lag <- c(1.0, head(exposure, -1))
  turnover <- abs(exposure - exposure_lag)
  
  cost_bps <- ifelse(trig_lag, stress_cost_bps, normal_cost_bps)
  cost <- turnover * cost_bps / 10000
  net_ret <- exposure * r - cost
  
  out <- xts(cbind(net_ret, exposure, turnover, cost), order.by = index(z))
  colnames(out) <- c("Net_Return", "Exposure", "Turnover", "Cost")
  out
}

state_cost_ret <- run_spy_overlay_state_dependent_cost(
  spy_ret, stress_regime, CFG$overlay_lambda, normal_cost_bps = 5, stress_cost_bps = 25
)[, "Net_Return"]

state_cost_stats <- ann_stats(state_cost_ret)

Table_7 <- rbind(
  cost_sens,
  data.frame(Cost_bps = "5 normal / 25 stress",
             Ann_Return = state_cost_stats["Ann_Return"],
             Ann_Vol = state_cost_stats["Ann_Vol"],
             Sharpe = state_cost_stats["Sharpe"],
             MaxDD = state_cost_stats["MaxDD"],
             stringsAsFactors = FALSE)
)
save_table(Table_7, "Table_7_transaction_cost_sensitivity_with_state_dependent_cost.csv")


# ---------------------------------------------------------
# 13.7 Additional PCA and strategy plots used in the article/supporting materials
# ---------------------------------------------------------
# These plots correspond to the PCA diagnostics and the strategy-comparison
# figures shown in the article/supporting material screenshots: scree plot,
# PC1--PC2 loading heatmap, PC1--PC2 factor-score time series, and the
# cumulative-wealth/drawdown figures including the residual-ranked long-short
# application. They are computed from the same sector-excess-return panel and
# the same strategy objects already created above.

if (is.null(CFG$save_figures_png)) CFG$save_figures_png <- TRUE
if (is.null(CFG$figures_dir)) CFG$figures_dir <- "figures"

get_full_sample_pca <- function() {
  X <- na.omit(sector_excess)
  if (NROW(X) < 20) return(NULL)
  pca_obj <- tryCatch(
    prcomp(coredata(X), center = TRUE, scale. = TRUE),
    error = function(e) NULL
  )
  if (is.null(pca_obj)) return(NULL)
  list(pca = pca_obj, X = X)
}

plot_pca_scree <- function() {
  obj <- get_full_sample_pca()
  if (is.null(obj)) {
    plot_message("PCA Scree Plot", "Not enough valid sector excess return observations.")
    return(invisible(NULL))
  }
  
  pca_obj <- obj$pca
  prop_var <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  pcs <- seq_along(prop_var)
  
  base_par(mar = c(4.4, 4.9, 3.2, 0.9) + 0.1)
  plot(
    pcs, prop_var,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "#0072B2",
    xlab = "Principal Component",
    ylab = "Proportion Variance Explained",
    main = "PCA Scree Plot (Sector Excess Returns)"
  )
  abline(v = CFG$K, col = "#CC6677", lty = 2, lwd = 1.2)
  add_grid()
}

plot_pca_loadings_heatmap <- function() {
  obj <- get_full_sample_pca()
  if (is.null(obj)) {
    plot_message("PCA Loadings Heatmap", "Not enough valid sector excess return observations.")
    return(invisible(NULL))
  }
  
  pca_obj <- obj$pca
  
  loadings_mat <- pca_obj$rotation[, 1:2, drop = FALSE]
  colnames(loadings_mat) <- c("PC1", "PC2")
  
  sectors <- rownames(loadings_mat)
  pcs <- colnames(loadings_mat)
  
  vals <- as.numeric(loadings_mat)
  zlim <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(zlim) || zlim <= 0) zlim <- 1
  
  pal <- colorRampPalette(c("#2c7bb6", "#f7f7f7", "#d95f02"))(101)
  val_to_col <- function(v) {
    idx <- round((v + zlim) / (2 * zlim) * 100) + 1
    idx <- pmax(1, pmin(101, idx))
    pal[idx]
  }
  
  base_par(mar = c(6.6, 4.5, 3.1, 1.0) + 0.1)
  
  plot(
    NA,
    xlim = c(0.5, length(sectors) + 0.5),
    ylim = c(0.5, length(pcs) + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = "PCA Loadings Heatmap (PC1 and PC2)",
    bty = "n"
  )
  
  for (j in seq_along(sectors)) {
    for (i in seq_along(pcs)) {
      v <- loadings_mat[j, pcs[i]]
      rect(
        xleft = j - 0.5,
        ybottom = i - 0.5,
        xright = j + 0.5,
        ytop = i + 0.5,
        col = val_to_col(v),
        border = "white"
      )
    }
  }
  
  axis(1, at = seq_along(sectors), labels = sectors, las = 2)
  axis(2, at = seq_along(pcs), labels = pcs, las = 1)
  box()
}

plot_pca_factor_scores <- function() {
  obj <- get_full_sample_pca()
  if (is.null(obj)) {
    plot_message("PCA Factor Scores", "Not enough valid sector excess return observations.")
    return(invisible(NULL))
  }
  
  pca_obj <- obj$pca
  X <- obj$X
  scores_xts <- xts(pca_obj$x[, 1:2, drop = FALSE], order.by = index(X))
  colnames(scores_xts) <- c("PC1", "PC2")
  
  dates <- index(scores_xts)
  y1 <- safe_numeric(scores_xts$PC1)
  y2 <- safe_numeric(scores_xts$PC2)
  
  base_par(mar = c(4.4, 4.9, 3.2, 0.9) + 0.1)
  plot(
    dates, y1,
    type = "l",
    lwd = 1.45,
    col = "#CC79A7",
    xaxt = "n",
    xlab = "Date",
    ylab = "Score",
    main = "PCA Factor Scores (PC1 and PC2)",
    ylim = finite_range(c(y1, y2), 0.10)
  )
  nice_date_axis(dates)
  lines(dates, y2, lwd = 1.45, col = "#00A6D6")
  abline(h = 0, col = "#0072B2", lty = 3)
  
  onset_aligned <- onset[index(onset) %in% dates]
  onset_dates <- index(onset_aligned)[as.logical(onset_aligned[, 1])]
  if (length(onset_dates) > 0) {
    abline(v = onset_dates, col = adjustcolor("orange", 0.16), lwd = 0.7)
  }
  
  add_grid()
  legend(
    "topleft",
    legend = c("PC1", "PC2"),
    col = c("#CC79A7", "#00A6D6"),
    lty = 1,
    lwd = 1.8,
    bty = "n",
    cex = 0.85
  )
}

# Override Figure 4 and Figure 5 so they include the residual-ranked long-short
# application line, matching the strategy-comparison plots in the article/supporting
# material screenshots. If ls_app is not available, the functions fall back to the
# three SPY-based series.
plot_figure_4 <- function() {
  if (exists("ls_app") && NROW(ls_app) > 2) {
    z <- na.omit(merge(
      calc_nav(spy_bh_ret),
      calc_nav(stress_overlay_ret),
      calc_nav(vol_overlay_ret),
      calc_nav(ls_app[, "LS_Return"]),
      all = FALSE
    ))
    labels <- c("SPY_Return", "SPY_Overlay_Resid_ex_0.3", "SPY_Overlay_SPYVol_ex_0.3", "APT_LS_Excess")
    cols <- c("black", "#E69F00", "#009E73", "#0072B2")
  } else {
    z <- na.omit(merge(calc_nav(spy_bh_ret), calc_nav(stress_overlay_ret), calc_nav(vol_overlay_ret), all = FALSE))
    labels <- c("SPY_Return", "SPY_Overlay_Resid_ex_0.3", "SPY_Overlay_SPYVol_ex_0.3")
    cols <- c("black", "#E69F00", "#009E73")
  }
  
  if (NROW(z) < 2) {
    plot_message("Figure 4. Cumulative wealth", "No aligned strategy return data.")
    return(invisible(NULL))
  }
  
  colnames(z) <- labels
  base_par(mar = c(4.4, 4.8, 3.1, 0.8) + 0.1)
  matplot(
    index(z),
    coredata(z),
    type = "l",
    lty = 1,
    lwd = rep(1.8, NCOL(z)),
    col = cols,
    xaxt = "n",
    xlab = "Date",
    ylab = "Cumulative Wealth",
    main = "Cumulative Wealth (Strategies)"
  )
  nice_date_axis(index(z))
  abline(h = 1, col = "#7fb3d5", lty = 3)
  add_grid()
  legend("topleft", legend = colnames(z), col = cols, lty = 1, lwd = 1.8, bty = "n", cex = 0.78)
}

plot_figure_5 <- function() {
  if (exists("ls_app") && NROW(ls_app) > 2) {
    z <- na.omit(merge(
      calc_drawdown(spy_bh_ret),
      calc_drawdown(stress_overlay_ret),
      calc_drawdown(vol_overlay_ret),
      calc_drawdown(ls_app[, "LS_Return"]),
      all = FALSE
    ))
    labels <- c("SPY_Return", "SPY_Overlay_Resid_ex_0.3", "SPY_Overlay_SPYVol_ex_0.3", "APT_LS_Excess")
    cols <- c("black", "#E69F00", "#009E73", "#0072B2")
  } else {
    z <- na.omit(merge(calc_drawdown(spy_bh_ret), calc_drawdown(stress_overlay_ret), calc_drawdown(vol_overlay_ret), all = FALSE))
    labels <- c("SPY_Return", "SPY_Overlay_Resid_ex_0.3", "SPY_Overlay_SPYVol_ex_0.3")
    cols <- c("black", "#E69F00", "#009E73")
  }
  
  if (NROW(z) < 2) {
    plot_message("Figure 5. Drawdowns", "No aligned drawdown data.")
    return(invisible(NULL))
  }
  
  colnames(z) <- labels
  base_par(mar = c(4.4, 4.8, 3.1, 0.8) + 0.1)
  matplot(
    index(z),
    coredata(z),
    type = "l",
    lty = 1,
    lwd = rep(1.8, NCOL(z)),
    col = cols,
    xaxt = "n",
    xlab = "Date",
    ylab = "Drawdown",
    main = "Drawdowns (Strategies)"
  )
  nice_date_axis(index(z))
  abline(h = 0, col = "#7fb3d5", lty = 3)
  add_grid()
  legend("bottomleft", legend = colnames(z), col = cols, lty = 1, lwd = 1.8, bty = "n", cex = 0.78)
}

show_pca_figures <- function() {
  render_figure("Figure_6_PCA_scree_plot", plot_pca_scree)
  render_figure("Figure_7_PCA_loadings_heatmap", plot_pca_loadings_heatmap)
  render_figure("Figure_8_PCA_factor_scores_PC1_PC2", plot_pca_factor_scores)
  invisible(NULL)
}

save_plot_png <- function(filename, plot_fun, width = 1400, height = 900, res = 150) {
  dir.create(CFG$figures_dir, showWarnings = FALSE, recursive = TRUE)
  png(file.path(CFG$figures_dir, filename), width = width, height = height, res = res)
  tryCatch(
    plot_fun(),
    error = function(e) {
      plot_message(paste("Plot failed:", filename), e$message)
    },
    finally = dev.off()
  )
  invisible(NULL)
}

save_all_figures_png <- function() {
  save_plot_png("Figure_1_combined_residual_stress_drawdown.png", plot_figure_1_combined, width = 1400, height = 1050)
  save_plot_png("Figure_2_event_study_residual_stress.png", plot_figure_2, width = 1400, height = 900)
  save_plot_png("Figure_3_combined_ROC_PR.png", plot_figure_3_combined, width = 1500, height = 800)
  save_plot_png("Figure_4_cumulative_wealth_strategies.png", plot_figure_4, width = 1400, height = 900)
  save_plot_png("Figure_5_drawdowns_strategies.png", plot_figure_5, width = 1400, height = 900)
  save_plot_png("Figure_6_PCA_scree_plot.png", plot_pca_scree, width = 1400, height = 800)
  save_plot_png("Figure_7_PCA_loadings_heatmap.png", plot_pca_loadings_heatmap, width = 1400, height = 700)
  save_plot_png("Figure_8_PCA_factor_scores_PC1_PC2.png", plot_pca_factor_scores, width = 1400, height = 800)
  if (exists("plot_event_overlap")) save_plot_png("Diagnostic_event_overlap_bar.png", plot_event_overlap, width = 1200, height = 800)
  if (exists("plot_ls_application_nav")) save_plot_png("Application_residual_ranked_long_short_wealth.png", plot_ls_application_nav, width = 1400, height = 800)
  if (exists("plot_ls_application_turnover")) save_plot_png("Application_residual_ranked_long_short_turnover.png", plot_ls_application_turnover, width = 1400, height = 800)
  invisible(NULL)
}

# ---------------------------------------------------------
# 13.7 Redraw helpers for all article figures and tables
# ---------------------------------------------------------

show_core_figures <- function() {
  render_figure("Figure_1_article_combined_residual_stress_and_SPY_drawdown", plot_figure_1_combined)
  render_figure("Figure_2_event_study_residual_stress", plot_figure_2)
  render_figure("Figure_3_combined_ROC_PR_residual_stress", plot_figure_3_combined)
  render_figure("Figure_4_cumulative_wealth_strategies", plot_figure_4)
  render_figure("Figure_5_drawdowns_strategies", plot_figure_5)
  show_pca_figures()
  invisible(NULL)
}

print_all_tables <- function() {
  print_table(Table_1_Full, "Table 1. Early-warning metrics across residual stress and benchmark signals")
  print_table(Table_1_Full[Table_1_Full$Signal %in% c("Residual stress", "SPY realized volatility"), ],
              "Table 2 Panel A. Standalone early-warning metrics")
  print_table(Table_2B, "Table 2 Panel B. Conditional onset probability by joint regimes")
  print_table(event_diag$overlap, "Table 2 Panel C. Event-level overlap")
  print_table(event_diag$lead, "Table 2 Panel D. Lead time to onset for paired events")
  print_table(Table_3, "Table 3. Subperiod stability")
  print_table(Table_4, "Table 4. Robustness to PCA window, factor selection, and residual scaling")
  print_table(Table_5, "Table 5. Baseline and broader-sample evidence")
  print_table(Table_6, "Table 6. Overlay performance")
  print_table(Table_7, "Table 7. Transaction-cost sensitivity")
  invisible(NULL)
}

RESULT_TABLES <- c(
  RESULT_TABLES,
  list(
    Table_1_Full = Table_1_Full,
    Table_3 = Table_3,
    Table_4 = Table_4,
    Selected_K_Summary = Selected_K_Summary,
    Table_5 = Table_5,
    Table_7_StateDependent = Table_7
  )
)

cat("\n====================================================\n")
cat("FULL OUTPUT READY\n")
cat("====================================================\n")
cat("The script now prints all named combined tables and displays all named article figures.\n")
cat("Run show_core_figures() to redraw Figures 1-5 in combined format.\n")
cat("Run print_all_tables() to reprint Tables 1-7 in combined format.\n")
cat(".\n")
cat("Done.\n")

# Draw the combined figures once at the end.
show_core_figures()
if (isTRUE(CFG$save_figures_png)) {
  save_all_figures_png()
  cat("Saved PNG figures to: ", CFG$figures_dir, "\n", sep = "")
}
print_all_tables()
