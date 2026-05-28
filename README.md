# Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring

This repository provides the R replication code, result tables, and output figures for the article:

**Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring**

## Overview

The code implements a leakage-safe residual-stress framework for drawdown-risk monitoring using SPY and 11 U.S. sector ETFs.

The workflow includes:

* Downloading adjusted price data from Yahoo Finance
* Constructing sector excess returns relative to SPY
* Estimating a leakage-safe PCA residual-stress signal
* Computing train-only rolling quantile thresholds
* Constructing drawdown-onset and early-warning labels
* Evaluating ROC-AUC and PR-AUC performance
* Comparing residual stress with volatility, correlation, dispersion, downside semivolatility, breadth, and VIX benchmarks
* Conducting conditional regime analysis
* Computing event-overlap and lead-time diagnostics
* Running subperiod and robustness checks
* Evaluating stress-managed SPY overlay strategies
* Producing output figures for the main analysis and diagnostics

## Asset Universe

The main analysis uses SPY and the following 11 U.S. sector ETFs:

SPY, XLC, XLY, XLP, XLE, XLF, XLV, XLI, XLK, XLB, XLU, XLRE.

## Requirements

The code uses the following R packages:

* quantmod
* xts
* zoo
* PerformanceAnalytics
* pROC
* PRROC
* scales

The main script automatically checks for missing packages and installs them when needed.

## How to Run

Open RStudio, set the working directory to this repository folder, and run:

```r
rm(list = ls())
graphics.off()
source("replication_code.R")
```

The script prints result tables in the R console and displays figures directly in the RStudio Plots pane.

The script also saves PNG versions of the generated figures into the `figures/` folder.

## Main Outputs

Running `replication_code.R` produces the main empirical tables and figures used in the analysis.

The code prints the following result tables in the R console:

* Table 1: Early-warning metrics across residual stress and benchmark signals
* Table 2 Panel A: Standalone early-warning metrics
* Table 2 Panel B: Conditional onset probability by joint regimes
* Table 2 Panel C: Event-level overlap
* Table 2 Panel D: Lead-time diagnostics
* Table 3: Subperiod stability
* Table 4: Robustness to PCA window, factor selection, and residual scaling
* Table 5: Baseline and broader-sample evidence
* Table 6: Overlay performance
* Table 7: Transaction-cost sensitivity

The generated figures include:

* Residual stress with train-only threshold
* SPY drawdown and drawdown-onset events
* Event-study of residual stress around drawdown onsets
* ROC and precision-recall curves
* Overlay cumulative wealth
* Overlay drawdown
* PCA scree plot
* PCA loading heatmap
* PCA factor-score plot
* Event-overlap diagnostic plot
* Residual-ranked long-short application plots
* Additional robustness and diagnostic figures

## Figures Folder

The `figures/` folder contains saved PNG output figures generated from the replication code.

These files are included so that readers can directly inspect the main visual outputs without rerunning the full script.

## Data Source

The financial price data are publicly available from Yahoo Finance and are downloaded directly by the R script using `quantmod::getSymbols()`.

Because Yahoo Finance may revise adjusted-price histories over time, small numerical differences may arise when the code is rerun at a later date.

## Reproducibility Notes

The code follows a leakage-safe timing protocol:

* PCA mappings are estimated using information available only through time `t - 1`.
* Residual stress is computed out-of-sample at time `t`.
* Rolling thresholds are computed using train-only historical windows and shifted forward by one trading day.
* Portfolio and overlay applications use lagged signals to avoid look-ahead bias.

The repository is intended to reproduce the main empirical workflow, tables, figures, and qualitative conclusions. Exact numerical values may differ slightly if Yahoo Finance updates adjusted-price histories.

## Repository Files

* `replication_code.R`: Main R replication script
* `requirements.R`: Package installation helper
* `README.md`: Repository description and instructions
* `LICENSE`: Academic-use license
* `figures/`: Saved PNG output figures

## License

This repository is released for academic reproducibility.
