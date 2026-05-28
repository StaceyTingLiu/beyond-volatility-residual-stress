≈
# Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring

This repository provides the R replication code for the article:

**Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring**

## Overview

The code implements a leakage-safe residual-stress framework for drawdown-risk monitoring using SPY and 11 U.S. sector ETFs.

The workflow includes data download from Yahoo Finance, sector excess return construction, leakage-safe PCA residual-stress construction, train-only rolling thresholds, drawdown-onset labels, ROC/PR evaluation, conditional regime analysis, event-overlap and lead-time diagnostics, stress-managed overlay analysis, and robustness checks.

## Asset Universe

SPY, XLC, XLY, XLP, XLE, XLF, XLV, XLI, XLK, XLB, XLU, XLRE.

## Requirements

The code uses the following R packages:

- quantmod
- xts
- zoo
- PerformanceAnalytics
- pROC
- PRROC
- scales

## How to Run

Open RStudio and run:

```r
rm(list = ls())
graphics.off()
source("replication_code.R")
