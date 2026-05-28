# Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring

This repository provides the R replication code for the article:

**Beyond Volatility: A Leakage-Safe Residual-Stress Signal for Drawdown Risk Monitoring**

## Overview

The code implements a leakage-safe residual-stress framework for drawdown-risk monitoring using SPY and 11 U.S. sector ETFs.

The workflow includes:

- Downloading adjusted price data from Yahoo Finance
- Constructing sector excess returns relative to SPY
- Estimating a leakage-safe PCA residual-stress signal
- Computing train-only rolling quantile thresholds
- Constructing drawdown-onset and early-warning labels
- Evaluating ROC-AUC and PR-AUC performance
- Conducting conditional regime analysis
- Computing event-overlap and lead-time diagnostics
- Evaluating stress-managed SPY overlay strategies
- Running robustness checks

## Asset Universe

The analysis uses SPY and the following 11 U.S. sector ETFs:

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

The main script automatically checks for missing packages and installs them when needed.

## How to Run

Open RStudio, set the working directory to this repository folder, and run:

```r
rm(list = ls())
graphics.off()
source("replication_code.R")
