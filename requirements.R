packages <- c(
  "quantmod",
  "xts",
  "zoo",
  "PerformanceAnalytics",
  "pROC",
  "PRROC",
  "scales"
)

install.packages(setdiff(packages, rownames(installed.packages())))
