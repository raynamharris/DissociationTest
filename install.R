# script to install R packages needed for pipeline
# must be run from analysis_pipeline directory
# argument to script is path to R library for install

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) .libPaths(args[1])


install.packages(c("highr", "caTools", "knitr", "rmarkdown", "dplyr", "plyr", "reshape2"), repos=NULL, type="source")