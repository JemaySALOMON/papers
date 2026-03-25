# install.R

if (!requireNamespace("svglite", quietly = TRUE)){
  install.packages("svglite")}
#Bio-conductor
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("BiocVersion")

pkgs  <- c(
  "MASS", "dplyr", "gtools", "ggplot2", "future.apply",
  "mvtnorm", "MixMatrix", "Matrix", "seriation", "TMB",
  "corpcor", "data.table", "doParallel", "agricolae",
  "grid", "tidyr", "ggExtra", "here", "scrm",
  "plotly", "data.table", "lme4", "car", "remotes",
  "devtools", "ggpubr")
pkgs <- unique(pkgs)

# Install packages and handle errors
if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}
lapply(pkgs, if_missing)

if (!requireNamespace("rutilstimflutre", quietly = TRUE)){
  devtools::install_github("timflutre/rutilstimflutre")}

