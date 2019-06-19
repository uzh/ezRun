
setwd("/scratch/gtan/dev/quickdev")
library(ezRun)
input <- "input_dataset.tsv"
input <- EzDataset(file=input, dataRoot="/srv/gstore/projects")

param <- list()
param$name <- "CellRangerAggr"
