# ezRun
An R meta-package for the analysis of Next Generation Sequencing data

## Installation of the development version of `ezRun` from github

```R
devtools::install_github("uzh/ezRun")
```

## Dependencies of R/Bioconductor packages
```R
packages <- c("testthat", "knitr", "gage", "goseq", "ChIPpeakAnno", 
              "DESeq2", "TEQC", "htSeqTools", "pathview", "reshape2", 
              "vsn", "Rsubread", "preprocessCore", "wesanderson",
              "RCurl", "caTools", "matrixStats", "Repitools", "DT", 
              "htmltools", "biomaRt", "grid", "gridExtra",
              "RColorBrewer", "WGCNA", "plyr", "pvclust", "parallel", 
              "Biostrings", "ReporteRs", "Rsamtools", "Hmisc", "XML", 
              "stringr", "GenomicAlignments", "GenomicFeatures",
              "GenomicRanges", "ShortRead", "Gviz", "gplots", "GO.db", 
              "GOstats", "annotate", "bitops", "edgeR", "limma", "S4Vectors",
              "VariantAnnotation", "rmarkdown", "plotly", "scran",
              "ReporteRsjars", "data.table", "kableExtra", "htmlwidgets",
              "RSelenium", "webshot", "clusterProfiler", "dupRadar",
              "taxize")
packages <- setdiff(packages, rownames(installed.packages()))
biocLite(packages)
```

## Dependencies of external software
* phantomjs
* bwa, bowtie, bowtie2, STAR, picard, sambamba
* lsof