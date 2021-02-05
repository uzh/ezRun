# ezRun
An R meta-package for the analysis of Next Generation Sequencing data.

The version of `ezRun` package is bound to the release of Bioconductor development branch.

## Dependencies of Python packages
```Python
pip3 install velocyto magic-impute
pip3 install multiqc
```

## Dependencies of R/Bioconductor packages
```R
packages <- c("testthat", "knitr", "goseq", "ChIPpeakAnno", 
              "DESeq2", "TEQC", "pathview", "reshape2", 
              "vsn", "Rsubread", "preprocessCore", "wesanderson",
              "RCurl", "caTools", "matrixStats", "Repitools", "DT", 
              "htmltools", "biomaRt", "grid", "gridExtra",
              "RColorBrewer", "WGCNA", "plyr", "pvclust", "parallel", 
              "Biostrings", "Rsamtools", "Hmisc", "XML", 
              "stringr", "GenomicAlignments", "GenomicFeatures",
              "GenomicRanges", "ShortRead", "Gviz", "gplots", "GO.db", 
              "GOstats", "annotate", "bitops", "edgeR", "limma", "S4Vectors",
              "VariantAnnotation", "rmarkdown", "plotly", "scran",
              "data.table", "kableExtra", "htmlwidgets",
              "webshot", "clusterProfiler", "dupRadar", "pheatmap",
              "taxize", "SingleCellExperiment", "SummarizedExperiment",
              "scater", "DropletUtils", "shiny", "heatmaply", "readxl",
              "readr", "dplyr", "shinycssloaders", "shinyjs", "slingshot",
              "Rmagic", "reticulate", "viridis", "Seurat", "tidyverse",
              "httr", "jsonlite", "xml2", "writexl", "zip")
packages <- setdiff(packages, rownames(installed.packages()))
BiocManager::install(packages)

remotes::install_github("velocyto-team/velocyto.R")
```

## Dependencies of external software
* bwa, bowtie, bowtie2, STAR, picard, sambamba, samtools, igvtools
* lsof


## Installation of the development version of `ezRun` from github
```R
remotes::install_github("uzh/ezRun")
```

## Development of `ezRun` package at FGCZ environment
Always at the conda environment `ezRun` during the development. The conda environment contains the necessary external tools/software.


## Coding style
Do follow the guidelines in [CodingStyle.md]

