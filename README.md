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
Do follow the guidelines in [CodingStyle.md](CodingStyle.md)

## ScMultiOmics extension — gotchas

Captured while iterating on the `scmultiomics` branch (single-cell multi-omics
report layered on top of an annotated ScSeurat output). For deployment and
submission notes see
[`sushi/master/lib/ScMultiOmicsApp_DEPLOY.md`](https://github.com/uzh/sushi/blob/scmultiomics/master/lib/ScMultiOmicsApp_DEPLOY.md).

### scRepertoire `combineExpression` — proportion-based default cloneSize bins
- Default `clone.size = c(Rare = 1e-4, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded = 1)` is **unreachable** for any group with < 1,000 TCR cells. With 824 TCR cells, every singleton lands at `1/824 ≈ 0.0012` → already in `Medium`, so `Single`/`Small`/`Rare` are mathematically empty.
- Defaults are unchanged between scRepertoire 2.4 and 2.7.3 — upgrading does not help.
- Fix: pass count-based bins explicitly. `processVDJ()` now uses
  `proportion = FALSE, cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500)`.
- `cloneSizePalette()` keys must match the bin formula scRepertoire stamps onto factor levels (`"Single (0 < X <= 1)"`, …).

### scRepertoire 2.7+ `vizGenes()` API change
- The `scale = TRUE` argument was removed; new API uses
  `summary.fun = c("percent", "proportion", "count")`.
- The `safe_print` wrapper swallows the resulting "unused argument" error, so
  V/J panels go silently blank — symptom is empty plots, not a render failure.

### Seurat v5 `AverageExpression()` returns a data.frame
- `AverageExpression(obj, assays = "ADT", layer = "data")$ADT` is a data.frame
  in Seurat 5, which breaks `t()` / `scale()` with `argument is not a matrix`.
- Wrap the return in `as.matrix()` before any matrix algebra (see
  `adt_top_per_group()` in `_scMultiOmics_adt.Rmd`).

### `Seurat::DotPlot()` + `coord_flip()` puts `features[1]` at the BOTTOM
- Counter-intuitive but verified empirically: with `coord_flip()`, the FIRST
  element of the `features =` vector lands at the **bottom** of the y-axis,
  the LAST at the top.
- To put cluster 0's markers at the top of the plot, **reverse** the marker
  list before passing it to `DotPlot`. Both the ADT helpers
  (`adt_top_per_group()`) and the WNN per-modality DEG subtab apply `rev()`.
- Do **not** try to fix this with `scale_y_discrete(limits = ...)` — that
  refers to the *data* y-axis (the cluster ident, post-flip), not the visual
  y-axis, and silently filters every point out.

### `DimPlot(group.by, order = TRUE)` does not reliably foreground sparse highlights
- With ~70 highlighted cells against ~9,000 grey, the red points get visually
  drowned even with `order = TRUE` and a small `pt.size`.
- Build the plot manually: `geom_point(grey)` first, then `geom_point(red, size = 1.6)`
  on top. See `_scMultiOmics_vdj.Rmd` "Top 3 clones on UMAP".

### `pickCellTypeColumn()` finds nothing on standard PBMC ScSeurat output
- ezRun's `pickCellTypeColumn()` scans ~25 candidate columns
  (`predicted.celltype.l2`, `CyteType`, `azimuth_pan_human`, …). Default
  ScSeurat output has none, so the "by cell type" panels render the
  "No cell-type annotation column on the object" placeholder.
- For PBMC samples, run `Azimuth::RunAzimuth(obj, reference = "pbmcref")`
  upstream of the report. The `predicted.celltype.l2` column is then picked
  up automatically.

### ATAC annotation needs `refBuild`
- `processATAC()` resolves the EnsDb annotation from `refBuild`. If it's NULL,
  no annotation is attached → `TSSEnrichment` and `GeneActivity` are skipped
  silently with the "genome annotation could not be resolved" placeholder.
- `ezMethodScMultiOmics` falls back to `input$getColumn("refBuild")` when
  `param$refBuild` is empty (smoke-test scripts that bypass SUSHI usually
  populate the input dataset, not the param list).

### `ezInteractiveTableRmd` loses its dependencies in `results='asis'` chunks
- `print(htmltools::tagList(ezInteractiveTableRmd(mk, ...)))` inside an
  `asis` chunk emits the widget HTML but does **not** register its DT JS/CSS
  dependencies. Symptom: the rendered HTML has only ~2 occurrences of
  `DataTable` (just the title strings) instead of the ~20+ a working
  interactive table produces, and the table is invisible / stripped.
- Fix: put the dotplot + interactive table in a **child Rmd** (e.g.
  `_scMultiOmics_wnn_deg_subtab.Rmd`) and call it via
  `knitr::knit_child()` from a small `asis` loop in the parent. The child's
  chunks are non-asis, so the htmlwidget's deps load normally.

### `ScSeuratCombine` compatibility — `scData.qs2` symlink + dataset column
- `ScSeuratCombineApp.rb` reads `input$getColumn("SC Seurat")` and expects
  each path to end at a `scData.qs2`. `ezMethodScMultiOmics` saves the
  multi-modal object as `scMultiData.qs2`, so the combine job can't find
  anything by default.
- Fix landed on this branch:
  1. `app-ScMultiOmics.R` creates `scData.qs2 → scMultiData.qs2` (relative
     symlink) in the report dir after `makeRmdReport` returns.
  2. `ScMultiOmicsApp.rb#next_dataset` adds
     `'SC Seurat [Link]' => File.join(report_file, "scData.qs2")` so the
     output dataset is consumable by `ScSeuratCombine` /
     `ScSeuratCombinedLabelClusters` without user intervention.

### exploreSC silent crash on missing input
- If the `?data=` URL points to a file that does not exist at any
  `urlDataRoot` (`/srv/gstore/projects` or
  `/srv/GT/analysis/course_sushi/public/gstore/projects`),
  `execute_import()` triggers `show_error_modal("File Not Found")` +
  `stopApp()` and returns `NULL`. The user sees a stalled spinner /
  "Connection reset", with an empty error message in the outer trace.
- Fixed in shinyproxy_apps `613b407`:
  `parse_url_and_setup_dirs()` now logs `File not found. Tried: <both candidate paths>`,
  `execute_import()` has an explicit `is.na(dataDir)` guard, and the outer
  `tryCatch` skips `shiny.silent.error` instead of dressing it up as a real
  error.
- Diagnose via session stderr at
  `/srv/shinyproxy/logs/exploreSC_<proxy-id>_<date>_stderr.log`.

### `fgcz-shiny.uzh.ch` is on fgcz-c-051
- Distinct from `shiny.fgcz.uzh.ch` (private, port 8080) and
  `shiny-public.fgcz.uzh.ch` (public, port 8081). The
  `fgcz-shiny.uzh.ch` reverse proxy resolves to fgcz-c-051 — same server,
  different port and config (`/srv/shinyproxy/private.yml`).

### Tiny smoke fixtures may silently skip panels
- The smoke fixtures in
  `/srv/GT/analysis/p31662/Analyses_Paul/scmultiomics_smoke/` are deliberately
  downsampled (~2000 cells × 5000 genes). Tabs that gate on
  `nrow(obj[[assay]]) >= 2L` may silently skip on these. When something looks
  "missing" in a smoke render, first check the fixture dimensions before
  changing the template.

### `g-req copynow` is queued (eventually consistent)
- The CLI prints `successfully copied` immediately, but the file lands in
  `/srv/gstore/projects/...` ~15-30 s later (g-req drops the request into a
  queue processed by `gstore-list`). Don't `ls` for it right after the copy
  command — wait or poll. For overwriting an existing file, use
  `g-req copynow -f` (force).

### SUSHI submission gotchas
- **`sushi_fabric --class ScMultiOmicsApp` fails** with `class cannot be loaded`
  because the app's `initialize()` calls `ref_selector()` → `Rails.cache.fetch()`
  which is `nil` outside the full Rails environment. Workaround: skip
  `sushi_fabric` and submit a direct SBATCH that mirrors the script SUSHI
  would have generated, then call `EzAppScMultiOmics$new()$run(input, output, param)`.
- **Compute jobs run as `trxcopy`**, so `library(ezRun)` looks in
  `~trxcopy/R/x86_64-pc-linux-gnu-library/4.5/`. Installing into your own
  `~/R/...` is not enough — package the source, scp to fgcz-h-082, install
  with `R CMD INSTALL --library=/home/trxcopy/R/x86_64-pc-linux-gnu-library/4.5`
  as `trxcopy`.
- **`g-req -w copy <dir>` fails when the destination already exists**
  ("Destination path already exists"). Re-rendering on top of a previous
  output requires `g-req -w remove <dest>` first, then `g-req -w copy`. For
  individual files, `g-req copynow -f` (force) is fine.
- **exploreSC link in the rendered HTML** is built from
  `output$getColumn("Report")`. SUSHI populates this with the gstore-relative
  output path automatically; smoke-test drivers must set it themselves
  (otherwise the URL collapses to `pXXXXX/<basename>/scMultiData.qs2`, which
  404s in exploreSC).

