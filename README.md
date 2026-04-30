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

### Cluster colours must be homogenised across panels
- Each `pals::polychrome()` call returns the same vector, but `DimPlot`,
  `VlnPlot`, scRepertoire bars, and `clonalOverlap` all map colours to factor
  levels in their own way. If one tab uses `seurat_clusters` (factor levels
  `"0", "1", "10", ..., "9"`) and another uses `wsnn_res.0.5` (different
  level set), cluster 0 ends up two different colours in the same report.
- Build a single `cluster_pal` (named-by-level palette) once in the parent
  `ScMultiOmics.Rmd` setup, force the cluster column on the object to a
  numeric-sorted factor (`as.character(sort(as.integer(levels)))`), and pass
  `cols = cluster_pal` to every plot that fills/colors by cluster
  (`DimPlot`, `VlnPlot` weights/QC, `cloneOccupyFull`).
- `VlnPlot()` defaults to ggplot2 hue palette when `cols=` is omitted —
  silently breaks the shared palette. Always pass `cols = cluster_pal`.

### scRepertoire `clonalOverlap`: `order.by` doesn't reach the ggplot factor levels
- `clonalOverlap(..., order.by = cluster_levels)` reorders the underlying
  overlap matrix correctly, but the heatmap layer factorises axis names
  alphabetically anyway. Result: x/y axes render in character order
  (`0, 1, 10, 11, …, 2, 3, …`) regardless of `order.by`.
- Fix: append
  `+ scale_x_discrete(limits = cluster_levels) + scale_y_discrete(limits = cluster_levels)`
  to the returned ggplot. Same trick may apply to other scRepertoire
  heatmap-style outputs.

### Numeric cluster IDs need explicit `order.by` for scRepertoire bar plots
- scRepertoire defaults to character sorting for `group.by` factor levels,
  giving `0, 1, 10, 11, 2, 3, …` on the x-axis of `clonalQuant`,
  `clonalHomeostasis`, `clonalProportion`, `clonalDiversity`,
  `clonalLength` etc.
- Pass `order.by = cluster_levels` (numeric-sorted character vector) to
  every scRepertoire plotting function.

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

### ADTnorm fails on individual markers — never silently fall back
- ADTnorm's per-marker peak-fitting (`fda::smooth.morph` → `lnsrch_morph`)
  errors with `Initial slope not negative.` (or `subscript out of bounds`)
  on markers with near-flat or zero-heavy distributions. A single failing
  marker aborts the whole batch call — on the p31662 P8_AT_PBMCs panel,
  52/143 markers fail.
- ADTnorm has **no built-in multi-threading** — upstream `R/ADTnorm.R`
  runs a `for(adt_marker_index in adt_marker_index_list)` serial loop
  with a `TODO: pending parallel running setup` comment.
- Fix in `processADT()`: per-marker `parallel::mclapply` over `param$cores`
  workers. Each marker runs in its own worker with its own `tryCatch`;
  failed markers fall back to CLR for those rows only, the rest of the
  panel gets ADTnorm-aligned. Both speeds the call up and isolates failures.
- A previous version did
  `if (normMethod == "ADTnorm" && requireNamespace("ADTnorm"))` and
  silently dropped to CLR when the package was missing — that hid the
  ADTnorm-not-installed-on-trxcopy gap. Now it `stop()`s loudly.
- Surface the per-marker outcome in the rendered HTML
  (`@misc$adtnorm` → ADT QC chunk shows `n_ok / n_markers` plus the
  full failed-marker list). Silent fallback is bad UX.

### ADTnorm + CLR mixed-scale spliced matrix
- After per-marker mclapply, ADTnorm rows are on the
  `negPeak / valley / posPeak` 0-1-2 landmark scale; CLR-fallback rows
  are on the log-ratio scale (typically -2 to 5). Splicing them into one
  assay is fine for *within-marker* comparisons (DotPlots z-score per
  row), but cross-marker magnitude comparisons are meaningless. Make
  this clear in the report — the QC chunk does.
- Splicing must respect Seurat's `_` → `-` rewrite on feature names
  (`CreateAssay5Object` warns "Feature names cannot have underscores");
  match good-marker names with `gsub("_", "-", ...)` before indexing
  into the CLR baseline matrix or you hit `subscript out of bounds`.

### ADTnorm does NOT subtract isotype background
- ADTnorm aligns each marker's distribution onto the
  `negPeak / valley / posPeak` 0-1-2 landmark scale; the negative peak
  implicitly handles background, but isotype counts are **not used as a
  denominator or subtractor**. ADTnorm's `ADTnorm()` signature has no
  isotype parameter.
- For true isotype-aware normalization use `dsb`
  (Decontamination of Single-cell Background, `niaid/dsb`). Not currently
  wired into `processADT()` — would need a `normMethod = "dsb"` branch
  with an empty-droplets matrix as input.
- Our pipeline filters isotypes from RANKING (DotPlot top-3 helper drops
  `(?i)isotype` rows) but does no background subtraction.

### ADTnorm install on the SUSHI compute lib
- ADTnorm is not on CRAN or Bioconductor; install via
  `remotes::install_github("yezhengSTAT/ADTnorm")`.
- Its dep tree (~152 packages including `flowStats`, `flowCore`,
  `flowWorkspace`, `EMDomics`, `fda`) needs Fortran toolchain pieces
  trxcopy doesn't have — `urca` / `forecast` / `flowWorkspace` fail with
  `-lgfortran` errors when installed cold.
- Pragmatic workaround: install ADTnorm and its full dep chain into the
  user's lib first (where the Fortran symlinks are set up), then `rsync`
  the resolved binaries to `/home/trxcopy/R/x86_64-pc-linux-gnu-library/4.5/`.
  Same architecture + R version → binary-compatible.

### CellRanger Multi hashtag features pollute the Antibody Capture slot
- CellRanger Multi runs that demultiplex via TotalSeq hashtags emit each
  hashtag as a feature in the H5's `Antibody Capture` group alongside
  real ADT antibodies. Treating hashtags as an ADT modality (per-cell
  PCA, WNN) is meaningless — each cell expresses one hashtag and is
  near-zero for the others, producing a degenerate ADT PCA whose NaN
  scores break `FindMultiModalNeighbors`.
- `detectModalities()` now parses `config.csv`
  (`[samples] hashtag_ids` column), subtracts those IDs from the H5
  Antibody Capture features, and reports `hasADT = FALSE` if all that's
  left are hashtags. Without this, hashtag-only Multi runs incorrectly
  trigger the ADT branch and crash WNN.

### NFS cache lag after `R CMD INSTALL` on shared user lib
- Right after a fresh install of `ezRun` into `/home/<user>/R/...`,
  SBATCH jobs that immediately `library(ezRun)` on a compute node
  sometimes see the OLD package and error with
  `object 'EzAppScMultiOmics' not found` — the compute node's NFS cache
  hasn't refreshed.
- Fix: wait ~30 s after `R CMD INSTALL` returns, then submit. If the job
  hits the stale lib, just resubmit — the next attempt picks up the
  fresh files.

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

### BD Rhapsody `*_Seurat.rds` ships ATAC as `peaks`, not `ATAC`
- BD multiome RDS objects expose chromatin under `Assays(obj)` named `peaks`,
  `gene_activity`, and (sometimes) `chromvar`. The ATAC child Rmd keys off
  `obj[["ATAC"]]` and `obj@reductions$lsi`, so `loadBDRhapsody()` renames in
  place: `obj[["ATAC"]] <- obj[["peaks"]]; obj[["peaks"]] <- NULL`, same for
  `gene_activity` -> `GeneActivity`.
- BD also doesn't compute LSI / atac.umap; `loadBDRhapsody()` runs
  `Signac::RunTFIDF()` -> `RunSVD()` -> `RunUMAP(reduction = "lsi", reduction.name = "atac.umap")`
  when ATAC is present and Signac is installed.

### BD VDJ data lives in 3 places — use the AIRR file
- `*_Seurat.rds` carries dominant-chain wide metadata (`TCR_Alpha_Gamma_*`,
  `BCR_Heavy_*` …). Convertible but lossy.
- `*_VDJ_perCell.csv` mirrors the Seurat metadata.
- **`*_VDJ_Dominant_Contigs_AIRR.tsv`** is the canonical, scRepertoire-native
  source. `loadBDContigs(resultDir)` reads it via `read.delim()` then routes
  through `scRepertoire::loadContigs(input = list(df), format = "BD")` for a
  long contig table with TRA/TRB/TRG/TRD + IGH/IGK/IGL chains in one file.
- Caveat: `loadContigs()` directory-mode looks for files literally named
  `Contigs_AIRR.tsv` and ignores BD's `<sample>_VDJ_*` prefix. Always pass
  the file as a pre-loaded data frame.

### `Seurat::AverageExpression()` returns alphabetical columns regardless of factor levels
- `DotPlot()`'s x-axis follows `Idents()` factor-level order (often
  non-alphabetical from Azimuth / scType). When the top-features-per-group
  helper indexes columns of `AverageExpression()` (which IS alphabetical),
  `best_group=1` doesn't line up with the leftmost x-axis column and the
  diagonal collapses.
- Fix in `adt_top_per_group()`: reorder the avg matrix to `levels(Idents())`
  before computing z-scores. Drop low-count idents from the obj subset
  (`droplevels(Idents(obj))`) so they don't surface as empty x-axis columns.

### CellRanger v10+ Multi layout — no code changes needed
- CR Multi v10 keeps the `per_sample_outs/<sample>/sample_filtered_feature_bc_matrix.h5`
  + sibling `vdj_t/` / `vdj_b/` layout from earlier versions. Existing
  `detectModalities()` (looks for `*filtered_feature_bc_matrix\.h5$`) and
  `findVDJContigCsv()` work unchanged.
- The H5 carries both `Gene Expression` and `Antibody Capture` feature_types,
  so `h5HasAntibodyCapture()` lights up `mod$hasADT` automatically.

### Feature naming across modalities — prefix the assay name
- Protein names (CD3, CD4, …) routinely collide with gene symbols.
  `FetchData()`, `FeaturePlot()`, `DoHeatmap()`, and `FindMarkers()` only
  disambiguate when `DefaultAssay()` is set or the feature is prefixed with
  the assay key.
- Convention: prefix ADT features `ADT_<name>` (key `adt_`), ATAC peaks
  inherit `chr-start-end` (no collision risk, key `atac_`), GeneActivity
  features `GA_<gene>` (key `ga_`). Matches Seurat's
  `Key(obj[["ADT"]]) == "adt_"` so users can plot
  `FeaturePlot(obj, features = c("CD3E", "adt_CD3"))` without flipping
  `DefaultAssay`.
- Apply on assay creation, not after. Saved `scMultiData.qs2` objects from
  earlier runs predate this rule and need a one-shot rename script if
  re-rendered through a feature-name-sensitive helper.

### `Read10X_h5(unique.features = TRUE)` adds `.1` to ADT names that collide with RNA
- Make.unique runs across the *concatenated* feature list before the matrix
  is split by feature_type. So an ADT named "CD14" becomes "CD14.1" in the
  Antibody Capture sub-matrix when the same symbol exists in Gene
  Expression. The trailing suffix leaks into DotPlot / FeaturePlot labels.
- The Seurat Key system (`adt_`, `rna_`) only disambiguates inside
  `FetchData()` lookups; it never modifies raw `rownames(GetAssayData(...))`,
  which is what plot labels read.
- Fix in `readADTCounts()`: after extracting the ADT submatrix, strip
  `\.[0-9]+$` only from rows whose stripped name appears in the RNA feature
  set (those are the rows make.unique actually touched). Real protein names
  don't end in `.<digit>`, so the rule is safe.
- Saved `scMultiData.qs2` objects from before the fix retain the `.1`
  names; rebuild the assay in place with the helper script
  (`fix_adt_names_<sample>.R` template under `script/`) before re-rendering.

### CellRanger Multi hashtags pollute the ADT assay — parse `config.csv` to drop them
- CR Multi runs that demultiplex via TotalSeq hashtags expose those features
  in the per-sample H5 alongside (or instead of) real CITE-seq antibodies.
  Treating them as ADT in WNN is meaningless — each cell expresses one
  hashtag and is near-zero for the others, the resulting ADT PCA is
  degenerate, and `FindMultiModalNeighbors` produces a graph with NaN
  edges that `RunUMAP` rejects with `neighbor graph contains missing data`.
- Authoritative source: the `config.csv` next to (or 1-3 levels above) the
  count matrix has a `[samples]` section with `hashtag_ids` listing exactly
  which feature IDs are hashtags.
- Fix in `findCellRangerMultiHashtags()`: walk up six levels looking for
  `config.csv`, parse the `[samples]` section, return the IDs. Used by both
  `detectModalities()` (sets `hasADT = FALSE` when *all* antibodies are
  hashtags) and `readADTCounts()` (drops hashtag rows from the matrix).

### `processADT()` defensive NaN handling for tiny ADT panels
- With small ADT panels (<= 6 features, common when CR Multi mixes hashtags
  + a few real antibodies), CLR(margin = 2) on a cell with all-zero ADT
  counts returns -Inf, and `ScaleData()` divides by per-feature SD which
  is 0 for any constant feature, returning NaN in `scale.data`. Either
  path injects NaN into the ADT PCA scores and the WNN graph fails.
- Fix in `processADT()`: drop cells with `colSums(adtCounts) == 0` before
  the assay is created; replace any remaining non-finite values in the
  `data` and `scale.data` layers with 0 before `RunPCA()`.

### scRepertoire 2.7+ proportion-based cloneSize labels miss the count-based palette
- scRepertoire 2.6 emitted `cloneSize` levels like `"Hyperexpanded (100 < X <= 500)"`;
  2.7+ defaults to proportion thresholds → `"Hyperexpanded (0.1 < X <= 1)"`,
  `"Medium (0.001 < X <= 0.01)"`, etc., and renames the smallest non-zero
  bin from `Single` to `Rare`. Both schemes share the same conceptual
  ordering but the literal strings don't overlap.
- A palette keyed only by the count-based labels makes every cell fall
  through to `na.value = "grey85"` in `scale_color_manual()`, so the entire
  cloneSize UMAP renders grey.
- Fix: `cloneSizePalette()` now contains entries for *both* schemes
  (Single / Rare share the same colour). `cloneSizeMatchPalette(levels)`
  returns the slice that actually appears in the data so
  `scale_color_manual(limits = ...)` is driven by the real labels — same
  fix applied in `cloneOccupyFull()` and `cloneDimPlot()`.

### `cloneDimPlot()` layering — explicit "No clonotype" background + larger foreground dots
- After converting NA cloneSize to the explicit `"No clonotype"` level
  (required so exploreSC doesn't crash on NA factors), `order(!is.na(...))`
  no longer pushes those cells to the back — they're now non-NA and mix
  in front of real clonotypes, making the plot look uniformly grey.
- Fix: split the geom into two `geom_point()` layers — `bg_point_size = 0.3`
  alpha-0.5 for "No clonotype" cells (drawn first → behind), and
  `point_size = 1.4` for real clonotypes (drawn last → on top), with
  the foreground sub-sorted so Hyperexpanded ends up at the very top.
