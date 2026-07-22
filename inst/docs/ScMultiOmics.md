# ScMultiOmics Extension Gotchas

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
- **exploreSC link in the rendered HTML — cwd-stripping heuristic, NOT
  `output$getColumn("Report")`.** Every render before 2026-04-30 emitted a
  broken link: the old setup chunk did
  `rp <- file.path(proj, basename(rp))`, which throws away every intermediate
  segment between the project ID and the report dir. So
  `p31662/Analyses_Paul/scmultiomics_smoke_paul/<report_dir>` was written
  out as `p31662/<report_dir>` — exploreSC then 404'd silently
  (`Server module skipped: data load returned NULL`). The issue affected
  every run, not just ad-hoc smoke drivers; earlier wording in this README
  blaming smoke drivers was wrong.
- **Fix in `ScMultiOmics.Rmd` setup chunk (current):**
  1. Two render params take precedence — `params$explore_url` (literal URL)
     and `params$gstore_data_path` (gstore-relative dir without filename).
  2. Fallback now derives from `getwd()` by keeping the FULL remainder
     after `/pXXXXX/`, not just the basename:
     `regmatches(getwd(), regexec("/(p\\d+)/(.+)$", getwd()))[[1]]` →
     `paste0(proj, "/", subpath, "/scMultiData.qs2")`.
- **When patching an already-delivered HTML on gstore**, copy the file out,
  `sed` the broken URL to the working one, then `g-req copynow -f` it back.
  Cheaper than re-rendering the whole report just to fix the link.

### `ScMultiOmicsApp` output dataset must propagate `refFeatureFile` for `ScSeuratCombineApp`
- `ScSeuratCombineApp.rb` gates on
  `@required_columns = ['Name', 'Species', 'refBuild', 'refFeatureFile', 'Static Report']`.
  The original `ScMultiOmicsApp#next_dataset` emitted only
  `Name, Species, refBuild, Static Report [Link], Report [File], ScMultiOmics [Link], SC Seurat [Link]`
  — no `refFeatureFile`, so SUSHI refused to launch the combine job at the
  column-validation step before R was ever invoked. The
  `scData.qs2 → scMultiData.qs2` symlink fix only addressed half of the
  compatibility problem.
- Fix landed 2026-04-30 in `ScMultiOmicsApp.rb`: declared
  `@params['refFeatureFile'] = 'genes.gtf'`, inherited it in
  `set_default_parameters` from `@dataset[0]['refFeatureFile']` when present,
  and emitted it in `next_dataset`. `ScSeuratFilterClustersApp` and
  `ScSeuratLabelClustersApp` only require
  `Name, Species, refBuild, SC Seurat`, so they were already compatible via
  the symlink — only Combine was blocked.
- General lesson: when adding a new SUSHI app whose output is meant to feed
  any downstream app, scrape every downstream `.rb`'s `@required_columns`
  set and make sure `next_dataset` covers their union. Don't trust that
  inheriting `Factor` / `B-Fabric` tags fills the gap — those are user
  metadata, not pipeline-required columns.

### Input dataset shape: CellRanger column-name mismatch (RESOLVED 2026-05-02)
**Superseded by the ScSeurat-only input change below — kept for history.**
- `ScMultiOmicsApp.@required_columns` literally required `CountMatrix`, but
  **CellRanger Multi outputs do not emit a `CountMatrix` column** — they
  emit only `ResultDir [File,Link]` pointing at the per-sample dir
  (`p31662/o41361_CellRangerMulti_*/pUM07_LT1_dark`). CellRanger ARC,
  CellRangerCount, CellBender, and BDRhapsody DO emit `CountMatrix [Link]`
  natively, so the column gate only blocked CR Multi.
- Consequence at the time: the single-click Multi → ScMultiOmics flow
  never worked end-to-end; it always needed manual TSV editing to inject a
  `CountMatrix [Link]` column.
- Resolution: rather than patching the runtime to fall back to `ResultDir`,
  we made `ScMultiOmics` consume **ScSeurat output** instead of raw
  CellRanger. ScSeurat now propagates `CountMatrix [Link]` and
  `ResultDir [File]` in its `next_dataset`, so all 10x flavors (Multi /
  ARC / Count) end up with the canonical column shape downstream. See
  the "ScSeurat-only input" section below.

### Architectural decision (shipped 2026-05-02): ScSeurat-only input
- `ScMultiOmicsApp` now consumes **only** ScSeurat output (or BD Rhapsody
  directly), never raw CellRanger. Reason: ScSeurat has already filtered dead
  cells / mito high / ribo high / SoupX-corrected ambient RNA. Layering
  ADT/VDJ/ATAC on the raw CellRanger filtered matrix means clonotype counts,
  ADT QC, and WNN run on cells you'd discard anyway. QC'd cells are the only
  sensible base.
- What landed in this patch:
  1. **`ScSeuratApp.rb#next_dataset`** now propagates `CountMatrix [Link]` and
     `ResultDir [File]` (mirrors the prior `refFeatureFile` fix). Both are
     required inputs to ScSeurat (line 18 of `ScSeuratApp.rb`), so the values
     are guaranteed present in `@dataset` and only needed forwarding.
  2. **`ScMultiOmicsApp.@required_columns`** switched from `['Name', 'Species',
     'refBuild', 'CountMatrix']` to `['Name', 'Species', 'refBuild',
     'SC Seurat']`. `CountMatrix` is now propagated for free via ScSeurat.
  3. **`findScDataPath()`** in `app-ScMultiOmics.R` simplified to a one-liner
     gated on the `SC Seurat` column. The `SC Cluster Report` / `Report` /
     `Static Report` fallback walking is gone (dead branch — the SUSHI gate
     guarantees `SC Seurat`).
  4. The CR Multi column-name mismatch (above) is now irrelevant: ScSeurat
     normalises every 10x input (Multi, ARC, Count) to a canonical output
     shape with `SC Seurat`, `CountMatrix`, and `ResultDir` all populated.
- BD Rhapsody decision: BD bypasses ScSeurat at the runtime level
  (`SCDataOrigin = BDRhapsody` branch, loads `*_Seurat.rds` directly). BD's
  pipeline does its own cell filtering upstream, so re-QCing in ScSeurat
  would mostly be cosmetic. The SUSHI `@required_columns` gate now enforces
  `SC Seurat`, so BD datasets need either a hand-built `dataset.tsv` (with
  a `SC Seurat` placeholder column) or — cleaner long-term — a thin
  `BDRhapsodyToScSeurat` SUSHI app that wraps the BD RDS and writes a
  SUSHI-shaped `scData.qs2`. Everything would then funnel through
  `SC Seurat`.
- Tradeoff: this closes the door on directly running ScMultiOmics from a
  CR dataset (you must always go via ScSeurat). For the 99% case that's
  correct — the QC'd object is what the multi-omics report should be built on.

### Manual deployment via `scp + cp` vs git: untracked-file footgun
- Deploying a SUSHI app file via `scp ScMultiOmicsApp.rb fgcz-h-082:/tmp/`
  + `ssh trxcopy@fgcz-h-082 'cp /tmp/... /srv/sushi/production/master/lib/'`
  works for an immediate hot-fix, but it leaves the file as **untracked**
  in the production checkout. A subsequent `git pull` aborts with
  "The following untracked working tree files would be overwritten by
  merge: lib/ScMultiOmicsApp.rb" — even when the untracked content is
  byte-identical to what's incoming.
- Resolution: `rm` the untracked file before pulling. Confirm equivalence
  first with `diff <(git show origin/master:master/lib/X.rb) lib/X.rb`.
- Production checkout also accumulates phantom `M` flags on files whose
  local edits eventually got upstreamed via PRs; `git diff origin/master`
  shows zero, but `git status` still flags them. `git checkout -- <file>`
  clears the stat cache (no content change since they're already
  identical to origin/master).

### `update_sushi_prod` script: false success on pull failure
- The bash function in `~/.bash_aliases` on fgcz-c-053 runs:
  `ssh trxcopy@fgcz-h-082 'cd /srv/sushi/production && git pull'`
  followed by `ssh ... 'systemctl restart apache2'` and an unconditional
  `echo "SUSHI production updated successfully!"`. **None of the steps
  short-circuit on failure** — a failed `git pull` still triggers the
  Apache restart and the success line, leaving production running on
  whatever was on disk before the pull (which can be a stale or partial
  state if a cleanup was in progress).
- Patched 2026-04-30 to gate each step with `if ! ssh ...; then return 1; fi`
  and emit a stderr error line. The next `update_sushi_prod` after a pull
  failure now exits with code 1 *before* restarting Apache.
- Same pattern likely affects neighboring `restart_sushi` and any
  `update_sushi_course` / `_demo` helpers — apply the same guard.

### DEPLOY.md smoke checklist references fixtures, not SUSHI datasets
- The `[ ]` checklist lines in `ScMultiOmicsApp_DEPLOY.md` ("re-run on
  p31662 P8_AT_PBMCs", "p39132 IC104_Plate_6217_1_ARC", "p39179
  394581_1-CD34run1") are aspirational. P8_AT_PBMCs is a downsampled
  smoke fixture under `/srv/GT/analysis/p31662/Analyses_Paul/scmultiomics_smoke/`,
  not a SUSHI-registered dataset.tsv row. The IC104 / CD34 entries DO
  correspond to real SUSHI datasets but they were never actually run
  end-to-end through the production SUSHI submission path.
- For real production smoke testing, pick datasets that already have
  `CountMatrix [Link]` populated (CellRanger ARC, BDRhapsody) and a
  paired ScSeurat output if going through the planned ScSeurat-only
  flow. Examples that work today:
  - BD Rhapsody (no ScSeurat needed):
    `p39179/o39458_BDRhapsodySA_2025-09-03--15-54-00`
  - CR ARC (multiome, needs ScSeurat first):
    `p40259/o40270_CellRangerARCCount_2026-04-13--14-01-46` (5 mouse
    samples, no ScSeurat output exists yet — must run ScSeurat as
    prerequisite)
  - CR Multi (CITE-seq + TCR, needs both ScSeurat AND a TSV edit until
    the column-name mismatch is fixed):
    `p31662/o41361_CellRangerMulti_2026-03-31--13-10-24` paired with
    `p31662/o41361_ScSeurat_2026-04-28--16-11-43`

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

### `obj[[col, drop=TRUE]]` returns CHARACTER even when `Idents` is a factor with custom levels
- Subtle trap that broke the cell-type DotPlot diagonal even after the fix
  above was supposedly applied. The previous `adt_top_per_group()` did
  `groups <- obj[[group_col, drop=TRUE]]` and then
  `level_order <- if (is.factor(groups)) levels(groups) else sort(unique(...))`.
  For Azimuth's `predicted.celltype.l2` the meta.data column may be
  *character* even though `Idents(obj)` (assigned from the same column)
  is a factor with non-alphabetical levels — so the `is.factor(groups)`
  branch never fired and `level_order` came out alphabetical, NOT in the
  Idents/DotPlot x-axis order. The reorder block then produced an `avg`
  matrix that was alphabetical, and `best_group` indices pointed at the
  wrong columns of the rendered DotPlot.
- Symptom: the DotPlot looks "approximately diagonal" because the bright
  cell type for each marker is *correct* (z-score is order-invariant),
  but the FEATURE rows are placed at the wrong y-axis positions —
  rolling a clean diagonal into a scrambled one. Spearman(x_pos, y_pos)
  drops from ±1 to ~±0.4.
- Fix: read the canonical group order directly from `Idents()`, not from
  the metadata column:
  `groups <- droplevels(Idents(obj)); level_order <- levels(groups)`.
  After the fix, `colnames(avg)` matches `levels(Idents(obj))` exactly
  (modulo `_` → `-`), and Spearman correlation jumps to ±0.998.
- Diagnostic recipe (from the debug session): build the DotPlot, extract
  `p$data`, find each feature's max-`avg.exp.scaled` cell type, and
  compare its x-axis index in `levels(p$data$id)` to the feature's y-axis
  index in `levels(p$data$features.plot)`. Perfect diagonal ⇒ Spearman ±1.

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

### Smoke test drivers must source `scData.qs2` from a real ScSeurat run
- Hand-rolled drivers that build `scData.qs2` from a count matrix +
  `NormalizeData → PCA → FindClusters → RunUMAP` do NOT run Azimuth, so
  the resulting object has no `predicted.celltype.l2` (or any other
  `pickCellTypeColumn()` candidate). Every "by cell type" panel in the
  rendered report falls through to the "_no annotation column_"
  placeholder.
- For comparable output, point `SC Seurat [Link]` at the gstore-resident
  ScSeurat scData.qs2 (`/srv/gstore/projects/pXXXXX/.../*_SCReport/scData.qs2`)
  rather than rebuilding RNA-only on the fly. The CR v10+ Hbt_257 driver
  is fine for "does the pipeline run" smoke testing but generates reports
  without cell-type panels.

### CLR is fast; ADTnorm is the production default
- `ScMultiOmicsApp.rb` advertises `adtNorm = ['ADTnorm', 'CLR']` (first item
  is the SUSHI default). Smoke drivers in this repo set `adtNorm = "CLR"`
  for speed — fine for "does it run" tests, but the rendered report will
  read `normalization: CLR` in the QC tab, and the ADT UMAP separation
  will look noticeably worse than a SUSHI-driven run with ADTnorm.
- ADTnorm is single-threaded and slow on >50k cells (10-15 min on a typical
  CITE-seq panel); CLR is sub-second. Use CLR only for iteration loops
  where you don't care about the final ADT UMAP quality.

### `g-req copynow` races with concurrent SUSHI loops on the same gstore path
- A SUSHI app actively iterating on `pXXXXX/Analyses_Paul/<dir>/<report>`
  triggers `g-req` trash + copy cycles every few minutes. If a side
  re-render copies into the same path, it survives until the next SUSHI
  copy lands (or hits "Destination path already exists" and gets dropped).
- For ad-hoc smoke renders during a live SUSHI loop, copy to a parallel
  path (`pXXXXX/Analyses_Paul/<dir>_smoke_paul/<report>`) to avoid the
  race. `gstore-list` makes the contention visible.

### ADT modality correlation: literal `match()` misses ~80% of TotalSeq panels
- Original code at [`_scMultiOmics_adt.Rmd`](inst/templates/_scMultiOmics_adt.Rmd)
  did `match(sub("\\.\\d+$", "", adtFeats), rnaFeats)` — only resolves ADT
  names that are *literally* gene symbols. On the BioLegend TotalSeq-C
  Human Universal Cocktail v1 (143 antibodies, p31662), this matched
  26/143 = 18%; the remaining 117 markers showed no correlation panel
  even though most have an obvious gene counterpart.
- Markers that escape the literal rule include the clone-suffixed names
  (`CD4(RPA-T4)` → `CD4`), the parenthesized-payload form
  (`CD326-(EPCAM)` → `EPCAM`, `CD278-(ICOS)` → `ICOS`), markers whose
  gene differs from the antigen (`CD56` → `NCAM1`, `CD64` → `FCGR1A`,
  `CD16` → `FCGR3A`, `CD8` → `CD8A`, `CD45*` → `PTPRC`), Ig isotypes
  (`IgD` → `IGHD`, `Ig-light-chain-kappa` → `IGKC`), and TCR variable
  chains (`TCR-Valpha7.2` → `TRAV1-2`).
- Fix: `mapADTtoRNA()` in `multiOmicsUtils.R` driven by
  `inst/extdata/adt_to_gene_lookup.tsv` (~300 hand-verified rows
  covering TotalSeq Universal Cocktail v1/v2 + TotalSeq-A/B/C variants).
  Falls back to literal equality, then regex (`CDxxx-(GENE)` payload,
  `STEM(CLONE)` stem). Records `matchSource` so the QC tab can show how
  each marker resolved. Coverage on p31662 panel: 136/136 (100%) of
  non-isotype markers, all via lookup.
- The lookup-first / strip-second order is critical: stripping
  `\\.\\d+$` first turns `TCR-Valpha7.2` into `TCR-Valpha7` and
  guarantees a miss. Try the raw name against the lookup before
  applying the make.unique-collision strip.

### Existing protein↔gene resources (for future maintenance)
- **`AbNames` R package** by Helen Lindsay (UZH):
  https://github.com/HelenLindsay/AbNames — same problem, more
  comprehensive scope. Ships `data(totalseq)` (977 BioLegend antibodies
  with `Antigen`, `Clone`, `Cat_Number`, `HGNC_SYMBOL`, `ENSEMBL_ID`,
  `TotalSeq_Cat`), `data(citeseq)` (3,212 antibodies seen in published
  studies), `data(gene_aliases)` (HGNC alias table), plus matching
  helpers `searchTotalseq`, `getCommonName`, `renameADT`, `formatIg`,
  `formatTCR`, `replaceGreekSyms`. Aiming for Bioconductor; not on CRAN
  yet.
- We don't take a hard dep on AbNames because it has gaps for several
  textbook markers (CD16, CD56, CD8 all have `HGNC_SYMBOL = NA` in
  `totalseq` as of 2026-04 — the same markers we get right via the
  curated TSV). Once AbNames lands on Bioconductor and fills those
  gaps, switch to it as the primary data source and reduce our TSV to
  FGCZ-specific overrides only.
- BioLegend itself publishes per-panel **feature-reference CSVs** at
  `biolegend.com/Files/Images/BioLegend/totalseq/<PANEL>_<CAT#>_Antibody_reference_UMI_counting*.csv`,
  but the `name` column is the antibody label, NOT the gene symbol — so
  these CSVs are useful for `cellranger --feature-ref` but not directly
  for ADT↔RNA correlation lookup.
