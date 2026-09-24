---
app: EzAppCellRangerMulti
defined_in: R/app-cellRangerMulti.R (ezMethodCellRangerMulti, prepareFastqData, buildMultiConfigFile)
report_templates: []
called_by:
  - SUSHI master/lib/CellRangerMultiApp.rb
last_checked:
  date: 2026-09-24
  ezrun_commit: dea3ffc9
  sushi_commit: 78f7007b
checked_against:
  - ezRun R/app-cellRangerMulti.R (ezMethodCellRangerMulti, prepareFastqData, getSampleMultiplexFiles, getCellRangerMultiData, buildMultiConfigFile, EzAppCellRangerMulti)
  - ezRun R/app-cellRanger.R (getCellRangerGEXReference, getCellRangerVDJReference, link10xFastqPaths, subsample, cellRangerAnnotatableRef)
  - ezRun R/util.R (tarExtract)
  - ezRun R/01classes.R (write_methods, which passes config.csv to the methods writer)
  - SUSHI master/lib/CellRangerMultiApp.rb
  - installed CellRanger 10.1.0 mro/rna/sc_multi_cs.mro (which internal stages exist)
  - FGCZ reference layout under /srv/GT/reference (index folder naming)
  - cellranger-fgcz skill, section on local cell annotation
verified_with_run: >-
  No. Written from the code, the CellRanger pipeline definition and the
  cellranger-fgcz skill. Not yet checked against the logs and outputs of a real run.
---

# EzAppCellRangerMulti

## Purpose

Runs 10x Genomics `cellranger multi` on one 10x library, or one pool of libraries,
per job. ezRun's part is preparation: it stages the FASTQs, makes sure a CellRanger
reference exists, and writes the multi configuration CSV. `cellranger multi` then
does all of the analysis itself, inside one vendor pipeline. Use it for any 10x
single-cell design that needs more than plain gene expression: immune repertoire,
antibody capture, sample multiplexing, or Flex (fixed RNA). For a plain gene
expression library, `EzAppCellRanger` (`cellranger count`) is the simpler route.

## Accepted inputs

Library types are chosen in `TenXLibrary` (several at once):

| Library type | What it is | Input columns (tar / FASTQ) |
|---|---|---|
| `GEX` | 3' or 5' gene expression | `RawDataDir` / `Read1` + `Read2` |
| `VDJ-T` | T-cell receptor | `VdjTDataDir` / `VdjT Read1` + `VdjT Read2` |
| `VDJ-B` | B-cell receptor | `VdjBDataDir` / `VdjB Read1` + `VdjB Read2` |
| `FeatureBarcoding` | Antibody capture (CITE-seq / ADT) | `FeatureDataDir` / `Feature Read1` + `Feature Read2` |
| `Multiplexing` | Several samples in one library | depends on `MultiplexingType`, below |
| `fixedRNA` | 10x Flex, probe-based fixed RNA profiling | GEX columns; GEX is added automatically |

Multiplexing types (`MultiplexingType`):

- `ocm`: on-chip multiplexing. The sample information is in the GEX reads; no extra library.
- `antibody`: hashtag antibodies (HTO), read from the Feature library.
- `cmo`: 10x CellPlex lipid barcodes, read from the `MultiDataDir` library.
- With `fixedRNA`, multiplexing uses Flex probe barcodes, and `MultiplexingType` does not apply.

Other input facts:

- Each library comes either as a tar archive (legacy) or as FASTQ paths.
- `RawDataDir` may list several tars, comma-separated. They are pools from one chip
  and are combined into one CellRanger run.
- Multiplexed runs also need `<prefix>_Sample2Barcode.csv` files in the order's
  `o<order>_metaData/` folder on gStore. These map each sample to its barcodes. The
  file prefix must match the start of the dataset sample `Name`.
- The input dataset usually has to be curated by hand for anything beyond plain GEX.
- One job per input row. A row is one library or pool, and it may contain many
  biological samples when multiplexed.

## Stages

In data-flow order. "Trace" says where, if anywhere, the run record shows the stage.
The run record for this app has three parts: the job script, the job logs, and the
`config.csv` in the result folder, which is CellRanger's copy of the configuration
it actually received.

1. **Input staging.** In FASTQ mode, the FASTQ files are symlinked into a folder
   layout CellRanger expects. In tar mode, the tars are extracted. When a tar's
   inner folder name differs from the dataset `Name`, the folder and files are
   renamed. Trace: the extraction and linking are not logged; renames appear as
   `EXECUTED CMD: mv ...` / `rename ...`; multi-pool runs log a
   `Multi-pool mode: combining ...` message.
2. **Read subsampling (conditional).** Tar mode only, and only when `nReads` is
   greater than 0. `nReads` is not on the SUSHI form; it can only be set through
   `specialOptions`. Reads are subsampled with seqtk. Trace: one
   `EXECUTED CMD: seqtk sample ...` line per file, with the seed and read number on
   the command line.
3. **Gene expression reference (for GEX / fixedRNA).** Two cases:
   - *Shared reference (default).* A CellRanger index kept next to the selected
     annotation, one per combination of transcript types. It is reused when it
     exists. On first use it is built: the annotation is first reduced to the
     transcript types chosen in `transcriptTypes`, then `cellranger mkref` builds
     the index.
   - *Per-run custom reference.* Used when `controlSeqs`, `secondRef` or
     `extendThreePrime` is set. A new reference is built for this run: the genome
     plus the extra sequences and their annotation, optionally with gene 3' ends
     extended, then `cellranger mkref`. It is deleted after the run.

   Trace: when an index is built, an `EXECUTED CMD: cellranger mkref ...` line. When
   an existing shared index is reused there is no log line; the reference path in
   `config.csv` is the only record. See "Recorded settings" for what that path
   tells you. The reduction to transcript types is never logged.
4. **VDJ reference (for VDJ-T / VDJ-B).** A shared CellRanger VDJ index built from
   the same genome and the full annotation. Reused when it exists, otherwise built
   with `cellranger mkvdjref`. Trace: the `EXECUTED CMD: cellranger mkvdjref` line
   when built; otherwise only the `[vdj]` reference path in `config.csv`.
5. **Flex probe set preparation (for fixedRNA).** The selected probe set is reduced
   to probes whose gene ID and gene name both occur in the reference annotation.
   Probes from `customProbesFile` are added first, with gene and probe IDs given a
   `Gene_` prefix where missing. The result is written to a `_filtered.csv` file.
   Trace: the `probe-set` line in `config.csv` points to the filtered file. The
   filtering itself is not logged.
6. **Chemistry.** Written to the configuration only when set by the user, except for
   Flex. For Flex it is always written: taken from `chemistry` if set; otherwise
   inferred, for Flex v2 probe sets from the read length of the first R1 file, and
   for Flex v1 from whether the run is multiplexed. Trace: the chemistry line in
   `config.csv`; the Flex v2 inference also logs a `Flex v2 auto-detect: ...` message.
   Without a chemistry line, CellRanger detects the chemistry itself; the detected
   value is in CellRanger's own summary outputs, not in the configuration.
7. **Feature and multiplexing references.**
   - Antibody capture uses the uploaded `FeatureBarcodeFile`.
   - For HTO or CMO multiplexing, the chosen barcode set (`MultiplexBarcodeSet`) is
     reduced to the barcodes named in the Sample2Barcode files.
   - When antibody capture and HTO multiplexing are combined, both go into one
     feature reference, because CellRanger accepts only one.

   The sample-to-barcode assignment from the Sample2Barcode file goes into the
   `[samples]` section. Trace: the `[feature]`, `cmo-set` and `[samples]` entries in
   `config.csv`. The barcode filtering is not logged.
8. **cellranger multi.** CellRanger runs with the configuration file. Trace:
   `EXECUTED CMD: cellranger multi --id=... --localmem=... --localcores=... --csv=...`,
   plus anything from `cmdOptions`. The CellRanger version is the
   `module load Aligner/CellRanger/<version>` line in the job script. `config.csv`
   is authoritative for what CellRanger was given.
9. **Inside cellranger multi (vendor pipeline).** Depending on the library types:
   - read alignment and UMI counting per gene (GEX / Flex)
   - cell calling
   - antibody counting (Feature Barcoding)
   - V(D)J contig assembly, annotation and clonotype grouping (VDJ)
   - assignment of cells to their sample of origin (Multiplexing)
   - secondary analysis per sample: dimensionality reduction, clustering and
     differential expression between clusters

   ezRun offers no way to switch off the secondary analysis, so it always runs. Trace:
   that these steps ran is shown by CellRanger's outputs. How they were done is not
   recorded, except for what `config.csv` sets.
10. **Clean-up.** BAM files are deleted unless `keepBam` is on. When both `keepBam`
    and `secondRef` are set, each BAM is converted to CRAM with samtools and the BAM
    removed. Trace: `EXECUTED CMD: find ... -delete` or `EXECUTED CMD: samtools view ... -C`.
11. **Bookkeeping, not analysis.** An `expanded_dataset.tsv` with one row per
    demultiplexed sample is written to a shared folder and copied into the result
    folder with `g-req`.

## Stages that leave no trace

- **CellRanger's internal methods (stage 9).** Alignment, cell calling, clustering,
  dimensionality reduction, differential expression, V(D)J assembly, clonotype
  grouping and cell-to-sample assignment all happen inside `cellranger multi`. None
  of their methods or settings appear in the log. `config.csv` shows only the inputs
  it lists. Each step that applies to the chosen library types was performed and
  should be reported. Its settings are `[not recorded]`. They must not be filled in
  from CellRanger's documentation.
- **Transcript-type reduction of the annotation (stage 3).** Not logged. The chosen
  types are visible only in the `transcriptTypes` parameter and in the reference
  folder name in `config.csv`.
- **Flex probe-set filtering (stage 5)** and **barcode-set filtering (stage 7).**
  Not logged; only the resulting file paths appear in `config.csv`.
- **Input staging (stage 1).** Tar extraction and FASTQ linking are not logged.

## Steps that did not happen

These are listed so they are left out entirely, not described and not denied:

- **CellRanger's built-in cell-type annotation.** CellRanger runs its local annotation
  model only when the reference's declared genome name is one it recognises as
  human. FGCZ references declare their index folder name instead, so annotation is
  skipped, with no error. `EzAppCellRanger` works around this; this app does not. No
  `cell_types` output is produced.
- **Intronic counting for Flex.** The job script can show `includeIntrons` as true
  for a Flex run, but ezRun never passes it to CellRanger for `fixedRNA`. Whether
  introns were counted is decided by `config.csv`, not by the job script.

## Recorded settings: which ones affect results

On the `cellranger multi` command line:

- Housekeeping, no effect on the result: `--id` (name of the output folder),
  `--localmem` and `--localcores` (memory and threads), `--csv` (where the
  configuration file is).
- May affect the result: anything added through `cmdOptions`.

In `config.csv`, every entry affects the result except `create-bam`, which only
decides whether an alignment file is written. File paths in `config.csv` are
evidence, not text to reproduce. What each one tells you:

- **Gene expression reference path.** FGCZ references follow the layout
  `<organism>/<source>/<genome build>/Annotation/Release_<release>-<date>/Genes/<index>`.
  The path gives the organism, the annotation source (for example Ensembl or
  GENCODE), the genome build and the annotation release. The `-<date>` after the
  release number is when FGCZ set the reference up; it is not part of the annotation
  release. For a shared index, the index folder name
  `genes_10XGEX_SC_<types>_Index` lists the transcript types the annotation was
  reduced to. A per-run custom reference is a folder named `10X_customised_Ref`; its
  genome and annotation are the selected ones plus the extra sequences.
- **VDJ reference path.** Same layout; the index folder is `genes_10XVDJ_Index`.
- **Probe set** (Flex). The file name gives the 10x probe set and its version, with
  a `_filtered` suffix added by stage 5.
- **Feature reference.** An antibody panel supplied by the user. For HTO runs it
  also holds the hashtag barcodes in use.
- **`[samples]` section.** Which samples were multiplexed in the library, and with
  which barcodes.
- **`include-introns`, `expect-cells`, `chemistry`.** Report them as given, when present.

In the job script, the scheduler resources (cores, memory, scratch) and module lines
are housekeeping, except as the source of tool versions.

## Requested versus applied

The job script records what was requested. `config.csv` records what CellRanger
received. Where they differ, `config.csv` wins:

- `includeIntrons` is dropped for `fixedRNA`.
- `expectedCells` is written only when set. Otherwise CellRanger estimates the cell
  number itself.
- `chemistry` for plain GEX is written only when set. Otherwise CellRanger detects it.
- `keepBam`: the SUSHI form's value overrides the ezRun default. The BAM is always
  produced (`create-bam,true`) and then kept or deleted according to `keepBam`.

## Loaded but never run

The job script loads modules for seqtk and samtools. seqtk runs only in stage 2 and
samtools only in the CRAM case of stage 10. A loaded module is not evidence that the
tool was used.

## Parameters that change what runs

Every parameter is written to the job script as a `param[['<name>']] = '<value>'` line.

- `TenXLibrary`: which library types are processed, and so which CellRanger analyses run.
- `MultiplexingType`: HTO, CMO or OCM sample multiplexing.
- `MultiplexBarcodeSet`: barcode set for HTO / CMO.
- `probesetFile`, `customProbesFile`: Flex probe set, and added custom probes.
- `FeatureBarcodeFile`: antibody feature reference.
- `refBuild`, `refFeatureFile`, `transcriptTypes`: which reference is used or built.
- `controlSeqs`, `secondRef`, `extendThreePrime`: switch to a per-run custom
  reference. `extendThreePrime` is not on the SUSHI form; set it through
  `specialOptions`.
- `includeIntrons`, `chemistry`, `expectedCells`: passed into the configuration as
  described above.
- `nReads`: subsampling (stage 2). Tar mode only, and only through `specialOptions`.
- `keepBam`: keeps BAM files; with `secondRef` they are stored as CRAM.
- `cmdOptions`: appended as-is to the `cellranger multi` command.
- `CellRangerVersion`: which CellRanger is loaded. It also changes the output folder
  layout.

## Outputs

- `ResultDir [File,Link]`: `<result>/<Name>/`, CellRanger's `outs` folder renamed to
  the sample name. Contains `config.csv` and, per demultiplexed sample,
  `per_sample_outs/<sample>/` with its count matrices and web summary.
- `Report [Link]`: the per-sample web summary (CellRanger 9 and earlier) or the
  combined QC report (CellRanger 10 and later).
- One further row per demultiplexed sample, only for multiplexed runs. The rows are
  built from the Sample2Barcode file and carry `CountMatrix`,
  `UnfilteredCountMatrix`, `ResultDir` and `Report` for that sample. They are the
  usual input to downstream single-cell apps: `CellBender` takes them directly, and
  `ScSeurat` does too when the input dataset carries a `Condition` column.
- The main output row has no `CountMatrix` column. A run without multiplexing
  therefore produces no row that `ScSeurat` or `CellBender` can take directly.
- The layout of `per_sample_outs` changed between CellRanger 9 and 10. The column
  paths follow the selected `CellRangerVersion`.

## Limitations and gotchas

- The dataset `Name` must equal the GEX FASTQ file prefix. Otherwise CellRanger's
  preflight fails within a minute.
- Dropdown parameters left at their `select` placeholder break the run when it is
  submitted outside the web form. See the `cellranger-fgcz` skill.
- Sample2Barcode files are matched to samples by name prefix. A wrong or missing
  file stops multiplexed runs.
- An existing shared reference is trusted as complete. A run waits while another job
  is building the same reference, and gives up after a fixed timeout.
- For how to run, re-run and debug this app on FGCZ infrastructure, see the
  `cellranger-fgcz` and `sample2barcode-generation` skills.
