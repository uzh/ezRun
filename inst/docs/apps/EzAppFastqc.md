---
app: EzAppFastqc
defined_in: R/app-fastQC.R (ezMethodFastQC)
report_templates:
  - inst/templates/FastQC.Rmd (only with showNativeReports)
  - inst/templates/FastQC_overview.Rmd (only with showNativeReports)
called_by:
  - SUSHI master/lib/FastqcApp.rb
last_checked:
  date: 2026-09-24
  ezrun_commit: dea3ffc9
  sushi_commit: 78f7007b
checked_against:
  - ezRun R/app-fastQC.R (ezMethodFastQC, EzAppFastqc)
  - ezRun R/fastqIO.R (countReadsInFastq, ezMethodSubsampleFastq, subsampleFastqFile)
  - ezRun R/app-trim.R (ezMethodFastpTrim)
  - ezRun R/system.R (ezSystem)
  - ezRun R/00defaults.R (ezParam, parseOptions)
  - ezRun inst/extdata/EZ_PARAM_DEFAULTS.txt
  - ezRun inst/templates/FastQC.Rmd
  - SUSHI master/lib/FastqcApp.rb
  - installed FastQC 0.12.1 Configuration/limits.txt (which modules are switched off)
verified_with_run: >-
  One production run (2026-09-24, paired-end, both AI summaries on): job script,
  logs and methods.md. Not verified against a run: stage 1 (fastp) and stage 4
  (subsampling), which are described from the code only.
---

# EzAppFastqc

## Purpose

Read-level quality control of FASTQ files from bulk (non-10x) sequencing. FastQC is
run on every read file, and the results are aggregated into one MultiQC report. Two
optional steps add LLM-written summaries to that report. The app only produces
reports: it does not alter or filter the reads that later apps use. For 10x
single-cell reads, use `EzAppFastqc_10x` instead.

## Accepted inputs

- FASTQ, single-end or paired-end (`Read1`, optional `Read2`). SUSHI sets `paired`
  from whether the dataset has a `Read2` column. When `paired` is off, `Read2` is
  ignored.
- One sample may list several FASTQ files in a single cell, comma-separated. They are
  concatenated before anything else runs.
- Any organism and any assay. No reference genome is used.
- Unmapped BAM: conversion code exists but is disabled. Treat FASTQ as the only
  supported input.
- All samples of a dataset are processed in one job.

## Stages

In data-flow order. "Trace" says where, if anywhere, the run record shows the stage.

1. **Read preprocessing with fastp (conditional).** Runs only when `max_len1` is
   greater than 0. The SUSHI form does not offer `max_len1`; it can only be set
   through `specialOptions`. When it runs, it is not just length capping: it is the
   full ezRun fastp preprocessing using ezRun's global trimming settings (adapter
   trimming and the other fastp steps), and every later stage sees the trimmed reads.
   Trace: an `EXECUTED CMD: fastp ...` line with every flag it received. fastp's own
   report goes to a per-sample `_preprocessing.log` that is not kept in the result.
2. **Concatenation of multi-file samples (conditional).** Only when a sample lists
   several files. The files are joined into one FASTQ per read. Trace: partial. The
   `touch` commands are logged as `EXECUTED CMD`; the `cat` commands are not.
3. **Read counting (conditional).** Only when the input dataset has no `Read Count`
   column. Reads are counted with ShortRead. The count decides stage 4 and feeds the
   native report. Trace: none.
4. **Subsampling (conditional).** Only when the total read count across the dataset
   exceeds a threshold fixed in the code. Each sample is then subsampled with
   ShortRead to a fixed number of reads, using a random seed hard-coded in the
   function, and FastQC runs on the subsample. Trace: indirect only. The input file
   names on the FastQC command line end in `-subsample_R1.fastq.gz` (and `_R2`).
   The threshold, the target read number and the seed are code constants and do
   not appear in the run record.
5. **FastQC.** One FastQC run over all read files. FastQC assesses each file in a set
   of modules: per-base and per-sequence quality, per-tile quality, per-base sequence
   content, per-sequence GC content, per-base N content, sequence length
   distribution, sequence duplication, overrepresented sequences, and adapter
   content. Its k-mer module is switched off in the FGCZ FastQC installation (see
   Requested versus applied). Which modules produced results is visible in the
   record as MultiQC section names, in the `[AI-Sections]` log lines when stage 9
   ran. Trace: the full command line in the job's stdout log as
   `EXECUTED CMD: fastqc ...`. FastQC runs quietly and its own output files are
   deleted, so its version appears only in the job script's `module add` line
   (`QC/FastQC/<version>`).
6. **Native report rendering (conditional).** Only when `showNativeReports` is on.
   Renders one overview page per FastQC plot type across all samples, plus a main
   report with read counts per sample and a PASS/WARN/FAIL table per FastQC module.
   When off, the per-file FastQC report folders are deleted. Trace: the parameter
   value in the job script, and the output files.
7. **MultiQC aggregation.** MultiQC collects all FastQC results into
   `multi_FastQC/`. Trace: no command line. The call is made with logging switched
   off, on purpose, to keep the internal LLM endpoint out of the logs. What the log
   does show: MultiQC's own banner with its version, its per-module lines (for example
   how many FastQC reports it found), and ezRun's
   `[MultiQC] STARTED ... (generate_ai_summary=..., per_section_ai_summaries=..., model=...)`
   and `[MultiQC] FINISHED` lines.
8. **Whole-report AI summary (conditional, report annotation only).** Only when
   `generate_ai_summary` is on. MultiQC's built-in AI summary feature, pointed at the
   FGCZ-internal LLM, writes a summary at the top of `multiqc_report.html`. The
   summary is text to help read the report: it changes no data and produces no result.
   It is not part of the analysis and does not belong in a Methods description.
   Trace: the `[MultiQC] STARTED` line shows `generate_ai_summary=TRUE`; the prompts
   are saved in `multiqc_data/llms-full.txt`.
9. **Per-section AI summaries (conditional, report annotation only).** Only when
   `per_section_ai_summaries` is on. For each MultiQC section, ezRun sends that
   section's data table to the same LLM and inserts one or two bullet points into the
   report. Like stage 8, this is report annotation, not analysis, and does not belong
   in a Methods description. Trace: one `[AI-Sections] <i>/<N> <section name>` line
   per section; the prompts are saved in `multiqc_data/sushi_section_prompts.txt`.
   The section names in these lines are still useful evidence: they show which
   FastQC modules produced results (see stage 5).

## Stages that leave no trace

- **Read counting (stage 3).** No log line at all.
- **Subsampling (stage 4).** No log line. The only clue is the `-subsample_` file
  names in the FastQC command. The threshold, read number and seed are not in the
  record.
- **MultiQC settings (stage 7).** That MultiQC ran, and its version, are recorded. How
  it was configured is not, because its command line is deliberately kept out of the
  log. Report MultiQC with its version and mark its settings `[not recorded]`. Do not
  describe it as run "with default parameters": nothing in the record says so.
- **Concatenation (stage 2).** The joining commands themselves are not logged.
- **fastp's own report (stage 1).** The command line is logged, but fastp's report is
  not kept.

## Recorded settings: which ones affect results

The FastQC command line mixes settings that shape the result with housekeeping.

Affect the result:

- `-a <file>`: the list of adapter sequences FastQC searches for. It is an in-house
  FGCZ list. Describe it as such; the file path is not part of a description.
- Anything added through `cmdOptions`, which appears at the end of the command.
- Which files FastQC was given: original, concatenated (stage 2), trimmed (stage 1)
  or subsampled (stage 4). The file names on the command line show this.
- `paired`: whether `Read2` files were assessed.

Housekeeping, no effect on the result:

- `--extract` (unpack the report), `-o` and `--dir` (output and temporary folders),
  `-q` (quiet), `-t` (threads), and the redirection of FastQC's output to files.
- The scheduler resources in the job script (cores, memory, scratch).
- The module lines, except as the source of tool versions.

## Requested versus applied

- `--kmers` is on the FastQC command line, but FastQC's k-mer module is switched off
  in the FGCZ installation's configuration, so the setting has no effect. No k-mer
  analysis was performed; do not mention a k-mer size or k-mer content.
- `perLibrary` is declared in the app defaults, but its code path is disabled; it
  has no effect.

## Loaded but never run

The job script loads modules for Picard, samtools and fastp. This app never calls
Picard or samtools. fastp runs only in stage 1. A loaded module is not evidence
that the tool was used.

## Parameters that change what runs

Every parameter is written to the job script as a `param[['<name>']] = '<value>'` line.

- `paired`: whether `Read2` files are processed.
- `max_len1`: switches on stage 1 (full fastp preprocessing). Reachable only through
  `specialOptions`.
- `specialOptions`: free-form `key=value` pairs; any ezRun parameter can be set here,
  including ones the SUSHI form does not show.
- `cmdOptions`: appended as-is to the FastQC command; visible on its command line.
- `showNativeReports`: switches on stage 6 and keeps the per-file FastQC reports.
- `generate_ai_summary`: switches on stage 8.
- `per_section_ai_summaries`: switches on stage 9.

## Outputs

- `FastQC [File]`: folder of per-file FastQC results. The native HTML reports are
  kept only with `showNativeReports`.
- `MultiQC [File]` and `MultiQC Report [Link]`: `multi_FastQC/multiqc_report.html`,
  plus `multiqc_data/` holding the per-section tables and, when stages 8 or 9 ran,
  the AI prompt files.
- `FastQC Report [Link]`: the native main report; only with `showNativeReports`.
- The output dataset has one row, named `FastQC`, and carries the `Order Id` column
  over from the input. No other app takes it as input: FastQC is a final QC step.

## Limitations and gotchas

- Reads are assessed, never modified for later use. Even when stage 1 runs, the
  trimmed reads are not exported.
- After subsampling (stage 4), every QC statistic describes the subsample, not the
  full read set.
- Two input files with the same base name produce the same report name and stop
  the run.
- For 10x reads use `EzAppFastqc_10x`.
- Log lines such as `Unable to remove path: .../.multiqc_tmp/...` and
  `Couldn't remove tmp dir` are harmless cleanup noise, not an analysis step.
- The AI summaries (stages 8 and 9) are generated text inside the report, meant to
  help interpretation. They run as part of the job, but they are not part of the
  analysis: they change no data and produce no result.
