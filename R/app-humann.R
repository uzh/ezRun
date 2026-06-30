###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

# HUMAnN 4 alpha functional profiling. Auto-detects upstream profile vintage,
# runs HUMAnN in full (tier-1 + tier-2) or translated-only mode accordingly,
# CPM-renormalizes, optionally regroups to KEGG KO, and renders a per-sample
# HTML report. The conda env (`gi_qiime2-amplicon-2026.4`) is sourced by the
# SUSHI Ruby `commands` block before this script runs, so humann / python
# tools are already on PATH.
#
# HUMAnN 4 alpha caveats (per biobakery forum #8523):
#   - chocophlan pangenomes are pinned to mpa_vOct22_CHOCOPhlAnSGB_202403.
#     Newer MetaPhlAn DBs (vJan25 ...) break SGB matching, so they fall
#     back to translated-only mode.
#   - MetaPhlAn >= 4.2 has an API change HUMAnN 4 can't drive itself; we
#     always supply --taxonomic-profile externally (HUMAnN never auto-
#     invokes MetaPhlAn in this app).
#   - The UniRef DB is EC-filtered only; gene-family quantitation is
#     restricted to enzymes. MetaCyc pathways work as designed; KEGG KO
#     regroup covers the EC-mapped subset only.
#   - HUMAnN 4 renamed UNMAPPED -> READS_UNMAPPED in gene families, but
#     humann_renorm_table's --special filter doesn't yet recognize the
#     new name; we strip it manually before renorm.

ezMethodHUMAnN <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  sampleName <- input$getNames()

  ## ---- 1. Trim + (if paired) concat reads --------------------------------
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  read1 <- trimmedInput$getColumn("Read1")
  if (isTRUE(as.logical(param$paired))) {
    read2 <- trimmedInput$getColumn("Read2")
    concatFq <- paste0(sampleName, "_concat.fastq.gz")
    ezSystem(sprintf("cat %s %s > %s",
                     shQuote(read1), shQuote(read2), shQuote(concatFq)))
  } else {
    concatFq <- read1
  }

  ## ---- 2. Detect upstream profile column (exactly one) -------------------
  hasMpa     <- "MetaPhlAnProfile" %in% input$colNames
  hasBracken <- "BrackenReport"    %in% input$colNames
  hasKraken  <- "KrakenReport"     %in% input$colNames
  nSources   <- sum(c(hasMpa, hasBracken, hasKraken))
  if (nSources == 0) {
    stop("Input dataset has none of: MetaPhlAnProfile, BrackenReport, KrakenReport.")
  }
  if (nSources > 1) {
    stop("Input dataset carries multiple profile columns; ambiguous. ",
         "Pick a single-source dataset and re-run.")
  }

  ## ---- 3. Prepare HUMAnN-compatible profile + determine mode -------------
  humannProfile <- paste0(sampleName, "_humann_profile.tsv")
  mode <- "translated"
  detectedDb <- NA_character_
  profileSrc <- NA_character_

  if (hasMpa) {
    profileSrc <- "MetaPhlAnProfile"
    srcFile <- input$getColumn("MetaPhlAnProfile")
    detectedDb <- readLines(srcFile, n = 1, warn = FALSE)
    if (grepl("vOct22_CHOCOPhlAnSGB_202403", detectedDb, fixed = TRUE)) {
      mode <- "full"
      file.copy(srcFile, humannProfile, overwrite = TRUE)
    } else {
      ## Header-rewrite so HUMAnN accepts the file in translated-only mode
      ezSystem(sprintf("sed '1s/.*/#mpa_vOct22_CHOCOPhlAnSGB_202403/' %s > %s",
                       shQuote(srcFile), shQuote(humannProfile)))
    }
  } else {
    profileSrc <- if (hasBracken) "BrackenReport" else "KrakenReport"
    srcFile <- input$getColumn(profileSrc)
    detectedDb <- sprintf("(converted from %s)", profileSrc)
    converter <- system.file("python/bracken_to_mpa.py", package = "ezRun")
    if (!nzchar(converter) || !file.exists(converter)) {
      stop("Missing bracken_to_mpa.py in ezRun (inst/python/). Reinstall the package.")
    }
    ezSystem(sprintf("python3 %s %s -o %s",
                     shQuote(converter), shQuote(srcFile), shQuote(humannProfile)))
  }

  ## ---- 4. forceTranslatedOnly override -----------------------------------
  forced <- isTRUE(as.logical(param$forceTranslatedOnly))
  if (forced && mode == "full") {
    message("forceTranslatedOnly=TRUE — downgrading from full to translated-only.")
    mode <- "translated"
  }

  ## ---- 5. Run HUMAnN -----------------------------------------------------
  dbRoot <- param$humannDbRoot
  if (is.null(dbRoot) || !nzchar(dbRoot)) dbRoot <- "/srv/GT/databases/humann4_vOct22"
  protDb  <- file.path(dbRoot, "uniref")
  utilDb  <- file.path(dbRoot, "utility_mapping")
  chocoDb <- file.path(dbRoot, "chocophlan")
  for (d in c(protDb, utilDb)) {
    if (!dir.exists(d)) stop(sprintf("Required HUMAnN DB dir not found: %s", d))
  }
  if (mode == "full" && !dir.exists(chocoDb)) {
    stop(sprintf("ChocoPhlAn DB not found at %s (required for full mode).", chocoDb))
  }

  outDir <- paste0("humann_out_", sampleName)
  dir.create(outDir, showWarnings = FALSE)
  logFile <- paste0(sampleName, ".humann.log")

  humannCmd <- paste(
    "humann",
    "--input",             shQuote(concatFq),
    "--output",            shQuote(outDir),
    "--taxonomic-profile", shQuote(humannProfile),
    "--protein-database",  shQuote(protDb),
    "--utility-database",  shQuote(utilDb),
    "--threads",           ezThreads(),
    "--memory-use",        "maximum",
    "--remove-temp-output",
    "--resume",
    "--verbose",
    "--o-log",             shQuote(logFile)
  )
  if (mode == "full") {
    humannCmd <- paste(humannCmd,
                       "--nucleotide-database", shQuote(chocoDb))
    if (!isTRUE(as.logical(param$keepStratifiedOutput))) {
      humannCmd <- paste(humannCmd, "--remove-stratified-output")
    }
  } else {
    humannCmd <- paste(humannCmd,
                       "--bypass-nucleotide-search",
                       "--remove-stratified-output")
  }
  ezSystem(humannCmd)

  ## ---- 6. Per-sample post-processing -------------------------------------
  rawGf <- Sys.glob(file.path(outDir, "*_2_genefamilies.tsv"))
  rawRx <- Sys.glob(file.path(outDir, "*_3_reactions.tsv"))
  rawPa <- Sys.glob(file.path(outDir, "*_4_pathabundance.tsv"))
  if (length(rawGf) == 0 || length(rawRx) == 0 || length(rawPa) == 0) {
    stop("HUMAnN did not produce the expected three output files in ", outDir)
  }

  gfCpm <- sprintf("%s_genefamilies_cpm.tsv", sampleName)
  rxCpm <- sprintf("%s_reactions_cpm.tsv",    sampleName)
  paOut <- sprintf("%s_pathabundance.tsv",    sampleName)

  ## Strip READS_UNMAPPED (HUMAnN 4 alpha) before renorm — humann_renorm_table
  ## doesn't recognize the new name and would otherwise inflate the denominator.
  gfFiltered <- sprintf("%s_genefamilies.norm_in.tsv", sampleName)
  ezSystem(sprintf("grep -v '^READS_UNMAPPED\\b' %s > %s",
                   shQuote(rawGf[1]), shQuote(gfFiltered)))
  ezSystem(sprintf("humann_renorm_table -i %s -u cpm -s n -o %s",
                   shQuote(gfFiltered), shQuote(gfCpm)))
  unlink(gfFiltered)
  ezSystem(sprintf("humann_renorm_table -i %s -u cpm -s n -o %s",
                   shQuote(rawRx[1]), shQuote(rxCpm)))
  file.copy(rawPa[1], paOut, overwrite = TRUE)

  ## Optional KEGG KO regroup (non-fatal)
  koOut <- sprintf("%s_genefamilies_ko_cpm.tsv", sampleName)
  if (isTRUE(as.logical(param$regroupKEGGKO))) {
    koRaw <- sprintf("%s_genefamilies_ko_cpm.raw.tsv", sampleName)
    koOk <- tryCatch({
      ezSystem(sprintf("humann_regroup_table -i %s -g uniref90_ko -o %s",
                       shQuote(gfCpm), shQuote(koRaw)))
      tryCatch(
        ezSystem(sprintf("humann_rename_table -i %s -n kegg-orthology -o %s",
                         shQuote(koRaw), shQuote(koOut))),
        error = function(e) {
          message("humann_rename_table failed; keeping raw KO IDs.")
          file.rename(koRaw, koOut)
        }
      )
      unlink(koRaw)
      TRUE
    }, error = function(e) {
      message("KEGG KO regroup failed: ", conditionMessage(e))
      unlink(c(koRaw, koOut))
      FALSE
    })
  }

  ## Write a tiny mode marker so the SUSHI Ruby wrapper can surface the
  ## detected mode in dataset.tsv (Mode [Characteristic] column).
  modeFile <- sprintf("%s.mode.txt", sampleName)
  writeLines(mode, modeFile)

  ## ---- 7. Render the per-sample HTML report ------------------------------
  rmdSrc <- system.file("templates/HUMAnNReport.Rmd", package = "ezRun")
  if (!nzchar(rmdSrc) || !file.exists(rmdSrc)) {
    stop("Missing HUMAnNReport.Rmd in ezRun (inst/templates/). Reinstall the package.")
  }
  file.copy(rmdSrc, "HUMAnNReport.Rmd", overwrite = TRUE)
  rmarkdown::render(
    input = "HUMAnNReport.Rmd",
    envir = new.env(),
    output_dir = ".",
    output_file = htmlFile,
    quiet = TRUE,
    params = list(
      sample        = sampleName,
      mode          = mode,
      profileSource = profileSrc,
      detectedDb    = detectedDb,
      forced        = forced,
      humannLog     = file.path(outDir, logFile),
      gfCpmFile     = gfCpm,
      rxCpmFile     = rxCpm,
      paFile        = paOut,
      koFile        = if (file.exists(koOut)) koOut else NA_character_
    )
  )

  ## Clean up the concat fastq (large, regenerable)
  if (isTRUE(as.logical(param$paired)) && file.exists(concatFq)) {
    unlink(concatFq)
  }

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodHUMAnN()
##' @templateVar htmlArg )
##' @description Use this reference class to run HUMAnN 4 alpha functional profiling.
EzAppHUMAnN <-
  setRefClass(
    "EzAppHUMAnN",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodHUMAnN
        name <<- "EzAppHUMAnN"
        appDefaults <<- rbind(
          paired = ezFrame(
            Type = "logical", DefaultValue = "FALSE",
            Description = "Paired-end reads."),
          forceTranslatedOnly = ezFrame(
            Type = "logical", DefaultValue = "FALSE",
            Description = "Force translated-only mode even when a vOct22 MetaPhlAn profile is present."),
          humannDbRoot = ezFrame(
            Type = "character", DefaultValue = "/srv/GT/databases/humann4_vOct22",
            Description = "HUMAnN 4 reference DB root (chocophlan/, uniref/, utility_mapping/)."),
          keepStratifiedOutput = ezFrame(
            Type = "logical", DefaultValue = "TRUE",
            Description = "Keep SGB-stratified output rows in full mode."),
          regroupKEGGKO = ezFrame(
            Type = "logical", DefaultValue = "TRUE",
            Description = "Run UniRef90 -> KEGG KO regroup as an extra output table.")
        )
      }
    )
  )
