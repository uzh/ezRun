###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMetaPhlAn = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  sampleName <- input$getNames()
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  indexName <- as.character(param$metaphlanIndex)
  if (!nzchar(indexName)) {
    stop("metaphlanIndex is empty: no database selected.")
  }
  dbDir <- "/srv/GT/databases/metaphlan_databases"
  if (!file.exists(file.path(dbDir, paste0(indexName, ".pkl")))) {
    stop(sprintf("MetaPhlAn .pkl not found for index '%s' under %s.", indexName, dbDir))
  }

  outProfile <- paste0(sampleName, "_metaphlan.txt")
  outMapout <- paste0(sampleName, ".bowtie2.bz2")
  outLog <- paste0(sampleName, ".metaphlan.log")

  if (param$paired) {
    read1 <- trimmedInput$getColumn("Read1")
    read2 <- trimmedInput$getColumn("Read2")
    readArg <- paste0(read1, ",", read2)
  } else {
    readArg <- trimmedInput$getColumn("Read1")
  }

  ## Backwards-compatible: missing param is treated as TRUE so existing
  ## dataset definitions keep producing count-augmented profiles.
  estimateCounts <- if (is.null(param$estimateReadCounts)) TRUE
                    else isTRUE(as.logical(param$estimateReadCounts))
  analysisType   <- if (estimateCounts) "-t rel_ab_w_read_stats" else ""

  cmd <- paste(
    "metaphlan",
    readArg,
    "--db_dir", dbDir,
    "--index", indexName,
    "--input_type fastq",
    analysisType,
    "--nproc", ezThreads(),
    "--tmp_dir .",
    "--mapout", outMapout,
    "-o", outProfile,
    param$cmdOptions,
    "1>", outLog, "2>&1"
  )
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMetaPhlAn()
##' @templateVar htmlArg )
##' @description Use this reference class to run MetaPhlAn taxonomic profiling
EzAppMetaPhlAn <-
  setRefClass(
    "EzAppMetaPhlAn",
    contains = "EzApp",
    methods = list(
      ## MetaPhlAn defaults checked against `metaphlan --help` of Tools/MetaPhlAn/4.2.4,
      ## the module MetaPhlAnApp.rb loads.
      ## fastp + MetaPhlAn (Bowtie 2 inside MetaPhlAn) unconditional.
      citation = function() {
        c(
          "Blanco-Míguez, A. et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nature Biotechnology 41, 1633-1644 (2023). https://doi.org/10.1038/s41587-023-01688-w",
          "Langmead, B. & Salzberg, S.L. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012). https://doi.org/10.1038/nmeth.1923",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## ezMethodMetaPhlAn -> ezMethodFastpTrim (app-metaphlan.R:15)
          "Reads were preprocessed with fastp (ezRun FastpTrim with the trimming parameters of the job) and the trimmed reads, not the raw reads, were profiled with MetaPhlAn.",
          ## ezMethodFastpTrim adapter FASTA (app-trim.R:140-178)
          if (isTRUE(as.logical(param$trimAdapter)) && !isTRUE(as.logical(param$onlyAdapterFromDataset))) "fastp trimmed adapters given as a FASTA file of Illumina adapter sequences (the FGCZ Trimmomatic adapter set) plus any adapter named in the dataset.",
          ## ezMethodMetaPhlAn readArg (app-metaphlan.R:30-33)
          if (isTRUE(as.logical(param$paired))) "The two mates were given to MetaPhlAn as one comma-separated input and mapped as independent reads, which is how MetaPhlAn treats paired data.",
          ## ezMethodMetaPhlAn analysisType (app-metaphlan.R:40-42)
          if (isTRUE(as.logical(param$estimateReadCounts))) "MetaPhlAn ran with -t rel_ab_w_read_stats, so the profile also gives the estimated number of reads per clade." else if (isFALSE(as.logical(param$estimateReadCounts))) "MetaPhlAn ran in its default relative-abundance mode (-t rel_ab).",
          "ezRun left the MetaPhlAn profiling settings at their defaults unless cmdOptions set them: reads shorter than 70 nt ignored, minimum mapping quality 5 for short reads, and clade abundance as the marker average truncated at the 0.2 quantile (--stat tavg_g, --stat_q 0.2)."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMetaPhlAn
        name <<- "EzAppMetaPhlAn"
        appDefaults <<- rbind(
          metaphlanIndex = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "MetaPhlAn bowtie2 index basename under /srv/GT/databases/metaphlan_databases/ (e.g. mpa_vJan25_CHOCOPhlAnSGB_202503)."
          ),
          estimateReadCounts = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "When TRUE, pass -t rel_ab_w_read_stats so the profile carries the estimated_number_of_reads_from_the_clade column. Required by count-based downstream DA (ALDEx2 / ANCOM-BC2)."
          )
        )
      }
    )
  )
