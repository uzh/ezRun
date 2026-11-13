###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodNfCoreCutAndRun <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("CutAndRun_Result") |> basename()

  fullGenomeSize <- param$ezRef@refFastaFile %>%
    Rsamtools::FaFile() %>%
    GenomeInfoDb::seqlengths() %>%
    sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8) %>% round()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getCutAndRunSampleSheet(input, param)
  write_csv(nfSampleInfo, nfSampleFile)

  blackListFile <- getBlackListFile(input, param)
  prepNFCoreEnv()
  configFile <- writeNextflowLimits(param)
  cmd = paste(
    "nextflow run -c",
    configFile,
    "nf-core/cutandrun",
    "-ansi-log false",
    ## i/o
    "--input",
    nfSampleFile,
    "--outdir",
    outFolder,
    ## genome files
    "--fasta",
    param$ezRef@refFastaFile,
    "--gtf",
    param$ezRef@refFeatureFile,
    "--gene_bed",
    str_replace(
      param$ezRef@refAnnotationFile,
      basename(param$ezRef@refAnnotationFile),
      'genes.bed'
    ),
    "--blacklist",
    blackListFile,
    # if (param[['blacklist']] == "")  "" else paste0("--blacklist ", param[['blacklist']]),
    ## parameters
    "--dt_calc_all_matrix false",
    "--macs_gsize",
    effectiveGenomeSize,
    "--peakcaller",
    param[['peakCaller']],
    "--spikein_genome",
    param[['spikeinGenome']],
    "--normalisation_mode",
    param[['normalization']],
    if (param[['peakStyle']] == 'broad') "" else "--macs2_narrow_peak",
    ## configuration
    "-work-dir work",
    "-profile apptainer",
    "-r",
    param$pipelineVersion,
    param$cmdOptions
  )
  ezSystem(cmd)
  ezSystem(paste('mv', configFile, outFolder))
  baseUrl = file.path(PROJECT_BASE_URL, output$getColumn("CutAndRun_Result"))
  generateFastaFromBedFiles(outFolder, refFile = param$ezRef["refFastaFile"])
  generateAnnotatedPeaks(gtfFile = param$ezRef@refFeatureFile, outFolder)
  jsonFile = writeCutAndRunIgvSession(
    param,
    outFolder,
    jsonFileName = paste0(outFolder, "/igv_session.json"),
    bigwigRelPath = "/04_reporting/igv/",
    baseUrl
  )
  writeNfCoreIgvHtml(
    param,
    jsonFile,
    title = "NfCoreCutAndRun MultiSample Coverage Tracks",
    htmlTemplate = "templates/igvNfCoreTemplate.html",
    htmlFileName = paste0(outFolder, "/igv_session.html")
  )
  sampleNames <- paste(
    nfSampleInfo[["group"]],
    nfSampleInfo[["replicate"]],
    sep = "_R"
  )
  makeRmdReportWrapper(
    htmlPath = paste0(outFolder, "/04_reporting/"),
    rmdFile = "NfCoreCutAndRun.Rmd",
    reportTitle = "NfCoreCutAndRun",
    baseUrlPeaks = paste0(baseUrl, "/04_reporting/igv/", sep = "/"),
    sampleNames
  )

  if (ezIsSpecified(param$keepBams)) {
    keepBams <- param$keepBams
  } else {
    keepBams <- TRUE
  }
  cleanupCarOutFolder(outFolder, keepBams)

  return("Success")
}

EzAppNfCoreCutAndRun <- setRefClass(
  "EzAppNfCoreCutAndRun",
  contains = "EzApp",
  methods = list(
    ## nf-core/cutandrun (ezMethodNfCoreCutAndRun) and the tools of its CITATIONS.md/README that run
    ## with the options the app passes (skip_* left at the pipeline defaults):
    ## - nf-core and Nextflow: the pipeline framework, always.
    ## - Trim Galore/Cutadapt: trimming, always.
    ## - Bowtie 2: alignment to the target and the spike-in genome, always.
    ## - Picard MarkDuplicates: duplicate marking (removal in controls), always.
    ## - BEDTools: spike-in scaled bedGraphs and consensus peaks, always.
    ## - deepTools bamCoverage: the RPKM/CPM/BPM normalised coverage, only when normalization is
    ##   one of these (default Spikein uses BEDTools; None leaves the coverage unscaled).
    ## - SEACR / MACS2: the peak caller(s) in --peakcaller (param peakCaller, app default macs2).
    ## - ChIPpeakAnno: ezRun's nearest-TSS annotation of the peak BED files (generateAnnotatedPeaks).
    ## Not cited (generic BAM utilities, run QC and reports): SAMtools filters (the criteria are
    ## described instead), UCSC bedGraphToBigWig (format conversion), FastQC, preseq, deepTools QC
    ## plots and MultiQC.
    citation = function(param = list()) {
      callers <- tolower(if (ezIsSpecified(param$peakCaller)) param$peakCaller else "macs2")
      norm <- toupper(if (ezIsSpecified(param$normalization)) param$normalization else "SPIKEIN")
      c(
        "Ewels, P.A. et al. The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology 38(3), 276-278 (2020). https://doi.org/10.1038/s41587-020-0439-x",
        "Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35(4), 316-319 (2017). https://doi.org/10.1038/nbt.3820",
        "Trim Galore. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/",
        "Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17(1), 10-12 (2011). https://doi.org/10.14806/ej.17.1.200",
        "Langmead, B. & Salzberg, S.L. Fast gapped-read alignment with Bowtie 2. Nat Methods 9, 357-359 (2012). https://doi.org/10.1038/nmeth.1923",
        "Picard Toolkit. Broad Institute. https://broadinstitute.github.io/picard/",
        "Quinlan, A.R. & Hall, I.M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841-842 (2010). https://doi.org/10.1093/bioinformatics/btq033",
        if (norm %in% c("RPKM", "CPM", "BPM")) "Ramírez, F. et al. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research 44, W160-W165 (2016). https://doi.org/10.1093/nar/gkw257",
        if (grepl("seacr", callers)) "Meers, M.P., Tenenbaum, D. & Henikoff, S. Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics & Chromatin 12, 42 (2019). https://doi.org/10.1186/s13072-019-0287-4",
        if (grepl("macs2", callers)) "Zhang, Y. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biology 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137",
        "Zhu, L.J. et al. ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC Bioinformatics 11, 237 (2010). https://doi.org/10.1186/1471-2105-11-237"
      )
    },
    initialize = function() {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreCutAndRun
      name <<- "EzAppNfCoreCutAndRun"
      ## minimum nf-core parameters
      appDefaults <<- rbind(
        peakCaller = ezFrame(
          Type = "character",
          DefaultValue = "macs2",
          Description = "Select the peak caller for the pipeline"
        ),
        spikeinGenome = ezFrame(
          Type = "character",
          DefaultValue = "macs2",
          Description = "Select the reference for the spike-in genome"
        ),
        normalization = ezFrame(
          Type = "character",
          DefaultValue = "macs2",
          Description = "Select the target read normalization mode"
        ),
        peakStyle = ezFrame(
          Type = "character",
          DefaultValue = "broad",
          Description = "Run MACS2 in broadPeak mode, otherwise in narrowPeak mode"
        ),
        keepBams = ezFrame(
          Type = "logical",
          DefaultValue = FALSE,
          Description = "Should bam files be stored"
        ),
        pipelineVersion = ezFrame(
          Type = "character",
          DefaultValue = '3.2.2',
          Description = "specify pipeline version"
        )
      )
    }
  )
)

##' @description fetch nfcore blacklist files
getBlackListFile <- function(input, param) {
  baseBuildName = param$ezRef@refBuildName |> str_replace("\\.p.*", "")
  basePath = '/srv/GT/databases/nf-core/cutandrun/'
  blackListPath <- paste0(basePath, "/", baseBuildName, "-blacklist.bed")
  if (file.exists(blackListPath)) {
    return(blackListPath)
  } else {
    warning(
      "no blacklist file for: ",
      baseBuildName,
      " available in ",
      basePath
    )
    return("")
  }
}

##' @description get an nf-core/cutandrun-formatted csv file
getCutAndRunSampleSheet <- function(input, param) {
  ## TODO: does not yet support comma-separated file paths
  nfSampleInfo <- ezFrame(
    group = input$getColumn(param$grouping),
    replicate = ezReplicateNumber(input$getColumn(param$grouping)),
    fastq_1 = input$getFullPaths("Read1"),
    fastq_2 = input$getFullPaths("Read2"),
    control = input$getColumn(param$controlColumn)
  )

  return(nfSampleInfo)
}

##' @description clean up NfCoreCutAndRun_result directory
cleanupCarOutFolder <- function(outFolder, keepBams = TRUE) {
  if (!keepBams) {
    bamPath <- paste0(outFolder, "/02_alignment/")
    bamsToDelete <- dir(
      path = bamPath,
      pattern = "*.bam(.bai)?$",
      recursive = TRUE
    )
    file.remove(file.path(bamPath, bamsToDelete))
    ezLog("Deleted bam and bam.bai files form bwa directory.")
  }
  genomePath <- paste0(outFolder, "/04_reporting/igv/")
  filesToDelete <- dir(path = genomePath, pattern = "genome.fa(.fai)?$")
  file.remove(file.path(genomePath, filesToDelete))
  ezLog(
    "Deleted genome.fa and genome.fai files from 04_reporting/igv/ directory."
  )
}

##' @description write IGV session in json format
writeCutAndRunIgvSession <- function(
  param,
  outFolder,
  jsonFileName,
  bigwigRelPath,
  baseUrl
) {
  refBuildName = param$ezRef@refBuildName
  refUrlBase = file.path(REF_HOST, param$ezRef@refBuild)
  fastaUrl = sub(
    "Annotation.*",
    "Sequence/WholeGenomeFasta/genome.fa",
    refUrlBase
  )
  faiUrl = paste0(fastaUrl, ".fai")

  bigwigPath = file.path(outFolder, bigwigRelPath)
  bigwigFiles <- dir(path = bigwigPath, pattern = "*.bigWig$")
  trackNames <- bigwigFiles |> str_replace("\\..*", "")
  tracks <- list()
  tracks[[1]] <- list(type = "sequence")
  for (i in 1:length(bigwigFiles)) {
    tracks[[i + 1]] <- list(
      id = trackNames[i],
      url = paste0(baseUrl, file.path(bigwigRelPath, bigwigFiles[i])),
      format = "bigWig",
      name = trackNames[i]
    )
  }
  tracks[[i + 2]] <- list(
    id = "genes",
    url = file.path(
      REF_HOST,
      param$ezRef@refBuild,
      'Genes/transcripts.only.gtf'
    ),
    format = "gtf",
    type = "annotation",
    name = "genes"
  )

  tracks[[i + 3]] <- list(
    id = "exons",
    url = file.path(REF_HOST, param$ezRef@refBuild, 'Genes/genes.bed'),
    format = "bed",
    type = "annotation",
    name = "exons"
  )
  jsonLines <- list(
    version = "3.5.3",
    showSampleNames = FALSE,
    reference = list(id = refBuildName, fastaUrl = fastaUrl, indexURL = faiUrl),
    tracks = tracks
  )
  jsonFile <- rjson::toJSON(jsonLines, indent = 5, method = "C")
  write(jsonFile, jsonFileName)
  return(jsonFile)
}

##' @description write html wrapper for IGV session
writeNfCoreIgvHtml = function(
  param,
  jsonFile,
  title,
  htmlTemplate,
  htmlFileName
) {
  htmlLines = readLines(system.file(
    htmlTemplate,
    package = "ezRun",
    mustWork = TRUE
  ))
  htmlLines = gsub("TITLE", title, htmlLines)
  htmlLines = gsub("IGV_JSON_CONTENT", jsonFile, htmlLines)
  writeLines(htmlLines, htmlFileName)
}


##' @description generate fasta files from BED files
generateFastaFromBedFiles <- function(outFolder, refFile) {
  bedFilePath <- paste0(outFolder, "/04_reporting/igv/")
  bedFileNames <- dir(
    path = bedFilePath,
    pattern = ".bed$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (name in bedFileNames) {
    peakSeqFile = paste0(name, "_peaks.fa")
    cmd = paste(
      "bedtools",
      " getfasta -fi",
      refFile,
      "-bed",
      name,
      " -name -fo ",
      peakSeqFile
    )
    ezSystem(cmd)
  }
}

##' @description annotate peaks in BED files
generateAnnotatedPeaks <- function(gtfFile, outFolder) {
  require(ChIPpeakAnno)
  require(GenomicRanges)
  require(rtracklayer)
  gtf <- rtracklayer::import(gtfFile)
  if ('gene' %in% unique(gtf$type)) {
    idx = gtf$type == 'gene'
  } else if ('transcript' %in% unique(gtf$type)) {
    idx = gtf$type == 'transcript'
  } else if ('start_codon' %in% unique(gtf$type)) {
    idx = gtf$type == 'start_codon'
  } else {
    ezLog('gtf is incompatabible. Peak annotation skipped!')
    return(NULL)
  }
  gtf = gtf[idx]
  if (grepl('gtf$', gtfFile)) {
    names_gtf = make.unique(gtf$'gene_id')
  } else {
    names_gtf = make.unique(gtf$'ID')
  }
  names(gtf) = names_gtf

  bedFilePath <- paste0(outFolder, "/04_reporting/igv/")
  bedFileNames <- dir(
    path = bedFilePath,
    pattern = ".bed$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (bedFile in bedFileNames) {
    peakAnnFile = paste0(sub("\\.bed$", "", bedFile), ".xlsx")
    myPeaks = ezRead.table(bedFile, row.names = NULL, header = F)
    colnames(myPeaks)[1:6] <- c(
      "chrom",
      "start",
      "end",
      "name",
      "score",
      "strand"
    )
    peaksRD = makeGRangesFromDataFrame(
      myPeaks,
      keep.extra.columns = TRUE,
      start.field = "start",
      end.field = "end",
      seqnames.field = "chrom"
    )
    annotatedPeaks <- annotatePeakInBatch(
      peaksRD,
      AnnotationData = gtf,
      output = 'nearestStart',
      multiple = FALSE,
      FeatureLocForDistance = 'TSS'
    )
    annotatedPeaks <- as.data.frame(annotatedPeaks)
    annotatedPeaks <- annotatedPeaks %>%
      rename("feature_start" = "start_position", "feature_end" = "end_position")
    writexl::write_xlsx(annotatedPeaks, peakAnnFile)
  }
}

##' @description write HTML report
makeRmdReportWrapper <- function(
  htmlPath,
  rmdFile,
  reportTitle,
  baseUrlPeaks,
  sampleNames
) {
  cd = getwd()
  setwdNew(htmlPath)
  makeRmdReport(
    baseUrlPeaks = baseUrlPeaks,
    sampleNames = sampleNames,
    rmdFile = rmdFile,
    reportTitle = reportTitle,
    selfContained = TRUE
  )
  setwd(cd)
}
