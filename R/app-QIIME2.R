###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodQIIME2 = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  require(rmarkdown)
  require(data.table)
  require(here)
  require(tidyverse)
  require(qiime2R)
  dataset = input$meta
  sampleNames = input$getNames()
  isPaired <- param$paired
  ### read fastq files and prepare inputs
  file1PathInDataset <- input$getFullPaths("Read1")
  if (isPaired) {
    file2PathInDataset <- input$getFullPaths("Read2")
  }
  ###create sample metadata and manifest files
  if (!file.exists("sample_metadata.tsv")) {
    file.create("sample_metadata.tsv")
    sample_metadata <- data.frame(dataset[, grepl(
      param$grouping,
      names(dataset)
    )])
    sample_metadata$sample_id <- rownames(dataset)
    sample_metadata <- sample_metadata %>% select(sample_id, everything())
    setnames(sample_metadata, "sample_id", "sample-id")
    setnames(sample_metadata, colnames(sample_metadata)[2], "Group")
    write_tsv(sample_metadata, file = "sample_metadata.tsv")
  } else {
    print("The file exists")
  }

  if (!file.exists("manifest.tsv")) {
    file.create("manifest.tsv")
    manifest <- data.frame(dataset[, grepl("Read1", names(dataset))])
    manifest$sample_id <- rownames(dataset)
    manifest <- manifest[, c(2, 1)]
    #extra_col <- paste("/srv/gstore/projects", manifest[,2], sep="/")
    extra_col <- paste(
      "/srv/GT/analysis/course_sushi/public/gstore/projects",
      manifest[, 2],
      sep = "/"
    )
    manifest <- cbind(manifest, extra_col)
    manifest <- manifest[, c(1, 3)]
    setnames(manifest, "sample_id", "sample-id")
    setnames(manifest, colnames(manifest)[2], "absolute-filepath")
    write_tsv(manifest, file = "manifest.tsv")
    if (isPaired) {
      manifest1 <- data.frame(dataset[, grepl("Read1", names(dataset))])
      manifest1$sample_id <- rownames(dataset)
      manifest2 <- data.frame(dataset[, grepl("Read2", names(dataset))])
      manifest2$sample_id <- rownames(dataset)
      manifest <- merge(manifest1, manifest2, by = "sample_id")
      manifest <- manifest[, c(1, 2, 3)]
      #extra_col1 <- paste("/srv/gstore/projects", manifest[,2], sep="/")
      #extra_col2 <- paste("/srv/gstore/projects", manifest[,3], sep="/")
      extra_col1 <- paste(
        "/srv/GT/analysis/course_sushi/public/gstore/projects",
        manifest[, 2],
        sep = "/"
      )
      extra_col2 <- paste(
        "/srv/GT/analysis/course_sushi/public/gstore/projects",
        manifest[, 3],
        sep = "/"
      )
      manifest <- cbind(manifest, extra_col1, extra_col2)
      manifest <- manifest[, c(1, 4, 5)]
      setnames(manifest, "sample_id", "sample-id")
      setnames(manifest, colnames(manifest)[2], "forward-absolute-filepath")
      setnames(manifest, colnames(manifest)[3], "reverse-absolute-filepath")
      write_tsv(manifest, file = "manifest.tsv")
    }
  } else {
    print("The file exists")
  }

  updateBatchCmd1 <- paste0(
    "sed -e s/\"TRIM_LEFT\"/",
    param$trim_left,
    "/g",
    " -e s/\"TRUNC_LEN\"/",
    param$truncate_len,
    "/g",
    " -e s/\"SAMPLING_DEPTH\"/",
    param$sampling_depth,
    "/g ",
    " -e s/\"MAX_RAREFACTION_DEPTH\"/",
    param$max_rarefaction_depth,
    "/g ",
    " -e s/\"MIN_FREQ\"/",
    param$min_freq,
    "/g ",
    " -e s/\"MIN_SAMPLES\"/",
    param$min_samples,
    "/g ",
    " -e s/\"DB\"/",
    param$database,
    "/g ",
    " -e s/\"PRIMER1\"/",
    param$forward_primer,
    "/g ",
    " -e s/\"PRIMER2\"/",
    param$reverse_primer,
    "/g ",
    file.path(METAGENOMICS_ROOT, UNIFIED_QIIME2_WORKFLOW_SINGLEEND),
    " > ",
    UNIFIED_QIIME2_WORKFLOW_SINGLEEND
  )
  if (isPaired) {
    updateBatchCmd1 <- paste0(
      "sed -e s/\"TRIM_LEFT\"/",
      param$trim_left,
      "/g",
      " -e s/\"TRUNC_LEN\"/",
      param$truncate_len,
      "/g",
      " -e s/\"SAMPLING_DEPTH\"/",
      param$sampling_depth,
      "/g ",
      " -e s/\"MAX_RAREFACTION_DEPTH\"/",
      param$max_rarefaction_depth,
      "/g ",
      " -e s/\"MIN_FREQ\"/",
      param$min_freq,
      "/g ",
      " -e s/\"MIN_SAMPLES\"/",
      param$min_samples,
      "/g ",
      " -e s/\"DB\"/",
      param$database,
      "/g ",
      " -e s/\"PRIMER1\"/",
      param$forward_primer,
      "/g ",
      " -e s/\"PRIMER2\"/",
      param$reverse_primer,
      "/g ",
      file.path(METAGENOMICS_ROOT, UNIFIED_QIIME2_WORKFLOW_PAIREDEND),
      " > ",
      UNIFIED_QIIME2_WORKFLOW_PAIREDEND
    )
  }
  ezSystem(updateBatchCmd1)

  cmdQIIME2 = paste("sh", UNIFIED_QIIME2_WORKFLOW_SINGLEEND)
  if (isPaired) {
    cmdQIIME2 = paste("sh", UNIFIED_QIIME2_WORKFLOW_PAIREDEND)
  }
  ezSystem(cmdQIIME2)

  physeq <- qza_to_phyloseq(
    features = "table.qza",
    tree = "rooted-tree.qza",
    "taxonomy.qza",
    metadata = "sample_metadata.tsv"
  )
  dada2 <- read_qza("dada2_denoising_stats.qza")
  shannon <- read_qza("core-metrics-results/shannon_vector.qza")
  evenness <- read_qza("core-metrics-results/evenness_vector.qza")
  richness <- read_qza("core-metrics-results/observed_features_vector.qza")
  jaccard <- read_qza("core-metrics-results/jaccard_pcoa_results.qza")
  jaccardmatrix <- read_qza("core-metrics-results/jaccard_distance_matrix.qza")
  bray <- read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
  braymatrix <- read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
  unifrac <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")
  unifracmatrix <- read_qza(
    "core-metrics-results/weighted_unifrac_distance_matrix.qza"
  )
  metadata <- read.table("sample_metadata.tsv", header = TRUE)

  setwdNew(basename(output$getColumn("ResultDir")))
  ## Copy the style files and templates
  makeRmdReport(
    param = param,
    physeq = physeq,
    metadata = metadata,
    dada2 = dada2,
    shannon = shannon,
    evenness = evenness,
    richness = richness,
    jaccard = jaccard,
    jaccardmatrix = jaccardmatrix,
    bray = bray,
    braymatrix = braymatrix,
    unifrac = unifrac,
    unifracmatrix = unifracmatrix,
    output = output,
    rmdFile = "Qiime2Report.Rmd",
    reportTitle = "QIIME2 Static Report"
  )

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodQIIME2()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppQIIME2 <-
  setRefClass(
    "ezMethodQIIME2",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodQIIME2
        name <<- "ezMethodQIIME2"
        appDefaults <<- rbind(
          trim_left = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "Position at which sequences should be trimmed due to low quality"
          ),
          truncate_len = ezFrame(
            Type = "integer",
            DefaultValue = "150",
            Description = "Position at which sequences should be truncated due to decrease in quality"
          ),
          sampling_depth = ezFrame(
            Type = "integer",
            DefaultValue = "1000",
            Description = "Total frequency that each sample should be rarefied to"
          )
        )
      }
    )
  )
