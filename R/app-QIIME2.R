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
  require(whisker)
  dataset = input$meta
  sampleNames = input$getNames()
  isPaired <- param$paired

  file1PathInDataset <- input$getFullPaths("Read1")
  if (isPaired) {
    file2PathInDataset <- input$getFullPaths("Read2")
  }

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
    extra_col <- paste("/srv/gstore/projects", manifest[, 2], sep = "/")
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
      extra_col1 <- paste("/srv/gstore/projects", manifest[, 2], sep = "/")
      extra_col2 <- paste("/srv/gstore/projects", manifest[, 3], sep = "/")
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

  # Resolve database dropdown -> reference seq/tax artifacts. The dropdown
  # value is the bare <name> of a (<name>-seqs.qza, <name>-tax.qza) pair under
  # QIIME2_DB_ROOT, auto-detected by QIIME2App.rb#qiime2_db_choices. We resolve
  # both paths here and pass them straight to the batch as DB_SEQS / DB_TAX.
  db_seqs <- file.path(QIIME2_DB_ROOT, paste0(param$database, "-seqs.qza"))
  db_tax  <- file.path(QIIME2_DB_ROOT, paste0(param$database, "-tax.qza"))
  if (!file.exists(db_seqs) || !file.exists(db_tax)) {
    stop(
      "Reference DB pair not found for label '", param$database, "'.\n",
      "  expected: ", db_seqs, "\n",
      "  expected: ", db_tax
    )
  }

  classifier_path <- if (is.null(param$classifier_path) ||
    !nzchar(param$classifier_path)) {
    # NONE sentinel keeps `if [ "$CLASSIFIER_PATH" = "NONE" ]` simple in bash
    "NONE"
  } else {
    param$classifier_path
  }
  run_fastp_flag <- if (isTRUE(param$run_fastp)) "true" else "false"
  run_picrust2_flag <- if (isTRUE(param$run_picrust2)) "true" else "false"

  placeholder_values <- list(
    TRIM_LEFT_F = param$trim_left_f,
    TRIM_LEFT_R = param$trim_left_r,
    TRUNC_LEN_F = param$truncate_len_f,
    TRUNC_LEN_R = param$truncate_len_r,
    SAMPLING_DEPTH = param$sampling_depth,
    MAX_RAREFACTION_DEPTH = param$max_rarefaction_depth,
    MIN_FREQ = param$min_freq,
    MIN_SAMPLES = param$min_samples,
    DB_SEQS = db_seqs,
    DB_TAX  = db_tax,
    PRIMER1 = param$forward_primer,
    PRIMER2 = param$reverse_primer,
    CLASSIFIER_MIN_LEN = param$classifier_min_len,
    CLASSIFIER_MAX_LEN = param$classifier_max_len,
    CLASSIFIER_PATH = classifier_path,
    RUN_FASTP = run_fastp_flag,
    RUN_PICRUST2 = run_picrust2_flag,
    CORES = param$cores
  )

  # Render {{...}} mustache placeholders with whisker. We deliberately do not
  # use chained sed here: an unbounded `s/DB/.../g` would also rewrite DB_SEQS /
  # DB_TAX, and `s/CORES/16/g` would rewrite $CORES into $16. Mustache braces
  # make each placeholder unambiguous.
  tmpl_name <- if (isPaired) {
    "unifiedQIIME2Workflow.pairedend.batch"
  } else {
    "unifiedQIIME2Workflow.singleend.batch"
  }
  tmpl_path <- system.file("scripts", "QIIME2", tmpl_name, package = "ezRun")
  if (!nzchar(tmpl_path)) {
    stop("QIIME2 batch template not found in ezRun (inst/scripts/QIIME2/", tmpl_name, ")")
  }
  out_batch <- tmpl_name
  tmpl <- paste(readLines(tmpl_path), collapse = "\n")
  rendered <- whisker::whisker.render(tmpl, placeholder_values)
  writeLines(rendered, out_batch)
  ezSystem(paste("sh", out_batch))

  physeq <- qza_to_phyloseq(
    features = "table-filtered.qza",
    tree = "rooted-tree.qza",
    "taxonomy.qza",
    metadata = "sample_metadata.tsv"
  )
  physeq_unfiltered <- qza_to_phyloseq(
    features = "table.qza",
    tree = "rooted-tree.qza",
    "taxonomy.qza",
    metadata = "sample_metadata.tsv"
  )
  dada2 <- read_qza("dada2_denoising_stats.qza")
  table_unfiltered <- read_qza("table.qza")
  rep_seqs <- read_qza("rep-seqs-filtered.qza")
  asv_taxonomy <- read_qza("taxonomy.qza")
  asv_taxonomy_df <- asv_taxonomy$data
  tree <- read_qza("rooted-tree.qza")
  shannon <- read_qza("core-metrics-results/shannon_vector.qza")
  evenness <- read_qza("core-metrics-results/evenness_vector.qza")
  richness <- read_qza("core-metrics-results/observed_features_vector.qza")
  faith_pd <- read_qza("core-metrics-results/faith_pd_vector.qza")
  jaccard <- read_qza("core-metrics-results/jaccard_pcoa_results.qza")
  jaccardmatrix <- read_qza("core-metrics-results/jaccard_distance_matrix.qza")
  bray <- read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
  braymatrix <- read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
  unifrac <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")
  unifracmatrix <- read_qza(
    "core-metrics-results/weighted_unifrac_distance_matrix.qza"
  )
  metadata <- read.table("sample_metadata.tsv", header = TRUE)

  # Export rep-seqs FASTA so the Rmd can run a primer-orientation check
  rep_seqs_fasta <- "dada2_rep_seqs_filtered.fasta"
  rep_seqs_seqs <- rep_seqs$data
  if (inherits(rep_seqs_seqs, "DNAStringSet") ||
    inherits(rep_seqs_seqs, "XStringSet")) {
    Biostrings::writeXStringSet(rep_seqs_seqs, filepath = rep_seqs_fasta)
  } else if (is.character(rep_seqs_seqs)) {
    writeLines(
      c(rbind(paste0(">", names(rep_seqs_seqs)), unname(rep_seqs_seqs))),
      con = rep_seqs_fasta
    )
  }

  if (isTRUE(param$run_picrust2)) {
    picrust2_ko <- if (file.exists("ko_abundance.tsv")) {
      read.table("ko_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE)
    } else {
      NULL
    }
    picrust2_ec <- if (file.exists("ec_abundance.tsv")) {
      read.table("ec_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE)
    } else {
      NULL
    }
    picrust2_pathway <- if (file.exists("pathway_abundance.tsv")) {
      read.table("pathway_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE)
    } else {
      NULL
    }
  } else {
    picrust2_ko <- NULL
    picrust2_ec <- NULL
    picrust2_pathway <- NULL
  }

  fastp_reports <- if (isTRUE(param$run_fastp) && dir.exists("fastp_reports")) {
    list.files("fastp_reports", full.names = TRUE)
  } else {
    NULL
  }

  primer_check <- list(
    forward_primer = param$forward_primer,
    reverse_primer = param$reverse_primer,
    rep_seqs_fasta = rep_seqs_fasta
  )

  library_sizes_pre <- if (!is.null(table_unfiltered$data)) {
    colSums(as.matrix(table_unfiltered$data))
  } else {
    NULL
  }
  library_sizes_post <- phyloseq::sample_sums(physeq)

  read1_paths <- input$getFullPaths("Read1")
  raw_reads_path <- if (length(read1_paths) > 0) {
    unique(dirname(read1_paths))[1]
  } else {
    NULL
  }

  setwdNew(basename(output$getColumn("ResultDir")))
  makeRmdReport(
    param = param,
    physeq = physeq,
    physeq_unfiltered = physeq_unfiltered,
    metadata = metadata,
    dada2 = dada2,
    shannon = shannon,
    evenness = evenness,
    richness = richness,
    faith_pd = faith_pd,
    jaccard = jaccard,
    jaccardmatrix = jaccardmatrix,
    bray = bray,
    braymatrix = braymatrix,
    unifrac = unifrac,
    unifracmatrix = unifracmatrix,
    tree = tree,
    rep_seqs = rep_seqs,
    asv_taxonomy_df = asv_taxonomy_df,
    fastp_reports = fastp_reports,
    primer_check = primer_check,
    library_sizes_pre = library_sizes_pre,
    library_sizes_post = library_sizes_post,
    ko_abundance = picrust2_ko,
    ec_abundance = picrust2_ec,
    pathway_abundance = picrust2_pathway,
    run_picrust2 = isTRUE(param$run_picrust2),
    run_fastp = isTRUE(param$run_fastp),
    grouping = param$grouping,
    raw_reads_path = raw_reads_path,
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
          trim_left_f = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "DADA2 trim-left for forward reads (5' bases to remove; usually = forward primer length)"
          ),
          trim_left_r = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "DADA2 trim-left for reverse reads (5' bases to remove; usually = reverse primer length). Used only when paired"
          ),
          truncate_len_f = ezFrame(
            Type = "integer",
            DefaultValue = "150",
            Description = "DADA2 truncation length for forward reads"
          ),
          truncate_len_r = ezFrame(
            Type = "integer",
            DefaultValue = "150",
            Description = "DADA2 truncation length for reverse reads. Ignored when paired = false"
          ),
          sampling_depth = ezFrame(
            Type = "integer",
            DefaultValue = "6000",
            Description = "Total frequency that each sample should be rarefied to for alpha/beta diversity"
          ),
          max_rarefaction_depth = ezFrame(
            Type = "integer",
            DefaultValue = "4000",
            Description = "Maximum depth for alpha rarefaction curves"
          ),
          min_freq = ezFrame(
            Type = "integer",
            DefaultValue = "1",
            Description = "Minimum total frequency a feature must reach to be retained in differential abundance filtering"
          ),
          min_samples = ezFrame(
            Type = "integer",
            DefaultValue = "1",
            Description = "Minimum number of samples a feature must appear in to be retained for differential abundance"
          ),
          database = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Reference database label — bare <name> of a (<name>-seqs.qza, <name>-tax.qza) pair under QIIME2_DB_ROOT. Resolved to DB_SEQS/DB_TAX paths and passed to the batch"
          ),
          primer = ezFrame(
            Type = "character",
            DefaultValue = "V3-V4",
            Description = "16S region label that drives the default primer pair"
          ),
          forward_primer = ezFrame(
            Type = "character",
            DefaultValue = "CCTACGGGNGGCWGCAG",
            Description = "Forward primer used for classifier extract-reads"
          ),
          reverse_primer = ezFrame(
            Type = "character",
            DefaultValue = "GACTACHVGGGTATCTAATCC",
            Description = "Reverse primer used for classifier extract-reads"
          ),
          classifier_min_len = ezFrame(
            Type = "integer",
            DefaultValue = "350",
            Description = "qiime feature-classifier extract-reads --p-min-length"
          ),
          classifier_max_len = ezFrame(
            Type = "integer",
            DefaultValue = "550",
            Description = "qiime feature-classifier extract-reads --p-max-length"
          ),
          classifier_path = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Absolute path to a pre-built classifier .qza under CLASSIFIER_ROOT; empty/NONE forces on-the-fly training"
          ),
          run_fastp = ezFrame(
            Type = "logical",
            DefaultValue = "true",
            Description = "Run fastp adapter trimming before DADA2 (adapter contamination is frequent on amplicon runs)"
          ),
          run_picrust2 = ezFrame(
            Type = "logical",
            DefaultValue = "false",
            Description = "Run standalone picrust2_pipeline.py after taxonomy filtering"
          )
        )
      }
    )
  )
