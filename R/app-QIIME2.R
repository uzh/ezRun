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
  # Use bash, not sh: the fastp block needs ANSI-C $'\t' tab and `echo -e`,
  # neither supported by dash (Debian's /bin/sh).
  ezSystem(paste("bash", out_batch))

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
  # check.names = FALSE keeps "sample-id" as-is (else it becomes "sample.id")
  metadata <- read.table("sample_metadata.tsv", header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)

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
    # row.names = 1 puts the non-numeric feature-id column (e.g. "function"
    # with "ko:K00001") into rownames, so the returned data.frame is purely
    # numeric. Without this, as.matrix() downstream coerces everything to
    # character and ALDEx2/metagenomeSeq fail with "non-numeric argument".
    picrust2_ko <- if (file.exists("ko_abundance.tsv")) {
      read.table("ko_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
    } else {
      NULL
    }
    picrust2_ec <- if (file.exists("ec_abundance.tsv")) {
      read.table("ec_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
    } else {
      NULL
    }
    picrust2_pathway <- if (file.exists("pathway_abundance.tsv")) {
      read.table("pathway_abundance.tsv", header = TRUE, sep = "\t",
                 check.names = FALSE, row.names = 1)
    } else {
      NULL
    }
    # Absolute path to the stratified per-ASV contribution table. Must be
    # absolute: fastp_reports/contrib are resolved at render time, when the CWD
    # has moved into Results_Folder (setwdNew below), so a relative name would
    # no longer point at the scratch-root file. The Rmd reads this at render
    # time (scratch still alive) to build the "ASV contributors" table.
    picrust2_contrib <- if (file.exists("picrust2_contrib.tsv")) {
      normalizePath("picrust2_contrib.tsv")
    } else {
      NULL
    }
  } else {
    picrust2_ko <- NULL
    picrust2_ec <- NULL
    picrust2_pathway <- NULL
    picrust2_contrib <- NULL
  }

  fastp_reports <- if (isTRUE(param$run_fastp) && dir.exists("fastp_out")) {
    list.files("fastp_out", pattern = "\\.(json|html)$",
               full.names = TRUE, ignore.case = TRUE)
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
    picrust2_contrib = picrust2_contrib,
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
      ## Batch (inst/scripts/QIIME2/unifiedQIIME2Workflow.*.batch): QIIME 2, DADA2, q2-feature-classifier
      ## (extract-reads + naive Bayes, or a prebuilt classifier), MAFFT + FastTree (align-to-tree-mafft-
      ## fasttree) and the UniFrac and Faith's PD metrics of core-metrics-phylogenetic (weighted UniFrac
      ## PERMANOVA in the batch, Faith's PD tests in the report) unconditional; fastp unless run_fastp is false (step 1, default true); PICRUSt2 only
      ## when run_picrust2 is true (step 18). Report (Qiime2Report.Rmd), unconditional (needs >= 2 groups,
      ## data-dependent): phyloseq (filtering, rarefaction, richness), vegan adonis2 PERMANOVA + betadisper,
      ## pairwiseAdonis, ANCOM-BC2, DESeq2, edgeR, LEfSe via microbiomeMarker, ALDEx2; metagenomeSeq only
      ## in the PICRUSt2 section. Reference: SILVA / Greengenes from the database label, or from the
      ## prebuilt classifier's file name when classifier_path is set (the classifier then defines the
      ## taxonomy, the database is not used).
      citation = function(param = list()) {
        cp <- as.character(param$classifier_path)
        ref <- if (length(cp) == 1 && nzchar(cp) && cp != "NONE") basename(cp) else as.character(param$database)
        ref <- tolower(paste(ref, collapse = " "))
        picrust <- isTRUE(as.logical(param$run_picrust2))
        c(
          "Bolyen, E. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37, 852-857 (2019). https://doi.org/10.1038/s41587-019-0209-9",
          if (!isFALSE(as.logical(param$run_fastp))) "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560",
          "Callahan, B.J. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods 13, 581-583 (2016). https://doi.org/10.1038/nmeth.3869",
          "Bokulich, N.A. et al. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. Microbiome 6, 90 (2018). https://doi.org/10.1186/s40168-018-0470-z",
          if (grepl("silva", ref)) "Quast, C. et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research 41(D1), D590-D596 (2013). https://doi.org/10.1093/nar/gks1219",
          if (grepl("greengenes|^gg[-_]", ref)) "McDonald, D. et al. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. The ISME Journal 6, 610-618 (2012). https://doi.org/10.1038/ismej.2011.139",
          "Katoh, K. & Standley, D.M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30(4), 772-780 (2013). https://doi.org/10.1093/molbev/mst010",
          "Price, M.N., Dehal, P.S. & Arkin, A.P. FastTree 2 - approximately maximum-likelihood trees for large alignments. PLoS ONE 5(3), e9490 (2010). https://doi.org/10.1371/journal.pone.0009490",
          "Lozupone, C. & Knight, R. UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology 71(12), 8228-8235 (2005). https://doi.org/10.1128/AEM.71.12.8228-8235.2005",
          "Faith, D.P. Conservation evaluation and phylogenetic diversity. Biological Conservation 61(1), 1-10 (1992). https://doi.org/10.1016/0006-3207(92)91201-3",
          if (picrust) "Douglas, G.M. et al. PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38, 685-688 (2020). https://doi.org/10.1038/s41587-020-0548-6",
          "McMurdie, P.J. & Holmes, S. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE 8(4), e61217 (2013). https://doi.org/10.1371/journal.pone.0061217",
          "Anderson, M.J. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26(1), 32-46 (2001). https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x",
          "Oksanen, J. et al. vegan: Community Ecology Package. R package. https://CRAN.R-project.org/package=vegan",
          "Martinez Arbizu, P. pairwiseAdonis: Pairwise multilevel comparison using adonis. R package. https://github.com/pmartinezarbizu/pairwiseAdonis",
          "Lin, H. & Peddada, S.D. Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods 21, 83-91 (2024). https://doi.org/10.1038/s41592-023-02092-7",
          "Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8",
          "Robinson, M.D., McCarthy, D.J. & Smyth, G.K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26(1), 139-140 (2010). https://doi.org/10.1093/bioinformatics/btp616",
          "Segata, N. et al. Metagenomic biomarker discovery and explanation. Genome Biology 12, R60 (2011). https://doi.org/10.1186/gb-2011-12-6-r60",
          "Cao, Y. et al. microbiomeMarker: an R/Bioconductor package for microbiome marker identification and visualization. Bioinformatics 38(16), 4027-4029 (2022). https://doi.org/10.1093/bioinformatics/btac438",
          "Fernandes, A.D., Macklaim, J.M., Linn, T.G., Reid, G. & Gloor, G.B. ANOVA-like differential expression (ALDEx) analysis for mixed population RNA-Seq. PLoS ONE 8(7), e67019 (2013). https://doi.org/10.1371/journal.pone.0067019",
          "Fernandes, A.D. et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome 2, 15 (2014). https://doi.org/10.1186/2049-2618-2-15",
          if (picrust) "Paulson, J.N., Stine, O.C., Bravo, H.C. & Pop, M. Differential abundance analysis for microbial marker-gene surveys. Nature Methods 10, 1200-1202 (2013). https://doi.org/10.1038/nmeth.2658"
        )
      },
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
