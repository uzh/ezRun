###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodNfCoreSmRnaSeq <- function(input = NA, output = NA, param = NA) {
  refbuild = param$refBuild
  outFolder = output$getColumn("smRNASeq_Result") |> basename()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getSmRnaSeqSampleSheet(input)
  write_csv(nfSampleInfo, nfSampleFile)
  prepNFCoreEnv()
  configFile <- writeNextflowLimits(param)
  cmd = paste(
    "nextflow run nf-core/smrnaseq",
    ## i/o
    "--input",
    nfSampleFile,
    "--outdir",
    outFolder,
    ## genome
    "--genome",
    param[['referenceGenome']],
    ## parameters
    "--mirtrace_species",
    param[['mirtraceSpecies']],
    ## configuration
    "-work-dir work",
    "-profile apptainer",
    "-r",
    param$pipelineVersion,
    "-c",
    configFile,
    param$cmdOptions
  )
  ezSystem(cmd)
  ezSystem(paste('mv', configFile, outFolder))
  writePerSampleSmRnaCountFiles(
    nfSampleInfo,
    countDir = paste0(outFolder, "/mirna_quant/mirtop/")
  )

  return("Success")
}


EzAppNfCoreSmRnaSeq <- setRefClass(
  "EzAppNfCoreSmnRnaSeq",
  contains = "EzApp",
  methods = list(
    ## nf-core/smrnaseq (checked against the 2.4.0 code, the app's default pipelineVersion). nf-core,
    ## Nextflow, the pipeline, Bowtie 1 (mature/hairpin/genome mapping) and mirtop (the miRNA counts the
    ## app splits per sample) unconditional. Gated on the flags a user can add in cmdOptions (pipeline
    ## defaults otherwise): fastp unless --skip_fastp, miRBase (mature/hairpin FASTA and GFF3) unless
    ## --mirgenedb, miRDeep2 unless --skip_mirdeep (it needs the genome FASTA, which --genome always
    ## gives). MirGeneDB with --mirgenedb (nextflow.config mirgenedb_* params; the user supplies the
    ## files, so the release is data: MirGeneDB 3.0, the release offered for download since 2024, is
    ## cited). --with_umi runs UMICollapse and then UMI-tools extract (workflows/smrnaseq.nf, `if
    ## (params.with_umi)`). --filter_contamination runs the contaminant filter only on the databases
    ## given (--rrna/--trna/--cdna/--ncrna/--pirna/--other_contamination; none comes from iGenomes):
    ## Bowtie 2 for any of them, BLAT (removing hairpin-like sequences from the database) for cdna,
    ## ncrna, pirna or other_contamination (subworkflows/local/contaminant_filter). Not cited (read
    ## collapsing, run QC and reports): seqcluster collapse, samtools, FastQC, miRTrace, edgeR (the
    ## pipeline's QC plots) and MultiQC.
    citation = function(param = list()) {
      opts <- if (ezIsSpecified(param$cmdOptions)) param$cmdOptions else ""
      flagOn <- function(flag) {
        m <- regmatches(opts, regexec(paste0("(^|\\s)--", flag, "(\\s+|=|$)([^-\\s]\\S*)?"), opts, perl = TRUE))[[1]]
        length(m) > 0 && !tolower(m[4]) %in% "false"
      }
      contam <- flagOn("filter_contamination")
      blatDb <- contam && any(vapply(c("cdna", "ncrna", "pirna", "other_contamination"), flagOn, logical(1)))
      bowtie2Db <- blatDb || (contam && any(vapply(c("rrna", "trna"), flagOn, logical(1))))
      umi <- flagOn("with_umi")
      c(
        "Ewels, P.A. et al. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol 38(3), 276-278 (2020). https://doi.org/10.1038/s41587-020-0439-x",
        "Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. Nat Biotechnol 35(4), 316-319 (2017). https://doi.org/10.1038/nbt.3820",
        "nf-core/smrnaseq. nf-core. https://nf-co.re/smrnaseq",
        if (!flagOn("skip_fastp")) "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560",
        "Langmead, B., Trapnell, C., Pop, M. & Salzberg, S.L. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10, R25 (2009). https://doi.org/10.1186/gb-2009-10-3-r25",
        "Desvignes, T. et al. Unification of miRNA and isomiR research: the mirGFF3 format and the mirtop API. Bioinformatics 36(3), 698-703 (2020). https://doi.org/10.1093/bioinformatics/btz675",
        if (!flagOn("mirgenedb")) "Kozomara, A., Birgaoanu, M. & Griffiths-Jones, S. miRBase: from microRNA sequences to function. Nucleic Acids Research 47(D1), D155-D162 (2019). https://doi.org/10.1093/nar/gky1141",
        if (flagOn("mirgenedb")) "Clarke, A.W. et al. MirGeneDB 3.0: improved taxonomic sampling, uniform nomenclature of novel conserved microRNA families and updated covariance models. Nucleic Acids Research 53(D1), D116-D128 (2025). https://doi.org/10.1093/nar/gkae1094",
        if (umi) "Liu, D. Algorithms for efficiently collapsing reads with Unique Molecular Identifiers. PeerJ 7, e8275 (2019). https://doi.org/10.7717/peerj.8275",
        if (umi) "Smith, T., Heger, A. & Sudbery, I. UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research 27(3), 491-499 (2017). https://doi.org/10.1101/gr.209601.116",
        if (bowtie2Db) "Langmead, B. & Salzberg, S.L. Fast gapped-read alignment with Bowtie 2. Nature Methods 9(4), 357-359 (2012). https://doi.org/10.1038/nmeth.1923",
        if (blatDb) "Kent, W.J. BLAT—the BLAST-like alignment tool. Genome Research 12(4), 656-664 (2002). https://doi.org/10.1101/gr.229202",
        if (!flagOn("skip_mirdeep")) "Friedländer, M.R., Mackowiak, S.D., Li, N., Chen, W. & Rajewsky, N. miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades. Nucleic Acids Research 40(1), 37-52 (2012). https://doi.org/10.1093/nar/gkr688"
      )
    },
    initialize = function() {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreSmRnaSeq
      name <<- "EzAppNfCoreSmRnaSeq"
    }
  )
)

##' @description get an nf-core/smrnaseq-formatted csv file
getSmRnaSeqSampleSheet <- function(input) {
  nfSampleInfo <- ezFrame(
    sample = names(input$getFullPaths("Read1")),
    fastq_1 = input$getFullPaths("Read1")
  )

  return(nfSampleInfo)
}

##' @description split counts by sample
writePerSampleSmRnaCountFiles <- function(nfSampleInfo, countDir = ".") {
  sampleNames <- nfSampleInfo$sample
  sampleCountFiles <- paste0(countDir, "/", sampleNames, ".txt")
  annoColumnNames <- c("miRNA")
  x <- data.table::fread(file.path(countDir, "mirna.tsv"), data.table = FALSE)
  for (i in 1:nrow(nfSampleInfo)) {
    xSel <- x[, c(annoColumnNames, sampleNames[i])] |>
      dplyr::rename(
        "matchCounts" := !!sampleNames[i],
        "Identifier" := !!annoColumnNames
      )
    ezWrite.table(xSel, file = sampleCountFiles[i], row.names = FALSE)
  }
}
