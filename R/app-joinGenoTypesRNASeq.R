###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodJoinGenoTypesRNASeq <- function(input = NA, output = NA, param = NA) {
  library(VariantAnnotation)

  setwdNew(param[['name']])
  param[['genomeSeq']] = param$ezRef["refFastaFile"]
  param[['species']] = limma::strsplit2(param$ezRef['refBuild'], '/')[1]
  param[['javaCall']] = paste("java", "-Djava.io.tmpdir=.")
  param[['gatk']] = file.path(Sys.getenv("GATK"), 'gatk')

  dataset = input$meta
  dataset[['GVCF [File]']] = input$getFullPaths("GVCF")
  datasetCaseList = split(dataset, input$getColumn(param$grouping))
  results <- ezMclapply(
    names(datasetCaseList),
    runGatkPipelineRNASeq,
    param = param,
    datasetCaseList = datasetCaseList,
    mc.cores = param$cores
  )
  if (length(results) > 1) {
    results <- runGatkPipelineRNASeq(
      'allSamples',
      param = param,
      datasetCaseList = list(allSamples = dataset)
    )
  }
  vcfOutputFile = results[[1]]
  chromSizes = ezChromSizesFromVcf(vcfOutputFile)

  system(sprintf(
    "bcftools view -m2 -M2 -v snps %s -Oz -o tmp.snps.vcf.gz",
    shQuote(vcfOutputFile)
  ))
  system('tabix -f -p vcf tmp.snps.vcf.gz')

  tmp_vcf <- "tmp.snps.vcf.gz"
  out_vcf <- "snps.sample500k.vcf.gz"
  target_n <- 500000L
  seed <- 1L

  # Extract CHROM and POS
  cmd_pos <- sprintf(
    "bcftools query -f '%%CHROM\\t%%POS\\n' %s",
    shQuote(tmp_vcf)
  )
  pos_txt <- system(cmd_pos, intern = TRUE)
  if (length(pos_txt) == 0L) {
    stop("bcftools query returned no lines")
  }

  pos <- read.table(
    text = pos_txt,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(pos) <- c("CHROM", "POS")

  n_total <- nrow(pos)
  message("Total records: ", n_total)
  if (target_n > n_total) {
    stop("target_n > n_total")
  }

  # Sample exact N
  set.seed(seed)
  idx <- sample.int(n_total, target_n, replace = FALSE)
  pos_sub <- pos[idx, , drop = FALSE]

  # Sort (good hygiene)
  pos_sub <- pos_sub[order(pos_sub$CHROM, pos_sub$POS), , drop = FALSE]

  # Write regions file: 2 columns (CHROM POS), 1-based
  regions_file <- tempfile(fileext = ".regions.txt")
  write.table(
    pos_sub,
    file = regions_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # Subset
  cmd_view <- sprintf(
    "bcftools view -R %s -Oz -o %s %s",
    shQuote(regions_file),
    shQuote(out_vcf),
    shQuote(tmp_vcf)
  )
  status <- system(cmd_view)
  if (status != 0) {
    stop("bcftools view failed (exit ", status, ")")
  }

  # Index output
  status <- system(sprintf("tabix -f -p vcf %s", shQuote(out_vcf)))
  if (status != 0) {
    stop("tabix failed (exit ", status, ")")
  }

  # Verify
  n_out <- as.numeric(trimws(system(
    sprintf("bcftools index -n %s", shQuote(out_vcf)),
    intern = TRUE
  )))
  message("Output records: ", n_out)

  vcfOutputFile = 'snps.sample500k.vcf.gz'

  genotype = geno(readVcf(vcfOutputFile, genome = "genomeDummy"))
  gt = genotype$GT
  gt[genotype$DP < param$minReadDepth] = "lowCov" ## those calls will become NA in subsequent analyses

  makeRmdReport(
    input = input,
    output = output,
    param = param,
    chromSizes = chromSizes,
    gt = gt,
    rmdFile = "Mpileup.Rmd",
    reportTitle = 'GATK Joint Genotyping RNA-Seq Report'
  )
  return("Success")
}

runGatkPipelineRNASeq <- function(
  caseName,
  param = NA,
  datasetCaseList = NULL
) {
  gatk = param[['gatk']]
  datasetCase <- datasetCaseList[[caseName]]
  myLog = paste0('log_', caseName, '.txt')

  ##Create CombinedGVCF: (this is a genuine multi-sample gVCF, so .g.vcf is correct)
  combinedGVCF = NULL
  if (nrow(datasetCase) > 1) {
    combinedGVCF = paste0(caseName, ".g.vcf")
    cmd = paste(
      paste(gatk, 'CombineGVCFs'),
      "-R",
      param$genomeSeq,
      paste('--variant', datasetCase[['GVCF [File]']], collapse = ' '),
      "--output",
      combinedGVCF
    )
    ezSystem(paste(cmd, '2>', myLog))
    fileCmd = paste("--variant", combinedGVCF)
  } else {
    fileCmd = paste("--variant", datasetCase[['GVCF [File]']])
  }

  ## Joint genotyping produces a final (non-gVCF) VCF, hence .vcf (not .g.vcf)
  GenotypeGVCF = paste(gatk, 'GenotypeGVCFs')
  vcfFile = paste0(caseName, '.vcf')
  tmpVcf = paste0(caseName, '_temp.vcf')
  cmd = paste(GenotypeGVCF, "-R", param$genomeSeq, fileCmd, "--output", tmpVcf)
  ezSystem(paste(cmd, '2>', myLog))
  ezSystem(paste('mv', tmpVcf, vcfFile))
  ezSystem(paste('mv', paste0(tmpVcf, ".idx"), paste0(vcfFile, ".idx")))
  ezSystem(paste("bgzip", vcfFile))
  ezSystem(paste0("tabix -p vcf ", vcfFile, ".gz"))
  ## drop the large intermediate combined gVCF (no longer overwritten now that names differ)
  if (!is.null(combinedGVCF)) {
    ezSystem(paste("rm -f", combinedGVCF, paste0(combinedGVCF, ".idx")))
  }
  return(paste0(vcfFile, '.gz'))
}


##' @template app-template
##' @templateVar method ezMethodJoinGenoTypesRNASeq(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppJoinGenoTypesRNASeq <-
  setRefClass(
    "EzAppJoinGenoTypesRNASeq",
    contains = "EzApp",
    methods = list(
      methods_facts = function(param = list()) {
        grp <- if (ezIsSpecified(param$grouping)) paste0("the ", param$grouping, " column") else "the grouping column"
        dp <- if (ezIsSpecified(param$minReadDepth)) paste0(param$minReadDepth, " (minReadDepth)") else "minReadDepth"
        nSamples <- length(methodsInput(param, "Name"))
        c(
          ## ezMethodJoinGenoTypesRNASeq / runGatkPipelineRNASeq (app-joinGenoTypesRNASeq.R:19-33, 134-176)
          paste0("The gVCFs were genotyped per group of ", grp, ": GATK CombineGVCFs merged the gVCFs of a group with more than one sample and GenotypeGVCFs made one VCF per group, with no --dbsnp, interval or other extra option; when there was more than one group, all samples were also genotyped together and only this all-samples VCF was used for the report."),
          "No variant filtering (hard filters or VQSR) and no variant effect annotation was done; the per-group and all-samples VCFs are unfiltered.",
          ## ezMethodJoinGenoTypesRNASeq (app-joinGenoTypesRNASeq.R:36-116)
          "For the report, the VCF was reduced to biallelic SNPs (bcftools view -m2 -M2 -v snps) and exactly 500,000 of them were drawn at random (R sample.int, seed 1); the job stops when fewer than 500,000 biallelic SNPs are present.",
          ## ezMethodJoinGenoTypesRNASeq (app-joinGenoTypesRNASeq.R:118-120)
          paste0("Genotype calls with read depth (DP) below ", dp, " were set to missing for the report."),
          ## Mpileup.Rmd clustering chunk, eval=length(samples) > 3
          if (nSamples == 0 || nSamples > 3) "With more than 3 samples, samples were clustered by hierarchical clustering (hclust, method ward.D2) of Euclidean distances between genotype codes (0/0 = -1, 0/1 = 0, 1/1 = 1; any other genotype, including low-depth and phased calls, missing); the dendrogram is drawn only when every distance could be computed.",
          if (nSamples > 0 && nSamples <= 3) "No sample clustering was done because there were 3 or fewer samples.",
          ## Mpileup.Rmd chrom variants chunk
          "Genotypes were plotted along every chromosome whose name has no '.' or '_' and that is longer than 20,000 bp."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodJoinGenoTypesRNASeq
        name <<- "JoinGenoTypesRNASeq"
        appDefaults <<- rbind(
          minReadDepth = ezFrame(
            Type = "integer",
            DefaultValue = "20",
            Description = "use for clustering only SNV with coverage higher than"
          )
        )
      }
    )
  )
