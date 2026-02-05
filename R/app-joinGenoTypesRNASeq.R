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
  
  system('bcftools view -m2 -M2 -v snps allSamples.g.vcf.gz -Oz -o tmp.snps.vcf.gz')
  system('tabix -f -p vcf tmp.snps.vcf.gz')
  
  tmp_vcf  <- "tmp.snps.vcf.gz"
  out_vcf  <- "snps.sample500k.vcf.gz"
  target_n <- 500000L
  seed     <- 1L
  
  # Extract CHROM and POS
  cmd_pos <- sprintf("bcftools query -f '%%CHROM\\t%%POS\\n' %s", shQuote(tmp_vcf))
  pos_txt <- system(cmd_pos, intern = TRUE)
  if (length(pos_txt) == 0L) stop("bcftools query returned no lines")
  
  pos <- read.table(text = pos_txt, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(pos) <- c("CHROM", "POS")
  
  n_total <- nrow(pos)
  message("Total records: ", n_total)
  if (target_n > n_total) stop("target_n > n_total")
  
  # Sample exact N
  set.seed(seed)
  idx <- sample.int(n_total, target_n, replace = FALSE)
  pos_sub <- pos[idx, , drop = FALSE]
  
  # Sort (good hygiene)
  pos_sub <- pos_sub[order(pos_sub$CHROM, pos_sub$POS), , drop = FALSE]
  
  # Write regions file: 2 columns (CHROM POS), 1-based
  regions_file <- tempfile(fileext = ".regions.txt")
  write.table(pos_sub, file = regions_file, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Subset
  cmd_view <- sprintf(
      "bcftools view -R %s -Oz -o %s %s",
      shQuote(regions_file), shQuote(out_vcf), shQuote(tmp_vcf)
  )
  status <- system(cmd_view)
  if (status != 0) stop("bcftools view failed (exit ", status, ")")
  
  # Index output
  status <- system(sprintf("tabix -f -p vcf %s", shQuote(out_vcf)))
  if (status != 0) stop("tabix failed (exit ", status, ")")
  
  # Verify
  n_out <- as.numeric(trimws(system(sprintf("bcftools index -n %s", shQuote(out_vcf)), intern = TRUE)))
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

  ##Create CombinedGVCF:
  if (nrow(datasetCase) > 1) {
    GenotypeGVCF = paste(gatk, 'CombineGVCFs')
    myGVCF = paste0(caseName, ".g.vcf")
    cmd = paste(
      GenotypeGVCF,
      "-R",
      param$genomeSeq,
      paste('--variant', datasetCase[['GVCF [File]']], collapse = ' '),
      "--output",
      myGVCF
    )
    ezSystem(paste(cmd, '2>', myLog))
    fileCmd = paste("--variant", myGVCF)
  } else {
    fileCmd = paste("--variant", datasetCase[['GVCF [File]']])
  }

  GenotypeGVCF = paste(gatk, 'GenotypeGVCFs')
  gvcfFile = paste0(caseName, '.g.vcf')
  tmpGvcf = paste0(caseName, '_temp.vcf')
  cmd = paste(GenotypeGVCF, "-R", param$genomeSeq, fileCmd, "--output", tmpGvcf)
  ezSystem(paste(cmd, '2>', myLog))
  ezSystem(paste('mv', tmpGvcf, gvcfFile))
  ezSystem(paste('mv', paste0(tmpGvcf, ".idx"), paste0(gvcfFile, ".idx")))
  ezSystem(paste("bgzip", gvcfFile))
  ezSystem(paste0("tabix -p vcf ", gvcfFile, ".gz"))
  return(paste0(gvcfFile, '.gz'))
}


##' @template app-template
##' @templateVar method ezMethodJoinGenoTypesRNASeq(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppJoinGenoTypesRNASeq <-
  setRefClass(
    "EzAppJoinGenoTypesRNASeq",
    contains = "EzApp",
    methods = list(
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
