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
  param[['gatk']] = file.path(Sys.getenv("GATK"), 'gatk')
  param[['javaMem']] = paste0(max(4, floor(as.numeric(param$ram) * 0.8)), "g")

  dataset = input$meta
  dataset[['GVCF [File]']] = input$getFullPaths("GVCF")
  ## read.table turns a factor level like "T" into logical TRUE; use the
  ## grouping only to split the samples and keep the names file-system safe
  if (ezIsSpecified(param$grouping)) {
    groups = as.character(input$getColumn(param$grouping))
  } else {
    groups = rep("allSamples", nrow(dataset))
  }
  groups = gsub("[^[:alnum:]_.-]", "_", groups)
  datasetCaseList = split(dataset, groups)
  results <- ezMclapply(
    names(datasetCaseList),
    runGatkPipelineRNASeq,
    param = param,
    datasetCaseList = datasetCaseList,
    mc.cores = param$cores
  )
  if (length(results) > 1) {
    jointVcf <- runGatkPipelineRNASeq(
      'allSamples',
      param = param,
      datasetCaseList = list(allSamples = dataset)
    )
  } else {
    jointVcf <- results[[1]]
  }

  ## the joint VCF is the deliverable; name it after the dataset, not the group
  vcfOutputFile = paste0(param$name, ".vcf.gz")
  varSummary = postProcessRnaVcf(jointVcf, vcfOutputFile, param)
  chromSizes = ezChromSizesFromVcf(vcfOutputFile)

  gt = getReportGenotypes(vcfOutputFile, param, input$getNames())
  varSummary$counts = rbind(
    varSummary$counts,
    data.frame(Step = "Biallelic SNPs used for the report", Records = nrow(gt))
  )

  makeRmdReport(
    input = input,
    output = output,
    param = param,
    chromSizes = chromSizes,
    gt = gt,
    varSummary = varSummary,
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
    combinedGVCF = paste0(caseName, ".g.vcf.gz")
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

  ## Joint genotyping produces a final (non-gVCF) VCF, hence .vcf (not .g.vcf);
  ## GATK bgzips and tabix-indexes it itself because of the .vcf.gz suffix
  vcfFile = paste0(caseName, '.vcf.gz')
  tmpVcf = paste0(caseName, '_temp.vcf.gz')
  cmd = paste(
    gatk,
    'GenotypeGVCFs',
    "-R",
    param$genomeSeq,
    fileCmd,
    "--output",
    tmpVcf
  )
  if (ezIsSpecified(param$specialOptions)) {
    cmd = paste(cmd, param$specialOptions)
  }
  ezSystem(paste(cmd, '2>>', myLog))
  ezSystem(paste('mv', tmpVcf, vcfFile))
  ezSystem(paste('mv', paste0(tmpVcf, ".tbi"), paste0(vcfFile, ".tbi")))
  ## drop the large intermediate combined gVCF
  if (!is.null(combinedGVCF)) {
    ezSystem(paste("rm -f", combinedGVCF, paste0(combinedGVCF, ".tbi")))
  }
  return(vcfFile)
}

##' @title Filter, normalise and annotate a joint-genotyped RNA-seq VCF
##' @description Applies, each optional, the GATK RNA-seq hard filters
##' (VariantFiltration: FS > 30, QD < 2, clusters of 3 SNVs within 35 bp),
##' splitting of multiallelic records with left-alignment (bcftools norm) and
##' SnpEff effect/impact annotation. Writes \code{outFile} (all records, FILTER
##' column set) and, if hard filtering ran, a PASS-only subset next to it.
##' @param vcf the joint-genotyped, bgzipped and indexed VCF; it is consumed.
##' @param outFile the name of the final VCF (.vcf.gz).
##' @param param the app parameters.
##' @return a list with a \code{counts} table and, if SnpEff ran, the
##' \code{impact} table and the name of the SnpEff summary html.
postProcessRnaVcf <- function(vcf, outFile, param) {
  prefix = sub("\\.vcf\\.gz$", "", outFile)
  myLog = paste0('log_', prefix, '_postprocessing.txt')
  res = list(
    counts = data.frame(
      Step = "Joint genotyped (GenotypeGVCFs)",
      Records = countVcfRecords(vcf)
    )
  )
  current = vcf
  intermediates = character()
  nextFile = function(tag) {
    f = paste0(prefix, ".", tag, ".vcf.gz")
    intermediates <<- c(intermediates, f)
    f
  }

  if (param$hardFilter) {
    filtered = nextFile("filtered")
    cmd = paste(
      param$gatk,
      "--java-options",
      paste0('"-Xmx', param$javaMem, '"'),
      "VariantFiltration",
      "-R",
      param$genomeSeq,
      "-V",
      current,
      "--window 35 --cluster 3",
      '--filter-name FS --filter-expression "FS > 30.0"',
      '--filter-name QD --filter-expression "QD < 2.0"',
      "-O",
      filtered
    )
    ezSystem(paste(cmd, "2>>", myLog))
    current = filtered
    res$filterCriteria = c(
      "FS > 30.0 (Fisher strand bias)",
      "QD < 2.0 (quality by depth)",
      "SnpCluster: 3 or more variants within 35 bp"
    )
    res$counts = rbind(
      res$counts,
      data.frame(
        Step = "PASS after VariantFiltration",
        Records = countVcfRecords(current, passOnly = TRUE)
      )
    )
  }

  if (param$normalize) {
    normed = nextFile("norm")
    splitBcf = paste0(prefix, ".split.bcf")
    ezSystem(paste(
      "bcftools norm -m -any --check-ref w",
      "-f",
      param$genomeSeq,
      "--threads",
      param$cores,
      "-Ob -o",
      splitBcf,
      current,
      "2>>",
      myLog
    ))
    ## splitting leaves spanning-deletion (*) alleles as records of their own;
    ## they carry no variant and SnpEff cannot annotate them
    ezSystem(paste(
      "bcftools view -e 'ALT=\"*\"' -Oz -o",
      normed,
      splitBcf
    ))
    ezSystem(paste("rm -f", splitBcf))
    ezSystem(paste("tabix -f -p vcf", normed))
    current = normed
    res$counts = rbind(
      res$counts,
      data.frame(
        Step = "After splitting multiallelic records (bcftools norm)",
        Records = countVcfRecords(current)
      )
    )
  }

  if (param$snpEff) {
    db = prepareSnpEffDb(param)
    annotated = nextFile("ann")
    res$snpEffHtml = paste0(prefix, ".snpEff.html")
    snpEffCsv = paste0(prefix, ".snpEff.csv")
    ezSystem(paste(
      "java",
      paste0("-Xmx", param$javaMem),
      "-Djava.io.tmpdir=.",
      "-jar $SnpEff/snpEff.jar ann -noLog",
      "-c",
      db$config,
      "-dataDir",
      db$dataDir,
      "-s",
      res$snpEffHtml,
      "-csvStats",
      snpEffCsv,
      db$name,
      current,
      "2>>",
      myLog,
      "| bgzip -c >",
      annotated
    ))
    ezSystem(paste("tabix -f -p vcf", annotated))
    current = annotated
    res$snpEffDb = db$name
    res$impact = countSnpEffImpact(current)
  }

  ezSystem(paste("mv", current, outFile))
  ezSystem(paste("mv", paste0(current, ".tbi"), paste0(outFile, ".tbi")))
  toRemove = setdiff(c(vcf, intermediates), current)
  ezSystem(paste(
    "rm -f",
    paste(toRemove, paste0(toRemove, ".tbi"), collapse = " ")
  ))

  if (param$hardFilter) {
    passFile = paste0(prefix, ".PASS.vcf.gz")
    ezSystem(paste("bcftools view -f PASS -Oz -o", passFile, outFile))
    ezSystem(paste("tabix -f -p vcf", passFile))
    res$passFile = passFile
  }
  res
}

##' @title Build or locate a SnpEff database for the reference
##' @description SnpEff ships no databases matching our Ensembl/GENCODE
##' builds, so one is built from the reference GTF and genome and cached in
##' \code{<annotation>/Genes/snpEff/<db>}. The build runs in a temporary
##' directory that is renamed into place, so concurrent jobs never see a
##' half-built database; the loser of a race discards its copy. If the
##' reference is not writable, the database is built in the working directory.
##' @param param the app parameters with \code{ezRef} and \code{javaMem}.
##' @return a list with the database \code{name}, \code{config} and \code{dataDir}.
prepareSnpEffDb <- function(param) {
  gtfFile = param$ezRef@refFeatureFile
  fastaFile = param$ezRef@refFastaFile
  parts = strsplit(param$ezRef@refBuild, "/")[[1]]
  parts = parts[-(1:2)]
  parts = parts[parts != "Annotation"]
  gtfTag = sub("\\.gtf$", "", basename(gtfFile))
  if (gtfTag != "genes") {
    parts = c(parts, gtfTag)
  }
  dbName = gsub("[^[:alnum:]_.-]", "_", paste(parts, collapse = "_"))

  cacheDir = file.path(dirname(gtfFile), "snpEff")
  if (!dir.exists(cacheDir)) {
    dir.create(cacheDir, showWarnings = FALSE)
  }
  if (file.access(cacheDir, 2) != 0) {
    warning("cannot write to ", cacheDir, "; building SnpEff DB locally")
    cacheDir = file.path(getwd(), "snpEff")
    dir.create(cacheDir, showWarnings = FALSE)
  }
  dbDir = file.path(cacheDir, dbName)
  db = list(
    name = dbName,
    config = file.path(dbDir, "snpEff.config"),
    dataDir = file.path(dbDir, "data")
  )
  if (file.exists(file.path(db$dataDir, dbName, "snpEffectPredictor.bin"))) {
    return(db)
  }

  tmpDir = file.path(
    cacheDir,
    paste0(dbName, ".tmp.", Sys.info()[["nodename"]], ".", Sys.getpid())
  )
  genomeDir = file.path(tmpDir, "data", dbName)
  dir.create(genomeDir, recursive = TRUE)
  file.symlink(normalizePath(gtfFile), file.path(genomeDir, "genes.gtf"))
  file.symlink(normalizePath(fastaFile), file.path(genomeDir, "sequences.fa"))
  config = c(
    "data.dir = ./data/",
    paste0(dbName, ".genome : ", parts[1])
  )
  ## SnpEff's default codon table is the standard one; mitochondria need theirs
  chroms = ezRead.table(
    paste0(fastaFile, ".fai"),
    header = FALSE,
    row.names = NULL
  )[[1]]
  for (mt in intersect(chroms, c("MT", "M", "chrM", "chrMT"))) {
    config = c(
      config,
      paste0(dbName, ".", mt, ".codonTable : Vertebrate_Mitochondrial")
    )
  }
  writeLines(config, file.path(tmpDir, "snpEff.config"))
  ezSystem(paste(
    "java",
    paste0("-Xmx", param$javaMem),
    "-Djava.io.tmpdir=.",
    "-jar $SnpEff/snpEff.jar build -gtf22 -noCheckCds -noCheckProtein -noLog",
    "-c",
    file.path(tmpDir, "snpEff.config"),
    "-dataDir",
    file.path(tmpDir, "data"),
    "-v",
    dbName,
    "2>",
    paste0("log_snpEff_build_", dbName, ".txt")
  ))
  if (!file.rename(tmpDir, dbDir)) {
    unlink(tmpDir, recursive = TRUE)
    if (!file.exists(file.path(db$dataDir, dbName, "snpEffectPredictor.bin"))) {
      stop("could not install the SnpEff database in ", dbDir)
    }
  }
  db
}

countVcfRecords <- function(vcf, passOnly = FALSE) {
  cmd = paste("bcftools view -H -G", if (passOnly) "-f PASS", vcf, "| wc -l")
  as.integer(ezSystem(cmd, intern = TRUE))
}

## SnpEff sorts the annotations of a record by impact, so the first one is
## the most severe; count records, not transcript effects as SnpEff's own stats
countSnpEffImpact <- function(vcf) {
  x = data.table::fread(
    cmd = paste(
      "bcftools query -f '%FILTER\\t%INFO/ANN\\n'",
      vcf,
      "| cut -d, -f1 | cut -d'|' -f1,3"
    ),
    header = FALSE,
    sep = "\t",
    quote = "",
    col.names = c("FILTER", "ANN")
  )
  impact = ifelse(grepl("|", x$ANN, fixed = TRUE), sub(".*\\|", "", x$ANN), "none")
  impact = factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER", "none"))
  res = data.frame(
    Impact = levels(impact),
    `All records` = as.integer(table(impact)),
    check.names = FALSE
  )
  if (any(x$FILTER != ".")) {
    res$PASS = as.integer(table(impact[x$FILTER == "PASS"]))
  }
  res
}

##' @title Genotype matrix of biallelic SNPs for the report
##' @description Uses the PASS (if hard filtering ran) biallelic SNPs of the
##' final VCF, randomly subsampled to at most 500k for the report plots. Calls
##' with a depth below \code{param$minReadDepth} are set to "lowCov".
##' @return a character matrix of genotypes (variants x samples); it has zero
##' rows if there is no SNP.
getReportGenotypes <- function(vcf, param, sampleNames, targetN = 500000L) {
  snpVcf = "tmp.snps.vcf.gz"
  ezSystem(paste(
    "bcftools view -m2 -M2 -v snps",
    if (param$hardFilter) "-f PASS",
    "-Oz -o",
    snpVcf,
    vcf
  ))
  ezSystem(paste("tabix -f -p vcf", snpVcf))
  nTotal = countVcfRecords(snpVcf)
  message("Biallelic SNPs: ", nTotal)
  if (nTotal == 0) {
    ezSystem(paste("rm -f", snpVcf, paste0(snpVcf, ".tbi")))
    return(matrix(
      character(0),
      nrow = 0,
      ncol = length(sampleNames),
      dimnames = list(NULL, sampleNames)
    ))
  }

  if (nTotal > targetN) {
    pos = data.table::fread(
      cmd = paste("bcftools query -f '%CHROM\\t%POS\\n'", snpVcf),
      header = FALSE,
      col.names = c("CHROM", "POS"),
      colClasses = c("character", "integer")
    )
    set.seed(1L)
    pos = pos[sort(sample.int(nTotal, targetN)), ]
    regionsFile = "tmp.snps.regions.txt"
    data.table::fwrite(pos, regionsFile, sep = "\t", col.names = FALSE)
    subVcf = "tmp.snps.sub.vcf.gz"
    ezSystem(paste(
      "bcftools view -T",
      regionsFile,
      "-Oz -o",
      subVcf,
      snpVcf
    ))
    ezSystem(paste("mv", subVcf, snpVcf))
    ezSystem(paste("tabix -f -p vcf", snpVcf))
    unlink(regionsFile)
  }

  genotype = geno(readVcf(snpVcf, genome = "genomeDummy"))
  gt = genotype$GT
  ## those calls will become NA in subsequent analyses
  gt[!is.na(genotype$DP) & genotype$DP < param$minReadDepth] = "lowCov"
  ezSystem(paste("rm -f", snpVcf, paste0(snpVcf, ".tbi")))
  gt
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
            Description = "genotypes with a lower read depth are shown as low coverage in the report"
          ),
          hardFilter = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "GATK RNA-seq hard filters: FS > 30, QD < 2, SNP clusters (3 in 35 bp); also writes a PASS-only VCF"
          ),
          normalize = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "split multiallelic records and left-align indels with bcftools norm"
          ),
          snpEff = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "annotate variant effects and impact with SnpEff; the database is built from the reference annotation on first use"
          )
        )
      }
    )
  )
