###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodJoinGenoTypes = function(input = NA, output = NA, param = NA) {
  setwdNew(param[['name']])
  param[['genomeSeq']] = param$ezRef["refFastaFile"]
  param[['species']] = limma::strsplit2(param$ezRef['refBuild'], '/')[1]
  param[['knownSites']] = list.files(
    param$ezRef["refVariantsDir"],
    pattern = 'vcf.gz$',
    full.names = T
  )
  param[['dbsnpFile']] = param$knownSites[grep(
    'dbsnp.*vcf.gz$',
    param$knownSites
  )]
  param[['javaCall']] = paste("java", "-Djava.io.tmpdir=.")
  param[['gatk']] = file.path(Sys.getenv("GATK"), 'gatk')

  dataset = input$meta
  dataset[['GVCF [File]']] = input$getFullPaths("GVCF")
  datasetCaseList = split(dataset, input$getColumn(param$grouping))
  results <- ezMclapply(
    names(datasetCaseList),
    runGatkPipeline,
    param = param,
    datasetCaseList = datasetCaseList,
    mc.cores = param$cores
  )
  return("Success")
}

runGatkPipeline = function(caseName, param = NA, datasetCaseList = NULL) {
  gatk = param[['gatk']]
  datasetCase <- datasetCaseList[[caseName]]
  myLog = paste0('log_', caseName, '.txt')
  if (param$targetFile != '') {
    param$targetFile = file.path(TARGET_ENRICHMENT_DESIGN_DIR, param$targetFile)
    targetOption <- paste("-L", param$targetFile)
  } else {
    targetOption <- ""
  }

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

  ## Joint genotyping (and downstream recalibration/annotation) produces a final
  ## (non-gVCF) VCF, hence .vcf (not .g.vcf)
  GenotypeGVCF = paste(gatk, 'GenotypeGVCFs')
  vcfFile = paste0(caseName, '.vcf')
  tmpGvcf = paste0(caseName, '_temp.vcf')
  cmd = paste(
    GenotypeGVCF,
    "-R",
    param$genomeSeq,
    fileCmd,
    "--output",
    tmpGvcf,
    targetOption
  )
  if (length(param$dbsnpFile) == 1) {
    cmd <- paste(cmd, "--dbsnp", param$dbsnpFile)
  }
  ezSystem(paste(cmd, '2>', myLog))
  ezSystem(paste('mv', tmpGvcf, vcfFile))
  ezSystem(paste('mv', paste0(tmpGvcf, ".idx"), paste0(vcfFile, ".idx")))
  ## drop the large intermediate combined gVCF (no longer overwritten now that names differ)
  if (!is.null(combinedGVCF)) {
    ezSystem(paste("rm -f", combinedGVCF, paste0(combinedGVCF, ".idx")))
  }

  if (param$species == 'Homo_sapiens') {
    hapmapFile = param$knownSites[grep('hapmap_.*vcf.gz$', param$knownSites)]
    h1000G_omniFile = param$knownSites[grep(
      '1000G_omni.*vcf.gz$',
      param$knownSites
    )]
    h1000G_phase1File = param$knownSites[grep(
      '1000G_phase1.snps.*vcf.gz$',
      param$knownSites
    )]

    if (param$recalibrateVariants) {
      cmd = paste(
        gatk,
        'VariantRecalibrator',
        "-R",
        param$genomeSeq,
        "--variant",
        vcfFile,
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0",
        hapmapFile,
        "--resource:omni,known=false,training=true,truth=false,prior=12.0",
        h1000G_omniFile,
        "--resource:1000G,known=false,training=true,truth=false,prior=10.0",
        h1000G_phase1File,
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0",
        param$dbsnpFile,
        "-an QD -an FS -an ReadPosRankSum -an DP",
        "-mode SNP",
        "--output",
        paste0(caseName, '_raw.SNPs.recal'),
        "--tranches-file",
        paste0(caseName, '_raw.SNPs.tranches'),
        "--max-attempts 3 --max-gaussians 4",
        targetOption
      )
      ezSystem(paste(cmd, '2>>', myLog))

      #Apply RecalibrationOutput for SNPs:
      cmd = paste(
        gatk,
        'ApplyVQSR',
        "-R",
        param$genomeSeq,
        "--variant",
        vcfFile,
        "-mode SNP",
        "--recal-file",
        paste0(caseName, '_raw.SNPs.recal'),
        "--tranches-file",
        paste0(caseName, '_raw.SNPs.tranches'),
        "--output",
        tmpGvcf,
        "--truth-sensitivity-filter-level 99.0",
        targetOption
      )
      ezSystem(paste(cmd, '2>>', myLog))
      ezSystem(paste('mv', tmpGvcf, vcfFile))
      ezSystem(paste('mv', paste0(tmpGvcf, ".idx"), paste0(vcfFile, ".idx")))
      ezSystem(paste0('rm -f ', caseName, "*.tranches ", caseName, "*recal*"))
    }

    #2.Run VariantRecalibration for InDels:
    if (param$recalibrateInDels) {
      millsFile = param$knownSites[grep('Mills.*vcf.gz$', param$knownSites)]
      cmd = paste(
        gatk,
        'VariantRecalibrator',
        "-R",
        param$genomeSeq,
        "--variant",
        vcfFile,
        "-resource:mills,known=false,training=true,truth=true,prior=12.0",
        millsFile,
        "--max-gaussians 4 --minimum-bad-variants 500 --max-attempts 3", #might be suboptimal
        "-an QD -an FS -an ReadPosRankSum -an DP",
        "-mode INDEL",
        "--output",
        paste0(caseName, "_raw.InDels.recal"),
        "--tranches-file",
        paste0(caseName, "_raw.InDels.tranches"),
        targetOption
      )
      ezSystem(paste(cmd, '2>>', myLog))

      #Apply RecalibrationOutput for InDels:
      cmd = paste(
        gatk,
        'ApplyVQSR',
        "-R",
        param$genomeSeq,
        "--variant",
        vcfFile,
        "-mode INDEL",
        "--recal-file",
        paste0(caseName, '_raw.InDels.recal'),
        "--tranches-file",
        paste0(caseName, '_raw.InDels.tranches'),
        "--output",
        tmpGvcf,
        "--truth-sensitivity-filter-level 99.0",
        targetOption
      )
      ezSystem(paste(cmd, '2>>', myLog))
      ezSystem(paste('mv', tmpGvcf, vcfFile))
      ezSystem(paste('mv', paste0(tmpGvcf, ".idx"), paste0(vcfFile, ".idx")))
      ezSystem(paste0('rm -f ', caseName, "*.tranches ", caseName, "*recal*"))
    }
    #Add ExAc-Annotation:
    # runExAC = FALSE
    # if(runExAC){
    #     ExAcFile = param$knownSites[grep('ExAC.*vcf.gz$', param$knownSites)]
    #     VariantAnnotation = paste(gatk,'VariantAnnotator')
    #     cmd = paste(VariantAnnotation, "-R", param$genomeSeq,
    #                 "--variant", vcfFile,
    #                 "--output ",paste0(vcfFile, "_annotated.vcf"),
    #                 "--resource:ExAC",  ExAcFile,
    #                 "--expression ExAC.AF",
    #                 "--expression ExAC.AC",
    #                 "--expression ExAC.AN",
    #                 "--expression ExAC.AC_Het",
    #                 "--expression ExAC.AC_Hom")
    #     if(param$targetFile != ''){
    #         cmd = paste(cmd, "-L", param$targetFile) }
    #     ezSystem(paste(cmd,'2>>',myLog))
    # } else {
    #     ezSystem(paste('mv', vcfFile, paste0(vcfFile, "_annotated.vcf")))
    # }
    dbNsfpCheck <- dbNsfpBuildCheck(param$refBuild, param$dbNSFP_file)
    if (!dbNsfpCheck$use) {
      ezLog("Skipping SnpSift dbnsfp annotation: ", dbNsfpCheck$reason, level = "warn")
    }
  }
  if (param$species == 'Homo_sapiens' && dbNsfpCheck$use) {
    param$dbNSFP_file = file.path('/srv/GT/databases/dbNSFP', param$dbNSFP_file)
    cmd = paste(
      param$javaCall,
      "-jar",
      "$SnpEff/SnpSift.jar",
      "dbnsfp -f",
      param$dbNSFP_fields,
      "-v -db",
      param$dbNSFP_file,
      vcfFile,
      ">",
      tmpGvcf
    )
    #for hg38 upgrade - dbsfp starting from v4: gnomAD_exomes_AF,M-CAP_score,M-CAP_rankscore,M-CAP_pred,VindijiaNeandertal
    ezSystem(paste(cmd, '2>>', myLog))
    ezSystem(paste("mv", tmpGvcf, vcfFile))
    ezSystem(paste("gatk", "IndexFeatureFile -I", vcfFile))
  }

  #SnpEff:
  if (param$snpEffDB != '') {
    ## SnpEff databases use Ensembl contig names: strip a leading "chr" from the
    ## contig names only (CHROM column and ##contig header lines), not from INFO.
    headerLines <- system(paste("bcftools view -h", vcfFile), intern = TRUE)
    renameMap <- vcfChrRenameMap(headerLines)
    if (nrow(renameMap) > 0) {
      mapFile <- paste0(caseName, "_chrRename.txt")
      write.table(
        renameMap,
        mapFile,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
      )
      renamedVcf <- paste0(caseName, "_noChr.vcf")
      renameScript <- paste0(caseName, "_renameChrs.sh")
      writeLines(
        c(
          "set -euo pipefail",
          paste(
            "bcftools annotate --rename-chrs",
            mapFile,
            "-Ov -o",
            renamedVcf,
            vcfFile
          ),
          paste("mv", renamedVcf, vcfFile),
          paste("rm -f", paste0(vcfFile, ".idx"))
        ),
        renameScript
      )
      ezSystem(paste("bash", renameScript, "2>>", myLog))
      ezSystem(paste("gatk", "IndexFeatureFile -I", vcfFile, "2>>", myLog))
    }
    htmlOutputFile = paste0(caseName, '.html')
    protCodingTranscripts <- character(0)
    if (param$proteinCodingTranscriptsOnly) {
      gtfFile = param$ezRef@refFeatureFile
      gtf = rtracklayer::import(gtfFile)
      idx = gtf$type == 'transcript'
      gtf = data.frame(gtf[idx])
      biotypeCol <- intersect(
        c("transcript_type", "transcript_biotype"),
        colnames(gtf)
      )
      if (length(biotypeCol) > 0) {
        biotypes <- gtf[[biotypeCol[1]]]
        protCodingTranscripts <- unique(gtf$transcript_id[
          !is.na(biotypes) & biotypes == 'protein_coding'
        ])
      }
      if (length(protCodingTranscripts) == 0) {
        ezLog(
          "No protein-coding transcripts found in ",
          gtfFile,
          " (no transcript_type/transcript_biotype attribute); ",
          "running SnpEff on all transcripts",
          level = "warn"
        )
      }
    }
    if (length(protCodingTranscripts) > 0) {
      write.table(
        protCodingTranscripts,
        'protCodingTranscripts.txt',
        col.names = F,
        row.names = F,
        quote = F
      )
      cmd = paste(
        param$javaCall,
        "-jar",
        "$SnpEff/snpEff.jar",
        "-onlyTr protCodingTranscripts.txt",
        "-stats",
        htmlOutputFile,
        param$snpEffDB,
        "-v",
        vcfFile,
        ">",
        tmpGvcf
      )
    } else {
      cmd = paste(
        param$javaCall,
        "-jar",
        "$SnpEff/snpEff.jar",
        "-s",
        htmlOutputFile,
        param$snpEffDB,
        "-v",
        vcfFile,
        ">",
        tmpGvcf
      )
    }
    ezSystem(paste(cmd, '2>>', myLog))
    ezSystem(paste("mv", tmpGvcf, vcfFile))
    ezSystem(paste("gatk", "IndexFeatureFile -I", vcfFile))
  }

  ezSystem(paste("bgzip", vcfFile))
  ezSystem(paste0("tabix -p vcf ", vcfFile, ".gz"))
  return(vcfFile)
}


##' @template app-template
##' @templateVar method ezMethodJoinGenoTypes(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppJoinGenoTypes <-
  setRefClass(
    "EzAppJoinGenoTypes",
    contains = "EzApp",
    methods = list(
      ## GATK CombineGVCFs/GenotypeGVCFs + bgzip/tabix unconditional; SnpSift dbnsfp + dbNSFP only
      ## for Homo_sapiens with a dbNSFP file whose build matches refBuild (VQSR, also human-only,
      ## has no separate reference); SnpEff only when snpEffDB is set.
      citation = function(param = list()) {
        human <- isTRUE(grepl("^Homo_sapiens/", param$refBuild))
        useDbNsfp <- human &&
          dbNsfpBuildCheck(param$refBuild, param$dbNSFP_file)$use
        c(
          "McKenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research 20, 1297-1303 (2010). https://doi.org/10.1101/gr.107524.110",
          "Poplin, R. et al. Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv (2017) [preprint, not peer-reviewed]. https://doi.org/10.1101/201178",
          "DePristo, M.A. et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics 43, 491-498 (2011). https://doi.org/10.1038/ng.806",
          if (useDbNsfp) "Cingolani, P. et al. Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift. Frontiers in Genetics 3, 35 (2012). https://doi.org/10.3389/fgene.2012.00035",
          if (useDbNsfp) "Liu, X. et al. dbNSFP v4: a comprehensive database of transcript-specific functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome Medicine 12, 103 (2020). https://doi.org/10.1186/s13073-020-00803-9",
          if (ezIsSpecified(param$snpEffDB)) "Cingolani, P. et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. Fly 6, 80-92 (2012). https://doi.org/10.4161/fly.19695",
          "Bonfield, J.K. et al. HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience 10, giab007 (2021). https://doi.org/10.1093/gigascience/giab007",
          "Li, H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27, 718-719 (2011). https://doi.org/10.1093/bioinformatics/btq671"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodJoinGenoTypes
        name <<- "EzAppJoinGenoTypes"
        appDefaults <<- rbind(
          targetFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "restrict to targeted genomic regions"
          ),
          vcfFilt.minAltCount = ezFrame(
            Type = "integer",
            DefaultValue = 3,
            Description = "minimum coverage for the alternative variant"
          ),
          vcfFilt.minReadDepth = ezFrame(
            Type = "integer",
            DefaultValue = 20,
            Description = "minimum read depth"
          ),
          recalibrateVariants = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "recalibrateVariants"
          ),
          recalibrateInDels = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "recalibrate InDels"
          ),
          snpEffDB = ezFrame(
            Type = "character",
            DefaultValue = '',
            Description = "SnpEff DB name"
          ),
          proteinCodingTranscriptsOnly = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "annotate variants only with protCod variants"
          ),
          dbNSFP_file = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "name of dbNSFP file under /srv/GT/databases/dbNSFP"
          ),
          dbNSFP_fields = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "info fields for vcf annotation"
          )
        )
      }
    )
  )

##' @title Genome build of a reference or dbNSFP file name
##' @param x a refBuild string or a file name.
##' @return "GRCh38", "GRCh37" or NA if no build can be recognised.
genomeBuildFromName <- function(x) {
  if (!ezIsSpecified(x)) {
    return(NA_character_)
  }
  if (grepl("GRCh38|hg38", x, ignore.case = TRUE)) {
    return("GRCh38")
  }
  if (grepl("GRCh37|hg19|b37", x, ignore.case = TRUE)) {
    return("GRCh37")
  }
  NA_character_
}

##' @title Check whether a dbNSFP file can annotate a VCF of a reference build
##' @param refBuild the refBuild parameter.
##' @param dbNsfpFile the dbNSFP file name (relative to /srv/GT/databases/dbNSFP).
##' @return list with \code{use} (logical) and \code{reason} (character).
dbNsfpBuildCheck <- function(refBuild, dbNsfpFile) {
  if (!ezIsSpecified(dbNsfpFile)) {
    return(list(use = FALSE, reason = "no dbNSFP_file given"))
  }
  refGenomeBuild <- genomeBuildFromName(refBuild)
  dbBuild <- genomeBuildFromName(basename(dbNsfpFile))
  if (!is.na(refGenomeBuild) && !is.na(dbBuild) && refGenomeBuild != dbBuild) {
    return(list(
      use = FALSE,
      reason = paste0(
        "dbNSFP file ", dbNsfpFile, " is for ", dbBuild,
        " but refBuild ", refBuild, " is ", refGenomeBuild
      )
    ))
  }
  list(use = TRUE, reason = "")
}

##' @title Contig rename map that strips a leading "chr"
##' @param headerLines VCF header lines.
##' @return data.frame with columns old and new, one row per ##contig whose ID
##' starts with "chr".
vcfChrRenameMap <- function(headerLines) {
  contigLines <- grep("^##contig=<ID=", headerLines, value = TRUE)
  ids <- sub("^##contig=<ID=([^,>]+).*$", "\\1", contigLines)
  ids <- ids[grepl("^chr", ids)]
  data.frame(old = ids, new = sub("^chr", "", ids), stringsAsFactors = FALSE)
}
