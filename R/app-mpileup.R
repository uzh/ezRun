###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMpileup = function(input=NA, output=NA, param=NA){

  require("VariantAnnotation")
  
  reportDir = basename(output$getColumn("Report"))
  htmlFile = basename(output$getColumn("Html"))
  vcfOutputFile = output$getColumn("VCF")
  
  bamFiles = input$getFullPaths("BAM")
  bamDataset = input$meta
  genomeSeq = param$ezRef["refFastaFile"]
  nBamsInParallel = min(4, param$cores)
  bamFilesClean = ezMclapply(names(bamFiles), function(sampleName){
    javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", floor(param$ram/nBamsInParallel), "g")
    setwdNew(paste(sampleName, "proc", sep="-"))
    bf = bamFiles[sampleName]
    obf = file.path(getwd(), basename(bf))
    if (!is.null(param$region) && param$region != ""){
      bamParam = ScanBamParam(what=scanBamWhat())
      reg = splitRegion(param$region)
      if (is.na(reg$end)){
        seqLengths = ezBamSeqLengths(bf)
        reg$start = 1
        reg$end = seqLengths[reg$seq]
      }
      bamWhich(bamParam) = GRanges(seqnames=reg$seq, ranges=IRanges(start=reg$start, end=reg$end))
      filterBam(bf, "local.bam", param=bamParam)
      ezSystem(paste("samtools", "index", "local.bam"))
    } else {
      ezSystem(paste("cp", bf, "local.bam"))
      ezSystem(paste("cp", paste0(bf, ".bai"), "local.bam.bai"))
    }
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " AddOrReplaceReadGroups",
                " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "local.bam",
                " O=withRg.bam SORT_ORDER=coordinate",
                " RGID=RGID_", sampleName, " RGPL=illumina RGSM=", sampleName, " RGLB=RGLB_", sampleName, " RGPU=RGPU_", sampleName,
                " VERBOSITY=WARNING",
                " > addreplace.stdout 2> addreplace.stderr")
    ezSystem(cmd)
    file.remove("local.bam")
    
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " MarkDuplicates ",
                " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "withRg.bam",
                " O=", "dedup.bam",
                " REMOVE_DUPLICATES=false", ## do not remove, do only mark
                " ASSUME_SORTED=true",
                " VALIDATION_STRINGENCY=SILENT",
                " METRICS_FILE=" ,"dupmetrics.txt",
                " VERBOSITY=WARNING",
                " >markdup.stdout 2> markdup.stderr")
    ezSystem(cmd)
    file.remove("withRg.bam")
    
    #     ezSystem(paste("samtools", "index", "dedup.bam"))
    #     gatk = paste(javaCall, "-jar", "$GATK_jar")
    #     cmd = paste(gatk, "-T SplitNCigarReads", "-R", genomeSeq,
    #                 "-I", "dedup.bam",
    #                 "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
    #                 "-o", obf,
    #                 "> splitncigars.stdout 2> splitncigars.stderr") 
    #     ezSystem(cmd)
    #     file.remove("dedup.bam")
    ezSystem(paste("mv", "dedup.bam", obf))
    ezSystem(paste("samtools", "index", obf))
    setwd("..")
    return(obf)
  }, mc.cores=nBamsInParallel, mc.preschedule=FALSE)
  
  mpileupCmd = paste("bcftools", "mpileup","-Ou",
                     "-f", param$ezRef["refFastaFile"],
                     param$mpileupOptions,
                     ifelse(param$region == "", "", paste("--region", param$region)),
                     paste(bamFilesClean, collapse=" "))
  callCmd = paste("bcftools", "call",
                  "-Ou",
                  param$callOptions,
                  "-") ## read from stdin
  filterCmd = paste("bcftools", "filter",
                    "--output-type z",
                    "--output", basename(vcfOutputFile),
                    param$filterOptions,
                    "-") ## read from stdin
  ezSystem(paste(mpileupCmd, "|", callCmd, "|", filterCmd))
  indexTabix(basename(vcfOutputFile),format = "vcf")
  
  ## write an igv link
  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                  vcfUrls = paste(PROJECT_BASE_URL, vcfOutputFile, sep="/") )
  writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", vcfOutputFile),
               sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  
  
  ## create an html report
  setwdNew(reportDir)
  titles = list()
  titles[["VCF-Report"]] = paste("VCF-Report:", param$name)
  doc = openBsdocReport(titles[[length(titles)]])
  
  addDataset(doc, bamDataset, param)
  
  chromSizes = ezChromSizesFromVcf(file.path("..", basename(vcfOutputFile)))
  genotype = geno(readVcf(file.path("..", basename(vcfOutputFile)), genome="genomeDummy"))
  gt = genotype$GT
  gt[genotype$DP < param$minReadDepth] = "lowCov" ## those calls will become NA in subsequent analyses
  nSamples = nrow(bamDataset)
  conds = ezConditionsFromDataset(bamDataset, param=param)
  sampleColors = getSampleColors(conds, colorNames = names(conds))
  
  titles[["IGV"]] = "IGV"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file="igvSession.xml", vcfUrls = paste(PROJECT_BASE_URL, vcfOutputFile, sep="/") )
  writeIgvJnlp(jnlpFile="igv.jnlp", projectId = sub("\\/.*", "", output$getColumn("Report")),
               sessionUrl = paste0(PROJECT_BASE_URL, output$getColumn("Report"), "igvSession.xml"))
  addTxtLinksToReport(doc, "igv.jnlp", mime="application/x-java-jnlp-file")
  
  titles[["Sample Clustering based on Variants"]] = "Sample Clustering based on Variants"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  if (nrow(bamDataset) > 3){
    idxMat = ezMatrix(match(gt, c("0/0", "0/1", "1/1")) -2, rows=rownames(gt), cols=colnames(gt))
    d = dist(t(idxMat))
    if (all(!is.na(d))){
      hc=hclust(d, method="ward.D2" );
      hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
      hcd = colorClusterLabels(hcd, sampleColors)
      pngFile = "genotype-cluster.png"
      plotCmd = expression({
        plot(hcd, main="Cluster by Genotype", xlab="")
      })
      pngLink = ezImageFileLink(plotCmd, file=pngFile, width=800 + max(0, 10 * (nSamples-20))) # NOTEP: should the width be more sample-dependent?
      addParagraph(doc, pngLink)
    }
  }
  titles[["Variants by Chromosomes"]] = "Variants by Chromosomes"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  if (nrow(gt) > 0){
    chrom = sub(":.*", "", rownames(gt))
    pos = as.integer(sub("_.*", "", sub(".*:", "", rownames(gt))))
    isRealChrom = !grepl("[\\._]", names(chromSizes)) & chromSizes > 20000 ## TODO select chromosomes by name
    idxList = split(1:nrow(gt), chrom)
    snpColors = c("0/0"="blue", "0/1"="darkgrey", "1/1"="red")
    pngFiles = c()
    pngLinks = list()
    chromUse = sort(chromSizes[isRealChrom], decreasing = TRUE)
    for (ch in names(chromUse)){
      pngFiles[ch] = paste0("variantPos-chrom-", ch, ".png")
      plotCmd = expression({
        par(mar=c(4.1, 10, 4.1, 2.1))
        plot(0, 0, type="n", main=paste("Chromsome", ch), xlab="pos", xlim=c(1, chromSizes[ch]), ylim=c(0, 3*ncol(gt)),
             axes=FALSE, frame=FALSE, xaxs="i", yaxs="i", ylab="")
        axis(1)
        mtext(side = 2, at = seq(1, 3*ncol(gt), by=3), text = colnames(gt), las=2,
              cex = 1.0, font=2, col=sampleColors)
        idx = idxList[[ch]]
        xStart = pos[idx]
        nm  = colnames(gt)[1]
        for (i in 1:ncol(gt)){
          offSet = match(gt[idx ,i], names(snpColors))
          yTop = (i-1) * 3 + offSet
          rect(xStart, yTop - 1, xStart+1, yTop, col = snpColors[offSet], border=snpColors[offSet])
        }
        abline(h=seq(0, 3*ncol(gt), by=3))
      })
      pngLinks[ch] = ezImageFileLink(plotCmd, file=pngFiles[ch], height=150+25*ncol(gt), width=1000)
    }
    if (length(pngLinks) > 0){
      addFlexTable(doc, ezGrid(pngLinks))
    }
  }
  closeBsdocReport(doc, htmlFile, titles)
  setwd("..")
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMpileup(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMpileup <-
  setRefClass("EzAppMpileup",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMpileup
                  name <<- "EzAppMpileup"
                  appDefaults <<- rbind(region=ezFrame(Type="character",  DefaultValue="",  Description="should analysis be done on a region of the genom only chr:start-end"),
                                    mpileupOptions=ezFrame(Type="character",  DefaultValue="",	Description="options to mpileup command"),
                                    callOptions=ezFrame(Type="character", DefaultValue="", Description="options to call command"),
                                    filterOptions=ezFrame(Type="character", DefaultValue="", Description="options to filter command"),
                                    minReadDepth=ezFrame(Type="integer", DefaultValue="10", Description="use for clustering only SNV with coverage higher than"))
                }
              )
  )
