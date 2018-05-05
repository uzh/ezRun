###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGatkRnaHaplotyper = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
  bamDataset = input$meta
  
  reportDir = basename(output$getColumn("Report"))
  htmlFile = basename(output$getColumn("Html"))
  vcfOutputFile = output$getColumn("VCF [File]") #paste0(param$name, "-haplo.vcf.gz")
  
  bamFiles = input$getFullPaths("BAM")
  names(bamFiles) = rownames(bamDataset)
  genomeSeq = param$ezRef["refFastaFile"]
  nBamsInParallel = min(4, param$cores)
  bamFilesClean = ezMclapply(names(bamFiles), function(sampleName){
    javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", floor(param$ram/nBamsInParallel), "g")
    setwdNew(paste(sampleName, "proc", sep="-"))
    bf = bamFiles[sampleName]
    obf = file.path(getwd(), basename(bf))
    if (ezIsSpecified(param$seqNames)){
      bamParam = ScanBamParam(what=scanBamWhat())
      seqLengths = ezBamSeqLengths(bf)
      bamWhich(bamParam) = GRanges(seqnames=param$seqNames, ranges=IRanges(start=1, end=seqLengths[param$seqNames]))
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
    
    ezSystem(paste("samtools", "index", "dedup.bam"))
    gatk = paste(javaCall, "-jar", "$GATK_jar")
    cmd = paste(gatk, "-T SplitNCigarReads", "-R", genomeSeq,
                "-I", "dedup.bam",
                "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
                "-o", obf,
                "> splitncigars.stdout 2> splitncigars.stderr") 
    ezSystem(cmd)
    file.remove("dedup.bam")
    ezSystem(paste("samtools", "index", obf))
    setwd("..")
    return(obf)
  }, mc.cores=nBamsInParallel, mc.preschedule=FALSE)
  
  ########### haplotyping
  haplotyperCall = paste0("java -Djava.io.tmpdir=. -Xmx", param$ram, "g", " -jar ", "$GATK_jar", " -T HaplotypeCaller")
  cmd = paste(haplotyperCall, "-R", genomeSeq,
              "-nct", param$cores,
              paste("-I", bamFilesClean, collapse=" "),
              "--dontUseSoftClippedBases", ##--recoverDanglingHeads
              "-stand_call_conf 20 -stand_emit_conf 20",
              "--activeRegionOut", paste0(param$name, "-haplo-activeRegion.bed"),
              "--maxNumHaplotypesInPopulation", 4,
              #"--maxReadsInRegionPerSample", "10000",
              "--sample_ploidy", 2,
              "--min_mapping_quality_score 1",
              "--output_mode", "EMIT_VARIANTS_ONLY",   ## does not work: EMIT_ALL_CONFIDENT_SITES
              "-o", paste0(param$name, "-all-haplo.vcf"),
              ">", paste0(param$name, "-haplo.stdout"),
              "2>", paste0(param$name, "-haplo.stderr"))
  ezSystem(cmd)
  
  ## filter the vcf file
  require("VariantAnnotation")
  ezFilterVcf(vcfFile=paste0(param$name, "-all-haplo.vcf"), basename(vcfOutputFile), discardMultiAllelic=FALSE,
              bamDataset=bamDataset, param=param)
  gc()
  
  ## create an html report
  setwdNew(reportDir)
  doc = openBsdocReport(paste("VCF-Report:", param$name))
  addDataset(doc, bamDataset, param)
  
  chromSizes = ezChromSizesFromVcf(file.path("..", basename(vcfOutputFile)))
  genotype = geno(readVcf(file.path("..", basename(vcfOutputFile)), genome="genomeDummy"))
  gt = genotype$GT
  gt[genotype$DP < param$vcfCall.minReadDepth] = "lowCov"
  nSamples = nrow(bamDataset)
  
  titles[["IGV"]] = "IGV"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file="igvSession.xml", vcfUrls = paste(PROJECT_BASE_URL, vcfOutputFile, sep="/") )
  writeIgvJnlp(jnlpFile="igv.jnlp", projectId = sub("\\/.*", "", output$Report),
               sessionUrl = paste(PROJECT_BASE_URL, output$Report, "igvSession.xml", sep="/"))
  addTxtLinksToReport(doc, "igv.jnlp", mime = "application/x-java-jnlp-file")
  
  titles[["Sample Clustering based on Variants"]] = "Sample Clustering based on Variants"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  conds = ezConditionsFromDataset(bamDataset, param=param)
  sampleColors = getSampleColors(conds, colorNames = names(conds))
  idxMat = ezMatrix(match(gt, c("0/0", "0/1", "1/1")) -2, rows=rownames(gt), cols=colnames(gt))
  d = dist(t(idxMat))
  hc = hclust(d, method="ward.D2");
  hcd = as.dendrogram(hclust(d, method="ward.D2"), hang=-0.1)
  hcd = colorClusterLabels(hcd, sampleColors)
  pngFile = "genotype-cluster.png"
  plotCmd = expression({
    plot(hcd, main="Cluster by Genotype", xlab="")
  })
  pngLink = ezImageFileLink(plotCmd, file=pngFile, width=800 + max(0, 10 * (nSamples-20))) # nSamples dependent width
  addParagraph(doc, pngLink)
  
  titles[["Variants by Chromosomes"]] = "Variants by Chromosomes"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addParagraph(doc, "Genotype colors are: blue - homozygous reference; gray - heterozygous; red - homozygyous variant")
  chrom = sub(":.*", "", rownames(gt))
  pos = as.integer(sub("_.*", "", sub(".*:", "", rownames(gt))))
  isRealChrom = !grepl("[\\._]", names(chromSizes)) ## TODO select chromosomes by name
  idxList = split(1:nrow(gt), chrom)
  snpColors = c("0/0"="blue", "0/1"="darkgrey", "1/1"="red")
  pngFiles = c()
  pngLinks = character()
  for (ch in names(chromSizes)[isRealChrom]){
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
    pngLinks[ch] = ezImageFileLink(plotCmd, file=pngFiles[ch], height=200+30*ncol(gt), width=1200)
  }
  addFlexTable(doc, ezGrid(pngLinks))
  closeBsdocReport(doc, htmlFile, titles)
  setwd("..")
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkRnaHaplotyper(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppGatkRnaHaplotyper <-
  setRefClass("EzAppGatkRnaHaplotyper",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGatkRnaHaplotyper
                  name <<- "EzAppGatkRnaHaplotyper"
                  appDefaults <<- rbind(vcfFilt.minAltCount = ezFrame(Type="integer",  DefaultValue=10,  Description="minimum coverage for the alternative variant"),
                                        vcfCall.minReadDepth = ezFrame(Type="integer",  DefaultValue=10,  Description="minimum read deapth"))
                }
              )
  )
