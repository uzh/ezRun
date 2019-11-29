## 




# minReadsPerGene=3 <- is needed so that it works for iSeq runs

#' compute and plot the GC content and gene-length associated bias in RNA-seq data
#'
#' @param dsFile 
#' @param dsName 
#' @param param 
#' @param qcSummaryDir 
#' @param refBuildMap 
#' @param minReadsPerSample 
#' @param maxReadsPerSample 
#' @param minReadsPerGene 
#' @param minPresentFraction 
#'
#' @return
#' @export
#'
#' @examples
#' ## needs the following modules loaded module add Tools/samtools QC/Trimmomatic QC/Flexbar Aligner/kallisto
#' options(error=recover)
#' require(ezRun)
#' dsFile = "/srv/gstore/projects/p3082/ISeq_20190425_iSeq35_o5404_DataDelivery/dataset.tsv"
#' ezComputeBias(dsFile)

ezComputeBias = function(dsFile, dsName=NULL, param=NULL, qcSummaryDir="/srv/GT/analysis/p2220/RNA-seq-bias-results", refBuildMap=getRefBuildMap(), minReadsPerSample=30000,
                       maxReadsPerSample=5e6,
                       minReadsPerGene=3, minPresentFraction=0.2){
  
  inputMeta = ezRead.table(file=dsFile)
  inputMeta = inputMeta[ inputMeta$"Read Count" > minReadsPerSample, ]
  inputMeta$Name = rownames(inputMeta)

  outMeta = inputMeta
  outMeta$'Count [File]' = paste0(rownames(outMeta), "-count.txt")
  outMeta$'bootstrappedCount [File]' = paste0(rownames(outMeta), "-bootstrap.h5")
  outMeta$'runInfo [File]' = paste0(rownames(outMeta), "-runInfo.json")
  outMeta$'BAM [File]' = paste0(rownames(outMeta), "-tr.bam")
  outMeta$'CovStat [File]' = paste0(rownames(outMeta), "-genebody_coverage.rds")
  
  param = list()
  param[['cores']] = '8'
  param[['ram']] = '30'
  param[['scratch']] = '100'
  param[['node']] = ''
  param[['process_mode']] = 'SAMPLE'
  param[['paired']] = 'false'
  param[['strandMode']] = 'both'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['bootstrap-samples']] = '10'
  param[['seed']] = '42'
  param[['fragment-length']] = '150'
  param[['sd']] = '70'
  param[['bias']] = 'true'
  param[['pseudobam']] = 'true'
  param[['transcriptFasta']] = ''
  param[['transcriptTypes']] = 'protein_coding'
  param[['trimAdapter']] = 'true'
  param[['trimLeft']] = '0'
  param[['trimRight']] = '0'
  param[['minTailQuality']] = '15'
  param[['minAvgQuality']] = '20'
  #param[['mail']] = 'hubert.rehrauer@fgcz.ethz.ch'
  param[['dataRoot']] = '/srv/gstore/projects'
  param[['backgroundExpression']] = minReadsPerGene ## low numbers needed for the iSeq runs with
  param$sigThresh = minReadsPerGene
  param[['transcriptTypes']] = 'protein_coding'
  param$normMethod = "logMean"
  param$expressionName = "est_counts"
  param$nReads = maxReadsPerSample
  paramList = param
  paramList$refBuild = refBuildMap[inputMeta$Species[1]]
  
  if (is.na(paramList$refBuild)){
    stop(paste("refBuild not defined for species:", refBuildMap[inputMeta$Species[1]]))
  }
    
  if (is.null(dsName)){
    ## assumes the dsFile has the format: /srv/gstore/projects/<project>/<read folder>
    parentDirs = strsplit(dirname(dsFile), "/")[[1]]
    dsName = paste0(parentDirs[5], "_", parentDirs[6])
    plates = paste(unique(inputMeta$`PlateName [Characteristic]`), collapse="_")
    if (!is.null(plates) && plates != ""){
      dsName = paste0(plates, "_", dsName)
    }
  }
  message(dsName)
  resultDir = file.path(qcSummaryDir, dsName)
  setwdNew(resultDir)

  ## align and count with kallisto + coverage profiles
  for (i in 1:nrow(inputMeta)){
    message(i)
    covStatFile = file.path(resultDir, outMeta[i, "CovStat [File]"])
    countFile = file.path(resultDir, outMeta$`Count [File]`[i])
    if (file.exists(covStatFile) && file.exists(countFile)){
      next
    }
    scratchDir = paste0("/scratch/", dsName, "-", rownames(inputMeta)[i], "-", ezTime())
    setwdNew(scratchDir)
    EzAppKallisto$new()$run(input=inputMeta[i, ], output=outMeta[i, ], param=paramList)
    ezSystem(paste("mv", outMeta$`Count [File]`[i], resultDir))
    ga = ezReadGappedAlignments(outMeta[i, "BAM [File]"])
    tcList = coverage(ga)
    transcriptLengthTotal <- elementNROWS(tcList)
    sampledTranscriptCov <- ezMclapply(tcList,
                                       function(x){as.integer(x[round(seq(1, length(x), length.out=101))])},
                                       mc.preschedule=TRUE, mc.cores=min(as.integer(param$cores), 4))
    sampledTranscriptCov <- do.call(cbind, sampledTranscriptCov)
    trUse = colSums(sampledTranscriptCov) > 0
    sampledTranscriptCov = sampledTranscriptCov[ , trUse, drop=FALSE]
    trLength = transcriptLengthTotal[trUse]
    lengthClasses = ezCut(trLength, breaks=c(399, 1000, 1500, 4000), 
                          labels=c("less than 400nt", "400 to 999nt", 
                                   "1000 to 1500nt", "1501 to 4000", "above 4000nt"))
    genebody_coverage = list()
    for (lc in levels(lengthClasses)){
      isInLc = lengthClasses == lc
      if (sum(isInLc) > 50){
        ltc = sampledTranscriptCov[ , isInLc, drop=FALSE]
        avgCov = colMeans(ltc)
        relativeCov = ezScaleColumns(ltc, 1/colSums(ltc)) ## normalize so that every transcripts adds the same weight
        avgCovQuant = unique(quantile(avgCov, c(0.3, 0.95)))
        if (length(avgCovQuant) == 2){
          covClasses = ezCut(avgCov, breaks=avgCovQuant,
                             labels=c("low expressed", "medium expressed", 
                                      "high expressed"))
          genebody_coverage[[lc]] = list()
          for (cc in levels(covClasses)){
            genebody_coverage[[lc]][[cc]] = rowMeans(relativeCov[ , covClasses == cc, drop=FALSE ])
          }
        }
      }
    }
    saveRDS(genebody_coverage, covStatFile)
    setwd(resultDir)
    unlink(scratchDir, recursive = TRUE)
  }
    
  
  ## compute the stats from the 
  param = ezParam(paramList)
  setwd(qcSummaryDir)
  
    
  countDs = outMeta
  countDs$featureLevel = "isoform"
  countDs = countDs[file.exists(file.path(resultDir, countDs$`Count [File]`)), ]
  input = EzDataset(meta=countDs, dataRoot=resultDir)
  rawData = loadCountDataset(input, param)
  rawData = subset(rawData, select=colSums(assays(rawData)$counts >= minReadsPerGene) >= 100)
  counts = assays(rawData)$counts
  nGenes = colSums(counts >= param$sigThresh)
  meta = ezFrame(colData(rawData))
  assays(rawData)$signal = ezNorm(assays(rawData)$counts,
                                  presentFlag=assays(rawData)$counts >= 3,
                                  method=param$normMethod)
  seqAnno <- data.frame(rowData(rawData), row.names=rownames(rawData),
                        check.names = FALSE, stringsAsFactors=FALSE)
  x = assays(rawData)$signal
  stopifnot(!is.na(x))
  isPresent = rowMeans(x > 1) >= minPresentFraction
  logSig = log2(x + param$backgroundExpression)
  logRatio = logSig - rowMeans(logSig)
  isHighGc = seqAnno$gc > 0.6 & isPresent
  isLowGc = seqAnno$gc < 0.4 & isPresent
  gcEffect = apply(logRatio, 2, function(x){mean(x[isHighGc], na.rm =TRUE) - mean(x[isLowGc] ,na.rm = TRUE)})
  gcTypes = ezFrame(lowGc=isLowGc, highGc=isHighGc)
  widthTypes = data.frame("width < 800nt"=as.numeric(rowData(rawData)$featWidth) < 800, 
                          "width > 4000nt"=as.numeric(rowData(rawData)$featWidth) > 4000,
                          check.names=FALSE)
  widthEffect = apply(logRatio, 2, function(x){mean(x[widthTypes$`width > 4000nt`], na.rm =TRUE) - mean(x[widthTypes$`width < 800nt`] ,na.rm = TRUE)})
    
  if (!is.null(meta$`PlatePosition [Characteristic]`)){
    ptLabels = sub(".*_", "", meta$`PlatePosition [Characteristic]`)
  } else {
    ptLabels = 1:nrow(meta)
  }
  
  countDs$"Read Count"[countDs$"Read Count" > param$nReads] = param$nReads
    
  widthEffect = shrinkToRange(widthEffect, c(-1, 1))
  gcEffect = shrinkToRange(gcEffect, c(-1, 1))
    
  ## plot the 3'-enrichment
  biasScores = sapply(countDs$`CovStat [File]`, function(rdsFile){
    gbCov = readRDS(file.path(resultDir, rdsFile))
    midScore = sum(gbCov$`1000 to 1500nt`$`medium expressed`[80:101])
    longScore = sum(gbCov$`1501 to 4000`$`medium expressed`[80:101])
    c(midScore, longScore)
  })
  rownames(biasScores) = c("medium genes", "long genes")
  colnames(biasScores) = rownames(countDs)
  biasScores = t(biasScores)
  biasScores = shrinkToRange(biasScores, c(0.1, 0.5))
  
  pdf(file=paste0(dsName, "-bias.pdf"), width=16, height=8)
  par(mfrow=c(2,4))
  par(pty="s")
  
  smoothScatter(widthEffect, gcEffect, xlim=c(-1, 1), ylim=c(-1, 1), main=dsName, nrpoints = 0, nbin=40, bandwidth = 0.05,
                cex.main=0.8, colramp=colorRampPalette(c("white", blues9[1:6])))
  text(widthEffect, gcEffect, labels = ptLabels, cex=0.8)

  dsty = density(gcEffect, bw=0.05, from=-1, to=1)
  plot(dsty$y, dsty$x, type="l", xlab="density", ylab="", main=countDs$LibraryPrepKit[1], cex.main=0.8)
    
  plateRow = substr(ptLabels, start = 1, stop=1)
  plateCol = substr(ptLabels, start = 2, stop=3)
  mat = ezMatrix(NA, rows=toupper(letters[1:8]), cols=as.character(1:12))
  for (i in which(plateRow %in% rownames(mat) & plateCol %in% colnames(mat))) {
    mat[plateRow[i], plateCol[i]] = gcEffect[i]
  }
  lim = c(-1, 1)
  image(1:ncol(mat), 1:nrow(mat), t(shrinkToRange(mat, lim)), col=getBlueRedScale(), breaks = seq(from = lim[1], to = lim[2], length.out = 257), axes=FALSE,
        xlim=c(0.5, ncol(mat)+0.5), ylim=c(0.5, nrow(mat)+0.5), xlab="", ylab="")
  axis(1, 1:12, labels = colnames(mat))
  axis(2, 1:8, labels = rownames(mat), las=2)
    
  logWidthBreaks = c(9.5, 10.5, 11.5)
  gcBreaks = c(0.42, 0.48, 0.53, 0.57, 0.62)
  #gcBreaks = c(0.4, 0.43, 0.46, 0.49, 0.53, 0.57, 0.61, 0.64, 0.67)
  ctsCorr = ezCorrectBias(counts, seqAnno$gc, seqAnno$featWidth, logWidthBreaks = logWidthBreaks, gcBreaks = gcBreaks, widthOffset=100)
  lwClasses = ezCut(ctsCorr$logWidthShrinked, breaks = logWidthBreaks)
  gcClasses = ezCut(ctsCorr$gcShrinked, gcBreaks)
  mat = ctsCorr$binScoresFilled
  lim = c(-1.5, 1.5)
  image(1:nrow(mat), 1:ncol(mat), shrinkToRange(mat, lim), col=getBlueRedScale(), breaks = seq(from = lim[1], to = lim[2], length.out = 257), axes=FALSE,
        xlim=c(0.5, nrow(mat)+0.5), ylim=c(0.5, ncol(mat)+0.5), ylab="GC", xlab="width")
  

  ## second row of plots
  plot(1:length(gcEffect), gcEffect, type="n", ylim=c(-1, 1))
  text(1:length(gcEffect), gcEffect, labels = ptLabels, cex=0.8)
  
  plot(1:length(gcEffect), gcEffect, type="n", ylim=c(-1, 1), xlim=c(1, length(gcEffect)*1.1))
  text(1:length(gcEffect), gcEffect, labels = names(gcEffect), cex=0.8, pos = 4)
  
  smoothScatter(biasScores[ , "medium genes"], biasScores[ , "long genes"], xlim=c(0.1, 0.5), ylim=c(0.1, 0.5), 
                xlab="medium genes", ylab="long genes", main=dsName, nrpoints = 0, nbin=40, bandwidth = 0.05,
                cex.main=0.8, colramp=colorRampPalette(c("white", blues9[1:6])))
  abline(h=0.2, v=0.2, col="gray")
  text(biasScores[ , "medium genes"], biasScores[ , "long genes"], labels = ptLabels, cex=0.8)
    
  readCount = countDs[names(nGenes), "Read Count"] * exp(runif(length(nGenes), min=-0.1, max=0.1))
  plot(readCount,
       nGenes, log="xy",
       xlab="nReads", ylab="nGenes", main=dsName, type="n")
  text(readCount, 
       nGenes, labels=ptLabels, cex=0.8)
  dev.off()
}


getRefBuildMap = function(){
  refBuildMap = c(
    "Mus musculus (house mouse)" = "Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26",
    "Homo sapiens (human)"="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26",
    "human"="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26",
    "homo"="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26",
    "Equus caballus (horse)"="Equus_caballus/Ensembl/EquCab2/Annotation/Release_89-2017-06-07",
    "Canis familiaris"="Canis_familiaris/Ensembl/CanFam3.1/Annotation/Release_89-2017-06-07",
    "Oryctolagus cuniculus (rabbit)"="Oryctolagus_cuniculus/Ensembl/OryCun2.0/Annotation/Release_89-2017-06-08",
    "Escherichia coli K12"="Escherichia_coli/Ensembl/K12_MG1655_ASM584v2/Annotation/Release_2019-01-31",
    "E Coli"="Escherichia_coli/Ensembl/K12_MG1655_ASM584v2/Annotation/Release_2019-01-31",
    "Fusarium oxysporum"="Fusarium_oxysporum/Misc/Fo5176/Annotation/Release_5176-2019-03-18",
    "Canis lupus familiaris"="Canis_familiaris/Ensembl/CanFam3.1/Annotation/Release_89-2017-06-07",
    "Cricetulus griseus (Chinese hamster)"="Cricetulus_griseus/CHO/CHO-K1/Annotation/Version-2015-12-11/Genes",
    "Felis silvestris catus"="Felis_catus/Ensembl/Felis_catus_6.2/Annotation/Release_89-2017-06-07",
    "Candida glabrata"="Candida_albicans/CGD/SC5314_Assembly22/Annotation/Release_2018-11-13",
    "Toxoplasma gondii"="Toxoplasma_gondiiME49/ToxoDB/TgondiiME49_R33/Annotation/Release_33-2017-07-28",
    "Sus Linnaeus (Pig)"="Sus_scrofa/Ensembl/Sscrofa11.1/Annotation/Release_90-2017-11-17",
    "Canis familiaris"="Canis_familiaris/Ensembl/CanFam3.1/Annotation/Release_89-2017-06-07",
    "Sus scrofa"="Sus_scrofa/Ensembl/Sscrofa11.1/Annotation/Release_90-2017-11-17",
    "Bos taurus"="Bos_taurus/Ensembl/UMD_v3.1/Annotation/Release_92-2018-05-30",
    "Rattus norvegicus"="Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Release_98-2019-11-11",
    "Gallus gallus (chicken)"="Gallus_gallus/Ensembl/Gallus_gallus-5.0/Annotation/Release_89-2017-06-08",
    "Danio rerio (zebrafish)"="Danio_rerio/Ensembl/GRCz11/Annotation/Release_96-2019-04-22",
    "Drosophila melanogaster (fruit fly)"="Drosophila_melanogaster/Ensembl/BDGP6/Annotation/Release_89-2017-06-07",
    "Caenorhabditis elegans (nematode)"="Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Release_90-2017-08-31",
    "Saccharomyces cerevisiae (bakers yeast)"="Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Version-2015-08-04",
    "Oryza sativa (rice)"="Oryza_sativa_japonica/Ensembl/IRGSP_1.0/Annotation/Release_42-2019-01-28",
    "Arabidopsis thaliana (thale cress)"="Arabidopsis_thaliana/TAIR/TAIR10/Annotation/Release_01-2018-10-09",
    "Zea mays"="Sus_scrofa/Ensembl/Sscrofa11.1/Annotation/Release_90-2017-11-17",
    "Triticum aestivum"="Triticum_aestivum/URGI/IWGSCv3/Annotation/Version-2016-05-09",
    "Escherichia coli"="Escherichia_coli/Ensembl/K12_MG1655_ASM584v2/Annotation/Release_2019-01-31"
  )
  return(refBuildMap)
}
