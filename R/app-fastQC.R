###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Fast QC
##' @template htmlFile-template
##' @seealso \code{\link{EzAppFastqc}}
ezMethodFastQC = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  dataset = input$meta
  samples = rownames(dataset)
  files = c()
  for (sm in samples){
    files[paste0(sm, "_R1")] = input$getFullPaths(param,"Read1")[sm]
    if (!is.null(dataset$Read2)){
      files[paste0(sm, "_R2")] = input$getFullPaths(param,"Read2")[sm]    
    }
  }
  nFiles = length(files)
  
  ## guess the names of the report directories that will be creatd by fastqc
  reportDirs = sub(".fastq.gz", "_fastqc", basename(files))
  reportDirs = sub(".fq.gz", ".fq_fastqc", reportDirs)
  stopifnot(!duplicated(reportDirs))
  filesUse = files[!file.exists(reportDirs)]
  if (length(filesUse) > 0){
    cmd = paste(FASTQC, "--extract -o . -t", min(ezThreads(), 8), param$cmdOptions,
                paste(filesUse, collapse=" "),
                "> fastqc.out", "2> fastqc.err")
    result = ezSystem(cmd)
  }
  statusToPng = c(PASS="tick.png", WARN="warning.png", FAIL="error.png")
  
  ## collect the overview table
  plots = c("Per base sequence quality"="per_base_quality.png",
            "Per sequence quality scores"="per_sequence_quality.png",
            "Per tile sequence quality"="per_tile_quality.png",
            "Per base sequence content"="per_base_sequence_content.png",
            "Per sequence GC content"="per_sequence_gc_content.png",
            "Per base N content"="per_base_n_content.png",
            "Sequence Length Distribution"="sequence_length_distribution.png",
            "Sequence Duplication Levels"="duplication_levels.png",
            "Adapter Content"="adapter_content.png",
            "Kmer Content"="kmer_profiles.png")
  
  ## make for each plot type an html report with all samples
  plotPages = sub(".png", ".html", plots)
  for (i in 1:length(plots)){
    plotDoc = openBsdocReport(title=paste("FASTQC:", plotPages[i]))
    png = paste0("<img src=", reportDirs, "/Images/", plots[i], ">")
    tbl = ezFlexTable(ezMatrix(png, rows=names(files), cols=names(plots)[i]), 
                      header.columns=TRUE, add.rownames=TRUE, valign="middle")
    plotDoc = addFlexTable(plotDoc, tbl)
    closeBsdocReport(plotDoc, plotPages[i])
  }
  
  ## establish the main report
  titles = list()
  titles[["FastQC"]] = paste("FASTQC:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]], dataset=dataset)
  
  titles[["Read Counts"]] = "Read Counts"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  if (!is.null(dataset$"Read Count")){
    readCount = signif(dataset$"Read Count" / 1e6, digits=3)
    names(readCount) = rownames(dataset)
  } else {
    readCount = integer()
    for (i in 1:nFiles){
      x = ezRead.table(file.path(reportDirs[i], "fastqc_data.txt"), header=FALSE, nrows=7, fill=TRUE)
      readCount[names(files)[i]] = signif(as.integer(x["Total Sequences", 1]) / 1e6, digits=3)
    }
  }
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    barplot(readCount, las=2, ylab="Counts [Mio]", main="total reads")
  })
  readCountsLink = ezImageFileLink(plotCmd, file="readCounts.png", height=600,
                                   width=min(600 + (nFiles-10)* 30, 2000))
  doc = addParagraph(doc, readCountsLink)
  
  titles[["Fastqc quality measures"]] = "Fastqc quality measures"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  statusToPng = c(PASS="tick.png", WARN="warning.png", FAIL="error.png")
  for (i in 1:nFiles){
    smy = ezRead.table(file.path(reportDirs[i], "summary.txt"), row.names=NULL, header=FALSE)
    if (i == 1){
      rowNames = paste0("<a href=", reportDirs, "/fastqc_report.html>", names(files), "</a>")
      colNames = ifelse(smy[[2]] %in% names(plotPages),
                        paste0("<a href=", plotPages[smy[[2]]], ">", smy[[2]], "</a>"),
                        smy[[2]])
      tbl = ezMatrix("", rows=rowNames, cols=colNames)
    }
    href = paste0(reportDirs[i], "/fastqc_report.html#M", 0:(ncol(tbl)-1))
    img = paste0(reportDirs[i], 	"/Icons/", statusToPng[smy[[1]]])
    tbl[i, ] = paste0("<a href=", href, "><img src=", img, "></a>")
  }
  doc = addFlexTable(doc, ezFlexTable(tbl, header.columns=TRUE, add.rownames=TRUE, valign="middle"))
  
  titles[["Per Base Read Quality"]] = "Per Base Read Quality"
  addTitleWithAnchor(doc, titles[[length(titles)]], 2)
  qualMatrixList = ezMclapply(files, getQualityMatrix, mc.cores=ezThreads())
  pngMatrix = plotQualityMatrixAsHeatmap(qualMatrixList, isR2=grepl("_R2", names(files)))
  doc = addFlexTable(doc, ezGrid(pngMatrix))
  if(nrow(dataset) > 1){
    pngLibCons = character()
    pngLibCons["qPCR"] = plotReadCountToLibConc(dataset, colname='LibConc_qPCR [Characteristic]')
    pngLibCons["100_800bp"] = plotReadCountToLibConc(dataset, colname='LibConc_100_800bp [Characteristic]')
    if(length(pngLibCons) > 0){
      titles[["Correlation"]] = "Correlation between Library concentration measurements and ReadCounts"
      addTitleWithAnchor(doc, titles[[length(titles)]], 3)
      doc = addFlexTable(doc, ezGrid(matrix(pngLibCons, nrow=1)))
    }
  }
  closeBsdocReport(doc, htmlFile, titles)
  ezSystem(paste("rm -rf ", paste0(reportDirs, ".zip", collapse=" ")))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastQC()
##' @seealso \code{\link{ezMethodFastQC}}
EzAppFastqc <-
  setRefClass("EzAppFastqc",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFastQC
                  name <<- "EzAppFastqc"
                }
              )
  )

##' @describeIn ezMethodFastQC Plots \code{colname} from \code{dataset} against read counts in millions.
plotReadCountToLibConc = function(dataset, colname){
  if(colname %in% colnames(dataset) && nrow(dataset) > 1){
    if(!all(dataset[[colname]]==0)){
      dataset = dataset[order(dataset$'Read Count',decreasing = T),]
      dataset$'Read Count' = dataset$'Read Count'/10^6
      corResult = cor.test(dataset$'Read Count',dataset[[colname]],method = 'spearman')
      regressionResult = lm(dataset[[colname]]~dataset$'Read Count')
      label = sub(' \\[.*','',colname)
      plotCmd = expression({
        plot(dataset$'Read Count', dataset[[colname]], pch=c(18), cex=1.5, main=label,
             xlab='ReadCount in Mio', ylab=sub('\\[.*','', colname), xlim=c(0, max(dataset$'Read Count', na.rm=TRUE)*1.2), #min(dataset$'Read Count')*0.8
             ylim=c(min(dataset[[colname]], na.rm = TRUE) * 0.8, max(dataset[[colname]], na.rm=TRUE) * 1.2))
        legend("topright", paste('r=', round(corResult$estimate, 2)), bty="n") 
        abline(regressionResult, col='red',lty=c(2))
        text(dataset$'Read Count', dataset[[colname]], pos=2,
             labels=rownames(dataset), cex=0.7, col='darkcyan')
      })
      link = ezImageFileLink(plotCmd, file=paste0('ReadCount_', label, '.png'), width=500, height=500)
      return(link)
    }
  }
}

##' @describeIn ezMethodFastQC Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.
getQualityMatrix = function(inputFile){
  ## files could be fastq, or gziped fastq.gz
  library(ShortRead, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)  
  #job = ezJobStart(paste("start to collect quality matrix from", inputFile))
  
  ## get the qualCountMatrix
  qualCountMatrix = NULL
  maxReadLength = NULL
  subSample = 0.01
  fqs = FastqStreamer(inputFile, 1e6)
  while(length(x <- yield(fqs))){
    nSamples = round(subSample * length(x))
    if (nSamples == 0){
      nSamples = length(x)
    }
    qual = quality(x[sample.int(length(x), size=nSamples , replace=FALSE)]) ## this gives the integer quality as stringSets
    readLengths = width(qual)
    #qualMatrix = as(qual, "matrix") ## this should be now a matrix with integer quality values, only work with same read length
    if (is.null(maxReadLength)){
      maxReadLength = max(readLengths)
    }
    qualMatrix = ezMatrix(NA, cols=1:maxReadLength, rows=1:length(qual))
    availableLengths = unique(readLengths)
    useLengths = intersect(availableLengths, 1:maxReadLength)
    ## This is the fastest way to turn the quality into a integer matrix with NAs if some reads with shorted length
    for(l in useLengths){
      index = readLengths == l
      qualMatrix[index, 1:l] = as(qual[index], "matrix") 
    }
    
    if (is.null(qualCountMatrix)){
      maxQuality = max(qualMatrix, na.rm=TRUE)
      ## min quality is 0!!
      ## the qualCountMatrix has dimensions read quality * read length
      qualCountMatrix = ezMatrix(0, rows=0:maxQuality, cols=1:maxReadLength)
    }
    qualMatrix[qualMatrix > maxQuality] = maxQuality
    for (basePos in 1:ncol(qualMatrix)){
      qualCountByPos = table(qualMatrix[ , basePos])
      idx = names(qualCountByPos)
      qualCountMatrix[idx, basePos] = qualCountMatrix[idx, basePos] + qualCountByPos
    }
  }
  close(fqs)
  #ezWriteElapsed(job)
  return(qualCountMatrix)
}

##' @describeIn ezMethodFastQC Returns a png table of quality matrices interpreted as heatmaps.
plotQualityMatrixAsHeatmap = function(qualMatrixList, isR2=FALSE, xScale=1, yScale=1){
  pngFileNames = NULL
  colorsGray = gray((30:256)/256)
  minPercent = 0
  maxPercent = sqrt(40)
  minDiff = -5
  maxDiff = 5
  ## test if R2 exists
  index=list("R1"=which(!isR2))
  pngTable = data.frame(R1=character(0), stringsAsFactors=FALSE)
  if(any(isR2)){
    stopifnot(sum(isR2) == sum(!isR2))
    index[["R2"]]= which(isR2)
    pngTable[["R2"]] = character(0)
  }
  for(nm in names(index)){
    idx = index[[nm]]
    ## Plot the color key for the average quality heatmap R1_1
    colorKeyFile = paste0("averageReadsQuality-Key_", nm, ".png")
    by.label = 1
    at=seq(from=minPercent, to=maxPercent, by=by.label)
    plotCmd = expression({
      ezColorLegend(colorRange=c(minPercent, maxPercent), 
                    colors=colorsGray, vertical=FALSE, by.label=by.label,
                    at=at, labels=as.character(at^2))
    })
    colorKeyLink = ezImageFileLink(plotCmd, file=colorKeyFile, addPdfLink=FALSE, width=400*yScale, height=200*xScale)
    pngTable["Avg Qual Colors", nm] = colorKeyLink
    
    result = ezMatrix(0, dim=dim(qualMatrixList[[idx[1]]]))
    resultCount = result
    for(i in idx){
      qm = qualMatrixList[[i]]
      if (any(dim(qm) > dim(result))){
        oldResult = result
        result = ezMatrix(0, dim=dim(qm))
        result[1:nrow(oldResult), 1:ncol(oldResult)] = oldResult
        oldResultCount = resultCount
        resultCount = ezMatrix(0, dim=dim(qm))
        resultCount[1:nrow(oldResultCount), 1:ncol(oldResultCount)] = oldResultCount
      }
      result[1:nrow(qm), 1:ncol(qm)] = result[1:nrow(qm), 1:ncol(qm)] + qm
      resultCount[1:nrow(qm), 1:ncol(qm)] = resultCount[1:nrow(qm), 1:ncol(qm)] + 1
    }
    result = result / resultCount
    avgQual = signif(prop.table(result,2) * 100, digits=3)
    pngFileLink = plotQualityHeatmap(sqrt(avgQual), colorRange=c(minPercent, maxPercent), 
                                     pngFileName=paste0("averageReadsQuality-heatmap_", nm, ".png"),
                                     colors=colorsGray, main=paste("averageReadsQuality", nm, sep="_"), 
                                     xScale=xScale, yScale=yScale)
    pngTable["Average", nm] = pngFileLink
    
    ## plot the difference quality heatmap for R1_1
    colorKeyFile = paste0("diffReadsQuality-Key_", nm, ".png")
    at=seq(from=minDiff, to=maxDiff, by=by.label)
    plotCmd = expression({
      ezColorLegend(colorRange=c(minDiff, maxDiff), colors=getBlueRedScale(),
                    vertical=FALSE, by.label=by.label, at=at, labels=as.character(at))
    })
    colorKeyLink = ezImageFileLink(plotCmd, file=colorKeyFile, addPdfLink=FALSE, width=400*yScale, height=200*xScale)
    pngTable["Diff Qual Colors", nm] = colorKeyLink
    
    for(sampleName in names(qualMatrixList[idx])){
      qm = qualMatrixList[[sampleName]]
      diffResult = signif(prop.table(qm,2)*100, digits=3) - avgQual[1:nrow(qm), 1:ncol(qm)]
      pngFileLink = plotQualityHeatmap(diffResult, colorRange=c(minDiff, maxDiff), 
                                       pngFileName=paste0("diffReadsQuality-heatmap_", sampleName, ".png"),
                                       colors=getBlueRedScale(), main=paste("diffReadsQuality", sampleName, sep="_"),
                                       xScale=xScale, yScale=yScale)
      pngTable[sampleName, nm] = pngFileLink
    }
  }
  if(any(isR2)){
    pngTable = data.frame("R1"=na.omit(pngTable[,"R1"]), "R2"=na.omit(pngTable[,"R2"]), stringsAsFactors=FALSE)
  }
  
  return(pngTable)
}

##' @describeIn ezMethodFastQC Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.
plotQualityHeatmap = function(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1){
  ## some ugly controls of labels
  labCol = seq(0, ncol(result), by=10)
  labCol[1] = 1
  labRow = seq(0, nrow(result)-1, by=5)
  
  result[result > colorRange[2]] = colorRange[2]
  result[result < colorRange[1]] = colorRange[1]
  
  plotCmd = expression({
    image(1:ncol(result), 1:nrow(result), t(result), zlim=colorRange, col=colors, xlab="Read Position", ylab="Read Quality", main=main, axes=FALSE)
    axis(1, at=labCol, labels=labCol)
    axis(2, at=seq(1, nrow(result),by=5), labels=labRow)
  })
  pngFileLink = ezImageFileLink(plotCmd, file=pngFileName, height=500*yScale, width=700*xScale)
  return(pngFileLink)
}
