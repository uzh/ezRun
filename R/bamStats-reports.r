###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Plots the BAM statistics
##' @description Plots the BAM statistics and writes the report.
##' @param resultList a list of results.
##' @template dataset-template
##' @param param a list of parameters:
##' \itemize{
##'   \item{name}{ a character representing the name of the app.}
##'   \item{writeIgvSessionLink}{ a logical indicating whether to write an IGV session link.}
##'   \item{ezRef@@refBuild}{ a character containing the path of the reference build.}
##'   \item{sigThresh}{ an integer specifying the significance threshold.}
##'   \item{normMethod}{ a character specifying the normalization method to use.}
##'   \item{runGO}{ a logical indicating wheter to run the GO analysis.}
##' }
##' @template htmlFile-template
##' @template roxygen-template
plotBamStat = function(resultList, dataset, param, htmlFile=NULL){
  conds = ezConditionsFromDataset(dataset, param=param)
  samples = rownames(dataset)
  sampleColors = getSampleColors(conds, samples)
  files = dataset$BAM
  
  titles = list()
  titles[["BAM Statistics"]] = paste("BAM Statistics:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]], dataset=dataset)
  
## TODO: add later  
#   if (param$writeIgvSessionLink){
#     titles[["Genome Browser"]] = "Genome Browser"
#     addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
#     if (length(files) > 4){
#       idx = which(!duplicated(conds))
#       idx = idx[1:min(4, length(idx))]
#     } else {
#       idx = 1:length(files)
#     }
#     addIgvSessionLink(getIgvGenome(param), refBuild=param$ezRef["refBuild"], files[idx], doc, label="Open Integrative Genomics Viewer")
#   }
  
  titles[["Read Alignment Statistics"]] = "Read Alignment Statistics"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  titles[["Multi-Matching Reported in Bam File"]] = "Multi-Matching Reported in Bam File"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  doc = addParagraph(doc, "The table holds for each sample in column X the number of reads in Millions
                      that have X matches in the target and are reported in the file.")
  
  mmValues = integer()
  for (sm in samples){
    mmValues = union(mmValues, as.integer(names(resultList[[sm]]$multiMatchInFileTable)))
  }
  mmCounts = ezMatrix(0, rows=samples, cols=sort(mmValues))
  for (sm in samples){
    mm = resultList[[sm]]$multiMatchInFileTable
    mmCounts[sm, names(mm)] = mm
  }
  
  pngFile = "multiMatchInFile-barplot.png"
  pngLinks = makeAlignmentCountBarPlot(pngFile, mmCounts)
  doc = addParagraph(doc, pngLinks[["Counts"]])
  
  advancedTitles = list()
  advancedTitles[["Advanced Plots"]] = "Advanced Plots"
  advancedDoc = openBsdocReport(title=advancedTitles[[length(advancedTitles)]])
  advancedDoc = addParagraph(advancedDoc, pngLinks[["Relative"]])
  
  txtFile = "read-alignment-statistics.txt"
  colnames(mmCounts) = paste("#hits: ", colnames(mmCounts))
  colnames(mmCounts)[ncol(mmCounts)] = paste0(colnames(mmCounts)[ncol(mmCounts)], "+")
  ezWrite.table(mmCounts, file=txtFile, head="Sample")
  addTxtLinksToReport(doc, txtFile)
  
  for (nm in c("multiMatchTargetTypeCounts", "uniqueMatchTargetTypeCounts")){
    if (!is.null(resultList[[1]][[nm]])){
      readSet = switch(nm,
                       multiMatchTargetTypeCounts="Uniquely and multi-matching reads:",
                       uniqueMatchTargetTypeCounts="Uniquely matching reads:")
      titles[[paste(readSet, "Match Count Percentages")]] = paste(readSet, "Match Count Percentages")
      addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
      tct = getTypeCountTable(resultList, nm)
      ezWrite.table(tct, file=paste0(nm, ".txt"), digits=4)
      tpt = as.matrix(tct)
      for (cn in colnames(tpt)){
        tpt[ ,cn] = tct[ ,cn]/ tct["total", cn] * 100
      }
      minPercentage = 1
      rowsUse = setdiff(rownames(tpt)[apply(tpt, 1, max) > minPercentage], "total")
      tptUse = tpt[rowsUse, , drop=FALSE]
      if (nrow(tptUse) >= 2 && ncol(tptUse) >= 2){
        tptUseRel = log2(tptUse)
        tptUseRel = tptUseRel - rowMeans(tptUseRel)
        pngFile = paste0("typePercentage-", nm, "-heatmap.png")
        plotCmd = expression({
          ezHeatmap(tptUseRel, margins=c(10, 12), lim=c(-2, 2),
                    Rowv=FALSE, Colv=FALSE, main="Relative Prevalence [log2]")
        })
        heatmapLink = ezImageFileLink(plotCmd, file=pngFile, height=600, width=800)
      }
      
      doc = addFlexTable(doc, ezGrid(cbind(ifelse(nrow(tptUse) >= 2 && ncol(tptUse) >= 2, heatmapLink, NULL),
                                            as.html(ezFlexTable(signif(tptUse, digits=3), talign="right",
                                                                header.columns=TRUE, add.rownames=TRUE)))))
      
      titles[[paste(readSet, "Read Starts per Base")]] = paste(readSet, "Read Starts per Base")
      addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
      doc = addParagraph(doc, "Read Starts per Base is equivalent to Coverage divided by Read length.")
      tct = as.matrix(getTypeCoverageTable(resultList, nm))
      ezWrite.table(tct, file=paste0(nm, "-coverage.txt"), digits=4)
      if (nrow(tct) >= 2 && ncol(tct) >= 2){
        tctRel = log2(sweep(tct, 2, tct["total", ], FUN="/"))
        tctRel = tctRel[rowsUse, , drop=FALSE]
        pngFile = paste0("typeCoverage-", nm, "-heatmap.png")
        plotCmd = expression({
          ezHeatmap(tctRel, margins=c(10, 12),
                    lim=c(-5, 5),
                    Rowv=FALSE, Colv=FALSE, main="Coverage Enrichment")
        })
        heatmapLink = ezImageFileLink(plotCmd, file=pngFile, height=600, width=800)
      }
      doc = addFlexTable(doc, ezGrid(cbind(ifelse(nrow(tct) >= 2 && ncol(tct) >= 2, heatmapLink, NULL),
                                       as.html(ezFlexTable(signif(tct[rowsUse, ], digits=4), talign="right",
                                                           header.columns=TRUE, add.rownames=TRUE)))))
    }
  }
  ## use these two to record the generated png files
  plotBySamples = list()
  plotByStatistics = list()
  
  if (!is.null(resultList[[1]]$fragSizeHist)){
    pngFiles = character()
    for (sm in samples){
      fsh = resultList[[sm]]$fragSizeHist
      pngFiles[sm] = paste0(sm, "-fragSizeHist.png")
      ezWrite.table(cbind(Length=fsh$mids, Count=fsh$counts),
                     file=paste0(sm, "-fragSizeHist.txt"), row.names=FALSE)
      
      plotCmd = expression({
        plot(fsh, xlab="fragment size", main=paste(sm, "-- Length Histogram"), ylim=c(0, max(fsh$counts[-length(fsh$counts)]))) ## don't use the longest fragment size
      })
      try({
        pngFiles[sm] = ezImageFileLink(plotCmd, file=pngFiles[sm], width=600)
      })
      plotBySamples[[sm]] = c(plotBySamples[[sm]], pngFiles[sm])
      plotByStatistics[["Length distribution of fragments for paired reads"]] = c(plotByStatistics[["Length distribution of fragments for paired reads"]], pngFiles[sm])
    }
  }
  
  pngFiles = character()
  for (sm in samples){
    fsh = resultList[[sm]]$segmentCountHist
    pngFiles[sm] = paste0(sm, "-segmentCountHist.png")
    ezWrite.table(cbind(Length=fsh$mids, Count=fsh$counts),
                   file=paste0(sm, "-segmentCountHist.txt"), row.names=FALSE)
    
    plotCmd = expression({
      plot(fsh, xlab="# segments in alignment", main=paste(sm, "-- Histogram of Segments per Alignment"))
    })
    try({
      pngFiles[sm] = ezImageFileLink(plotCmd, file=pngFiles[sm], width=600)
    })
    plotBySamples[[sm]] = c(plotBySamples[[sm]], pngFiles[sm])
    plotByStatistics[["Histogram of aligned segments per read"]] = c(plotByStatistics[["Histogram of aligned segments per read"]], pngFiles[sm])
  }
  
  if (!is.null(resultList[[1]][["ErrorRates"]])){
    for (sm in samples){
      pngLinks = character()
      for (nm in names(resultList[[sm]][["ErrorRates"]])){
        errorRate = resultList[[sm]][["ErrorRates"]][[nm]]
        if (!is.null(errorRate)){
          pngFile = ezValidFilename(paste0(sm, "_", nm, ".png"))
          plotCmd = expression({
            plotPosSpecificErrorRate(errorRate, png=pngFile, main=paste(sm, nm))
          })
          pngLinks[nm] = ezImageFileLink(plotCmd, file=pngFile, width=1600)
        }
      }
      plotBySamples[[sm]] = c(plotBySamples[[sm]], pngLinks)
      plotByStatistics[["Read position specific error rate"]] = c(plotByStatistics[["Read position specific error rate"]], pngLinks)
    }
  }
  
  if (!is.null(resultList[[1]][["TranscriptsCovered"]])){
    pngFile = "TranscriptsCovered.png"
    minYlim = 0 #min(sapply(resultList, function(item){min(item[["TranscriptsCovered"]][["counts"]])}))
    minXlim = 0
    maxYlim = max(sapply(resultList, function(item){max(item[["TranscriptsCovered"]][["counts"]])}))
    maxXlim = 130 #max(sapply(resultList, function(item){max(item[["TranscriptsCovered"]][["mids"]])}))
    
    plotCmd = expression({
      plot(1, 1, xlim=c(minXlim, maxXlim), ylim=c(minYlim, maxYlim), xlab = "% length covered", ylab="number of isoforms", main="Isoforms Covered Histogram", type="n")
      for (sm in samples){
        transcript_covered = resultList[[sm]][["TranscriptsCovered"]]
        lines(transcript_covered[["mids"]], transcript_covered[["counts"]], col=sampleColors[sm])
      }
      legend("topright", samples, col=sampleColors[samples], cex=1.2, pt.cex=1.5, bty="o", pt.bg="white", lty=1)
    })
    pngLink = ezImageFileLink(plotCmd, file=pngFile, width=700)
    
    titles[["Transcripts covered plot"]] = "Transcripts covered plot"
    addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
    doc = addParagraph(doc, pngLink)
    
    pngFiles = character()
    pngLinks = character()
    for (sm in samples){
      tlc = resultList[[sm]][["TranscriptsCovered"]]
      pngFiles[sm] = paste0(sm, "-transcriptsCovered.png")
      ezWrite.table(cbind(Percents=tlc$mids, Count=tlc$counts),
                     file=paste0(sm, "-transcriptsCovered.txt"), row.names=FALSE)
      
      plotCmd = expression({
        cts = tlc$counts
        names(cts) = tlc$mids
        barplot(cts, xlab="% length covered", main=paste(sm, "-- Isoforms Covered Histogram"), ylab="number of isoforms")
      })
      pngLinks[sm] = ezImageFileLink(plotCmd, file=pngFiles[sm], width=700)
      plotBySamples[[sm]] = c(plotBySamples[[sm]], pngLinks[sm])
      plotByStatistics[["The fraction of isoform length covered"]] = c(plotByStatistics[["The fraction of isoform length covered"]], pngLinks[sm])
    }
  }
  
  if (length(resultList[[1]][["genebody_coverage"]]) != 0){ ## TODO this could be done better by searching alls results for a valide genebody_coverage element
    minYlim = 0 #min(sapply(resultList, function(item){min(item[["genebody_coverage"]])}))
    maxYlim = 0.08 ## this means we allow at most 10-fold enrichment at a percentile #max(sapply(resultList, function(item){max(item[["genebody_coverage"]])}))
    gbcTemplate = resultList[[1]][["genebody_coverage"]]
    pngMatrix = ezMatrix("", rows=names(gbcTemplate), cols=names(gbcTemplate[[1]]))
    pngLinks = character()
    for (rn in rownames(pngMatrix)){
      for (cn in colnames(pngMatrix)){
        ## skip all cases that are not medium expressed and skip all cases that are not above 4000 or 400-1000
        if (!grepl("medium", cn)) next
        if (grepl("less", rn) | grepl("1000 to", rn)) next
        pngMatrix[rn, cn] = ezValidFilename(paste0("genebody_coverage_", rn, "_", cn, ".png"))
        covValues = ezMatrix(0, cols=0:100, rows=samples)
        for (sm in samples){
          y = resultList[[sm]][["genebody_coverage"]][[rn]][[cn]]
          if (!is.null(y)){
            covValues[sm, ] = y
          } else {
            covValues[sm, ] = NA
          }
        }
        plotCmd = expression({
          plot(1, 1, xlim=c(0,100), ylim=c(minYlim, maxYlim), xlab="percentile of geneBody (5'->3')", ylab="relative coverage", 
               main=paste("Genebody coverage", rn, cn), type="n",
               axes=FALSE, frame=TRUE)
          axis(side=2)
          axis(side=1, at=seq(0, 100, by=10))
          for (sm in samples){
            y = resultList[[sm]][["genebody_coverage"]][[rn]][[cn]]
            if (!is.null(y)){
              covValues[sm, ] = y
              lines(0:100, y, col=sampleColors[sm])
            } else {
              covValues[sm, ] = NA
            }
          }
          #legend("topright", samples, col=sampleColors[samples], cex=1.2, pt.cex=1.5, bty="o", pt.bg="white", lty=1)
        })
        link = ezImageFileLink(plotCmd, file=pngMatrix[rn, cn], width=600)
        pngLinks = append(pngLinks, link)
        # pngLinks[rn, cn] = ezImageFileLink(plotCmd, file=pngMatrix[rn, cn], width=600)
        ezWrite.table(covValues, file=sub(".png$", ".txt", pngMatrix[rn, cn]), head="Name")
      }
    }
    titles[["Genebody coverage plot"]] = "Genebody coverage plot"
    addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
    
    plotCmd = expression({
      par(mar=c(0,0,0,0))
      plot(1,1, axes=FALSE, frame=FALSE, type="n", xlab="", ylab="")
      legend("topleft", samples, fill=sampleColors[samples], border=NA, bty="n", pt.bg="white", title="Sample Colors")
    })
    sampleLink = ezImageFileLink(plotCmd, file="sampleColors.png", height=length(samples)*15 + 20, width=300)
    
    doc = addParagraph(doc, sampleLink)
    doc = addFlexTable(doc, ezGrid(rbind(pngLinks)))
  }
  
  geneCounts = resultList[[1]][["geneCounts"]]
  if (!is.null(geneCounts)){
    counts = ezMatrix(0, rows=names(geneCounts), cols=names(resultList))
    for (i in 1:length(resultList)){
      counts[ ,i] = resultList[[i]][["geneCounts"]][rownames(counts)]
    }
    seqAnno = ezFeatureAnnotation(param, rownames(counts), "gene")
    rawData = list(counts=counts, isLog=FALSE,
                   presentFlag=counts > param$sigThresh, 
                   seqAnno=seqAnno, featureLevel="gene",
                   type="Counts", countName="multiMatchCounts")
    if (is.null(param$normMethod)){
      param$normMethod = "logMean"
    }
    if (!ezIsSpecified(param$runGO)){
      param$runGO = TRUE
    }
    setwdNew("Count_QC")
    rawData$signal = ezNorm(rawData$counts, presentFlag=rawData$presentFlag, method=param$normMethod)
    runNgsCountQC(dataset, "00index.html", param, rawData=rawData)
    setwd("..")
    doc = addParagraph(doc, pot("Count QC Report", hyperlink="Count_QC/00index.html"))
  }
  
  if (!is.null(names(resultList[[1]][["Junction"]]))){
    junctionMaxVal = numeric()
    for (nm in names(resultList[[1]][["Junction"]])){
      junctionMaxVal[nm] = 0
      for (sm in samples){
        junctionMaxVal[nm] = max(junctionMaxVal[nm], unlist(resultList[[sm]][["Junction"]][[nm]]))
      }
    }
    for (sm in samples){
      pngLinks = character()
      for (nm in names(resultList[[sm]][["Junction"]])){
        junctionPlot = resultList[[sm]][["Junction"]][[nm]]
        pngFile = ezValidFilename(paste0(sm, "_", nm, ".png"))
        
        plotCmd = expression({
          if (nm %in% c("splice_events", "splice_junction")){
            pie(junctionPlot, col=c(2,3,4), init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=nm, labels=paste(names(junctionPlot), paste0(round(junctionPlot), "%")))
          } else if (nm =="junctionSaturation"){
            x = as.numeric(names(junctionPlot[[1]])) * 100
            plot(1,1,xlab="percent of total reads", ylab='Number of splicing junctions (x1000)',type='o',
                 ylim=c(0, junctionMaxVal[nm]/1000), xlim=range(x))
            saturationColors = c("all junctions"="blue", "known junctions"="red", "novel junctions"="green")
            for (item in names(junctionPlot)){
              lines(x, junctionPlot[[item]]/1000, col=saturationColors[item], type="o")
            }
            legend("topleft", legend=names(saturationColors), col=saturationColors,lwd=1,pch=1)
          }
        })
        pngLinks[nm] = ezImageFileLink(plotCmd, file=pngFile, width=600)
      }
      plotBySamples[[sm]] = c(plotBySamples[[sm]], pngLinks)
      plotByStatistics[["Junction Plots"]] = c(plotByStatistics[["Junction Plots"]], pngLinks)
    }
    
    ## do the junction saturation plot in main page for all samples
    titles[["Junction saturation plot for all samples"]] = "Junction saturation plot for all samples"
    addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
    
    for (nm in names(resultList[[1]][["Junction"]][["junctionSaturation"]])){
      pngFile = ezValidFilename(paste0("junctionSaturation_", nm, ".png"))
      plotCmd = expression({
        plot(1,1,xlab="percent of total reads", ylab='Number of splicing junctions (x1000)',type='o', ylim=c(0, junctionMaxVal["junctionSaturation"]/1000), xlim=c(0, 130), main=nm)
        for (sm in samples){
          lines(x, resultList[[sm]][["Junction"]][["junctionSaturation"]][[nm]]/1000, col=sampleColors[sm], type="o")
        }
        legend("bottomright", legend=samples, col=sampleColors[samples], lwd=1, pch=1)
      })
      pngLink = ezImageFileLink(plotCmd, file=pngFile, width=700)
      
      if (grepl("known", nm)){
        doc = addParagraph(doc, pngLink)
      } else {
        plotByStatistics[["Junction Plots"]] = c(plotByStatistics[["Junction Plots"]], pngLink)
      }
    }
  }
  
  ## write the resulys by statistics per page to the main pages
  titles[["The result plots by each statistics"]] = "The result plots by each statistics"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  tableOfPages = ezMatrix("", rows=names(plotByStatistics), cols="Plots")
  # replace the space with undersocre if available, otherwise some errors will occur in the html
  plotByStatisticsPages = paste0(gsub("[[:space:]]+", "_", names(plotByStatistics)), ".html")
  for (i in 1:length(plotByStatistics)){
    tableOfPages[i, ] = as.html(pot(plotByStatisticsPages[i], hyperlink=plotByStatisticsPages[i]))
    subTitles = list()
    subTitles[[plotByStatisticsPages[i]]] = plotByStatisticsPages[i]
    subDoc = openBsdocReport(title=subTitles[[length(subTitles)]])
    subTitles[[names(plotByStatistics)[i]]] = names(plotByStatistics)[i]
    addTitle(subDoc, subTitles[[length(subTitles)]], 3, id=subTitles[[length(subTitles)]])
    subDoc = addFlexTable(subDoc, ezGrid(plotByStatistics[[i]]))
    closeBsdocReport(subDoc, plotByStatisticsPages[i], subTitles)
  }
  doc = addFlexTable(doc, ezFlexTable(tableOfPages, header.columns=TRUE, add.rownames=TRUE))
  
  ## write the reults by sample per page to the main pages
  titles[["The result plots by each sample"]] = "The result plots by each sample"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  plotBySamplesPages = paste0(gsub("[[:space:]]+", "_", names(plotBySamples)), ".html")
  tableOfPages = ezMatrix("", rows=names(plotBySamples), cols="Plots")
  for (i in 1:length(plotBySamples)){
    tableOfPages[i, ] = as.html(pot(plotBySamplesPages[i], hyperlink=plotBySamplesPages[i]))
    subTitles = list()
    subTitles[[plotBySamplesPages[i]]] = plotBySamplesPages[i]
    subDoc = openBsdocReport(title=subTitles[[length(subTitles)]])
    subTitles[[names(plotBySamples)[i]]] = names(plotBySamples)[i]
    addTitle(subDoc, subTitles[[length(subTitles)]], 3, id=subTitles[[length(subTitles)]])
    subDoc = addFlexTable(subDoc, ezGrid(plotBySamples[[i]]))
    closeBsdocReport(subDoc, plotBySamplesPages[i], subTitles)
  }
  doc = addFlexTable(doc, ezFlexTable(tableOfPages, header.columns=TRUE, add.rownames=TRUE))
  closeBsdocReport(advancedDoc, "advancedPlots.html", advancedTitles)
  doc = addParagraph(doc, pot("advancedPlots", hyperlink = "advancedPlots.html"))
  closeBsdocReport(doc, htmlFile, titles)
}

##' @describeIn plotBamStat Gets the table containing the counts of each type.
getTypeCountTable = function(resultList, name){
  tbl = data.frame(row.names=rownames(resultList[[1]][[name]]))
  for (sm in names(resultList)){
    counts = resultList[[sm]][[name]]
    tbl[sm] = signif(counts[ rownames(tbl), "count"] / 1e6, digits=4)
  }
  return(tbl)
}

##' @describeIn plotBamStat Gets the table containing the coverage of each type.
getTypeCoverageTable = function(resultList, name){
  tbl = data.frame(row.names=rownames(resultList[[1]][[name]]))
  for (sm in names(resultList)){
    counts = resultList[[sm]][[name]]
    tbl[sm] = counts[ rownames(tbl), "count"] / counts[ rownames(tbl), "width"]
  }
  return(tbl)
}

##' @describeIn plotBamStat Plots the alignment counts and returns the image file link.
makeAlignmentCountBarPlot = function(file, mmCounts){
  multiCount = as.integer(colnames(mmCounts))
  isSmall = multiCount <= 3
  if (any(!isSmall)){
    mmCounts = cbind(mmCounts[ , isSmall, drop=FALSE], 
                     ">3"=rowSums(mmCounts[ , !isSmall, drop=FALSE]))
  }
  ezWrite.table(mmCounts, file=sub(".png", ".txt", file))
  multiCountColors = c("0 hit(s)"="gray", "1 hit(s)"="blue", "2 hit(s)"="cyan",
                       "3 hit(s)"="green", ">3 hit(s)"="orange")
  colnames(mmCounts) = paste(colnames(mmCounts), "hit(s)")
  stopifnot(colnames(mmCounts) %in% names(multiCountColors))
  pngLinks = character()
  plotCmd = expression({
    par(mar=c(12, 4.1, 4.1, 2.1))
    x = mmCounts[ , rev(colnames(mmCounts))]
    barplot(t(x)/1e6, las=2, ylab="Counts [Mio]", main="total alignments", legend.text=TRUE, border=NA,
            col=multiCountColors[colnames(x)], xlim=c(0, nrow(x) +5))
  })
  pngLinks["Counts"] = ezImageFileLink(plotCmd, file=file, width=400 + nrow(mmCounts) * 10, height=480)
  
  plotCmd = expression({
    par(mar=c(12, 4.1, 4.1, 2.1))
    x = mmCounts[ , rev(colnames(mmCounts))]
    for (i in 1:nrow(x)) {
      x[i, ] = x[i, ]/sum(x[i,])
    }
    barplot(t(x), las=2, ylab="Counts [proportion]", main="alignment proportions", legend.text=TRUE, border=NA,
            col=multiCountColors[colnames(x)], xlim=c(0, nrow(x) +5))
  })
  pngLinks["Relative"] = ezImageFileLink(plotCmd, file="multiMatchInFile-barplot-rel.png", width=400 + nrow(mmCounts) * 10, height=480)
  
  return(pngLinks)
}

##' @describeIn plotBamStat Plots the position specific error rates.
plotPosSpecificErrorRate = function(errorRate, png, main="Per base mismatch rate", writeTxt=TRUE){
  par(mfrow=c(1,2))
  if (writeTxt){
    for (i in 1:length(errorRate)){
      ezWrite.table(errorRate[[i]], sub(".png$", paste0("-",names(errorRate)[i], ".txt"), png))
    }
  }
  trimmedClippedRateMatrix = rbind(errorRate[["trimmedRate"]], errorRate[["clippedRate"]])
  errorRateMatrix = errorRate[["errorRate"]] 
  colnames(trimmedClippedRateMatrix) = 1:ncol(trimmedClippedRateMatrix)
  barplot(trimmedClippedRateMatrix, col=c("darkorange", "blue"), space=0, border=NA,
          xlab="Base Position in Read", ylab="Trimmed-Clipped Rate", 
          main=paste(main, "trimmed-clipped rate"), ylim=c(0, 1), #max(0.1, 1.2*max(colSums(trimmedClippedRateMatrix)))), 
          xpd=FALSE)
  legend("top", c("trimmed", "clipped"), col=c("darkorange", "blue"), pch=20)
  yMax = max(0.1, max(errorRateMatrix, na.rm=TRUE))
  barCol = ifelse (yMax > 0.1, "red", "gray")
  barplot(errorRateMatrix, space=0, border=NA,
          xlab="Base Position in Read", ylab="Mismatch Rate",  
          main=paste(main, "mismatch rate"), ylim=c(0, yMax), col=barCol,
          xpd=FALSE)
  if (yMax > 0.1){
    warning("The error rate range is larger than 0.1")
  }
}
