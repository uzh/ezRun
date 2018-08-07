###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

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
    x = mmCounts[ , rev(colnames(mmCounts)), drop = F]
    barplot(t(x)/1e6, las=2, ylab="Counts [Mio]", main="total alignments", legend.text=TRUE, border=NA,
            col=multiCountColors[colnames(x)], xlim=c(0, nrow(x) +5),
            names.arg=ezSplitLongLabels(rownames(x)))
  })
  pngLinks["Counts"] = ezImageFileLink(plotCmd, file=file, width=min(600 + (nrow(mmCounts)-10)* 30, 2000)) # nSamples dependent width
  
  plotCmd = expression({
    par(mar=c(12, 4.1, 4.1, 2.1))
    x = mmCounts[ , rev(colnames(mmCounts)), drop = F]
    #for (i in 1:nrow(x)) {
    #  x[i, ] = x[i, ]/sum(x[i,])
    #}
    x <- sweep(x, MARGIN=1, STATS=rowSums(x), FUN="/")
    barplot(t(x), las=2, ylab="Counts [proportion]", main="alignment proportions", legend.text=TRUE, border=NA,
            col=multiCountColors[colnames(x)], xlim=c(0, nrow(x) +5),
            names.arg=ezSplitLongLabels(rownames(x)))
  })
  pngLinks["Relative"] = ezImageFileLink(plotCmd, file="multiMatchInFile-barplot-rel.png", width=min(600 + (nrow(mmCounts)-10)* 30, 2000)) # nSamples dependent width
  
  return(pngLinks)
}

alignmentCountBarPlot <- function(mmCounts, relative=FALSE,
                                  file=NULL){
  require(plotly)
  title <- ifelse(relative, "alignment proportions", "total alignments")
  multiCount = as.integer(colnames(mmCounts))
  isSmall = multiCount <= 3
  if (any(!isSmall)){
    mmCounts = cbind(mmCounts[ , isSmall, drop=FALSE], 
                     ">3"=rowSums(mmCounts[ , !isSmall, drop=FALSE]))
  }
  
  if(!is.null(file)){
    ezWrite.table(mmCounts, file=file)
  }
  
  multiCountColors = c("0 hit(s)"="gray", "1 hit(s)"="blue", "2 hit(s)"="cyan",
                       "3 hit(s)"="green", ">3 hit(s)"="orange")
  colnames(mmCounts) = paste(colnames(mmCounts), "hit(s)")
  stopifnot(colnames(mmCounts) %in% names(multiCountColors))
  
  x = mmCounts[ , rev(colnames(mmCounts)), drop = F]
  if(isTRUE(relative)){
    x <- sweep(x, MARGIN=1, STATS=rowSums(x), FUN="/")
  }
  x <- as.data.frame(x)
  x$"sample" <- rownames(x)
  m <- list(
    l = 80,
    r = 80,
    b = 200,
    t = 100,
    pad = 0
  )
  
  p <- plot_ly(x, x=~sample, y=as.formula(paste0("~`", colnames(x)[1], "`")), 
               type="bar", name = colnames(x)[1], 
               marker = list(color = multiCountColors[colnames(x)[1]]))
  for(i in 2:(ncol(x)-1)){
    p <- p %>% add_trace(y=as.formula(paste0("~`", colnames(x)[i], "`")), 
                         name=colnames(x)[i],
                         marker = list(color = multiCountColors[colnames(x)[i]]))
  }
  p <- p %>% plotly::layout(yaxis = list(title = 'Count'), barmode = 'stack',
                            title = title, margin=m)
  p
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
