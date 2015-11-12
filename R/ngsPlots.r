###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## REFAC, but function is currently unused.
fragLengthReadMidDensityPlot = function(fl, file=NULL, xlim=c(-500, 500), bw=5){
  
  if (!is.null(file)){
    pdf(file=file)
    on.exit(dev.off())
  }
  dens1 = density(fl$dist1, bw=bw)
  dens2 = density(fl$dist2, bw=bw)
  yMax = max(dens1$y, dens2$y)
  plot(dens1$x, dens1$y, col="black", ylim=c(0, yMax), xlim=xlim,
       main=paste(sm, "shift:", signif(fl$fragSizeMean, digits=4),
                  "+-", signif(fl$fragSizeSd, digits=4)),
       xlab="Distance to TSS", type="l")
  lines(dens2, col="lightblue")
  lines(dens2$x-fl$fragSizeMean + fl$readSizeMean, dens2$y, col="blue")
  legend("bottomright", c("Forward Reads", "Reverse Reads", "Shifted Reverse Reads"),
                          col=c("black", "lightblue", "blue"), pch=20)
}


## REFAC, but function is currently unused.
fragLengthCovPlot = function(fl, avgFragProf, file=NULL, xlim=c(-500:499), ylim=c(0,max(avgFragProf$forw, avgFragProf$rev)),...){
  
  if (!is.null(file)){
    pdf(file=file)
    on.exit(dev.off())
  }
  plot(-500:499, avgFragProf$forw, xlab="Distance to TSS", ylab="Avg coverage", type="l", ylim=ylim,
       col="black", ...)
  lines(-500:499, avgFragProf$rev,
       col="blue")
  abline(h=0)
}


## REFAC, but function is currently unused.
plotStrandedProfiles = function(gr, covList, gff=gff,
                              featColors = getFeatColors(unique(gff$type)),
                              sampleColors=getSampleColors(names(covList), names(covList)),
                              extend=1e3, minHeight=100, scalingFactor=NULL,
                              file=NULL){
  
  if (!is.null(file)){
      png(file=file, height=1000, width=1600)
      on.exit(dev.off())
  }
  
  if (is.null(scalingFactor)){
    scalingFactor = rep(1, length(covList))
    names(scalingFactor) = names(covList)
  }
  
  name = names(gr)[1]
  if (class(gr) == "GRangesList"){
    gr = unlist(gr)
  }
  #stopifnot(length(gr) == 1)
  start = min(start(gr))
  end = max(end(gr))
  chr = as.character(seqnames(gr)[1])
  xlim = c(max(1, start - extend), min(end + extend, covList[[1]]$seqLengths[chr], na.rm=TRUE))
  if (!chr %in% names(covList[[1]][[1]])){
    #### we have no coverage for this chromosome
    par(mfrow=c(5,1), cex=1.2)
    plotGenomicFeatures(chr, xlim[1], xlim[2], gff, featColors=featColors,  
                                   xlim=xlim,
                                   showLegend=TRUE,
                                   main=paste(chr, name, "-- no read coverage for this chrom"))
    return()
  }
  

  ### get the profiles
  modes = character()
  hasUnique = all(sapply(covList, function(x){!is.null(x[["cov unique plus"]])}))
  if (hasUnique){
    modes = c(modes, "unique")
  }
  hasMulti = all(sapply(covList, function(x){!is.null(x[["cov multi plus"]])}))
  if (hasMulti){
    modes = c(modes, "multi")
  }
  profiles = list()
  for (mode in modes){
    profiles[[mode]] = list()
    for (strand in c("plus", "minus")){
      profiles[[mode]][[strand]] = list()
      for (sm in names(covList)){
        nm = paste("cov", mode, strand)
        profiles[[mode]][[strand]][[sm]] = covList[[sm]][[nm]][[chr]][xlim[1]:xlim[2]] * scalingFactor[sm]
      }
    }
  }
        
  #ylimUnique = range(sapply(unlist(profiles[["unique"]]), range))
  #ylimUnique[2] = shrinkToRange(ylimUnique[2], c(50, 1000))
  #ylimMulti = range(sapply(unlist(profiles[["multi"]]), range))
  #ylimMulti[2] = shrinkToRange(ylimMulti[2], c(50, 1000))
  if (hasUnique){
    ylim = c(0, max(sapply(unlist(profiles[["unique"]]), function(x){max(x[xlim[1]:xlim[2] %in% start:end])})))
  } else {
    ylim = c(0, max(sapply(unlist(profiles[["multi"]]), function(x){max(x[xlim[1]:xlim[2] %in% start:end])})))
  }
  ylim[2] = ylim[2] * 1.3
  ylim[2] = shrinkToRange(ylim[2], c(50, 1000))
  

  par(mfrow=c(1 + 2 * length(modes),1), cex=1.2)
  par(mar=c(0, 4.1, 3, 3))
  if (hasMulti){
    profilePlot(x=xlim[1]:xlim[2], profiles[["multi"]][["plus"]], sampleColors[names(covList)], ylim=ylim, clampValues=TRUE, showLegend=TRUE,
                ylab="multi matches -- plus", main=paste(chr, name))
    par(mar=c(0, 4.1, 0.5, 3))
  }
  if (hasUnique){
    profilePlot(x=xlim[1]:xlim[2], profiles[["unique"]][["plus"]], sampleColors[names(covList)], ylim=ylim, clampValues=TRUE, showLegend=FALSE,
                ylab="unique matches -- plus")
    par(mar=c(0, 4.1, 0.5, 3))
  }
  plotGenomicFeatures(chr, xlim[1], xlim[2], gff, featColors=featColors,  
                                 xlim=xlim,
                                 showLegend=FALSE)
  if (hasUnique){
    profilePlot(x=xlim[1]:xlim[2], profiles[["unique"]][["minus"]], sampleColors[names(covList)], ylim=ylim, clampValues=TRUE, showLegend=FALSE,
                ylab="unique matches -- minus")
  }
  if (hasMulti){
    profilePlot(x=xlim[1]:xlim[2], profiles[["multi"]][["minus"]], sampleColors[names(covList)], ylim=ylim, clampValues=TRUE, showLegend=FALSE,
                ylab="multi matches -- minus")
  }
}



## @describeIn plotStrandedProfiles
profilePlot = function(x=1:length(profileList[[1]]), profileList,
                       profColors, legendLabels=names(profColors), file=NULL, 
                       ylim = range(unlist(profileList)),
                       xlim = range(x), addOn=NULL, clampValues=FALSE,
                       lwd=rep(2, length(profileList)), lty=NULL, smoothWindow=1, loessSpan=0,
                       xlab="Position", ylab="Coverage", main="", cex=1,
                       height=7, width=7, showLegend=TRUE, logVal=""){

  if (!is.null(file)){
    pdf(file, height=height, width=width )
    on.exit(dev.off())
  }
  if (is.null(lty)){
    lty = rep(1, length(profileList))
    names(lty) = names(profileList)
  }
  for (i in 1:length(profileList)){
    y = profileList[[i]]
    if (sum(!is.na(y))>0){
      if (smoothWindow >1){
        y = runmed(y, smoothWindow)
      }
      if (loessSpan > 0){
        y = predict(loess(y ~ x, span=loessSpan))
      }
    }
    profileList[[i]] = y
  }
      ## alternative smooth by lowess: ys = predict(loess(avgPelegBcs ~ x, span=0.3))
  plot(-1, -1, xlim=xlim, ylim=ylim,
       xlab=xlab, ylab=ylab,
       main=main, col="white", log=logVal)
  for (i in 1:length(profileList)){
    if (clampValues){
      y = shrinkToRange(as.numeric(profileList[[i]]), ylim)
    } else {
      y = profileList[[i]]
    }
    lines(x, y, lwd=lwd[i], lty=lty[names(profileList)[i]], col=profColors[names(profileList)[i]])
  }
  if (showLegend){
    legend("topright", legendLabels,
           col=profColors, pch=20, cex=cex) #, cex=0.9, pch=20, bty="o", pt.bg="white")
  }
  if (!is.null(addOn)){
    addOn
  }
}


## REFAC, but function is currently unused.
cumHistPlot = function(param, cts, png, colors, main="all transcripts"){

  if (!is.null(png)){
    png(filename=png, width=640, height=640)
    on.exit(dev.off())
  }
  useCts = cts > 0
  xlim = c(1, max(cts))
  ylim = c(1, max(colSums(useCts)))
  plot(1, 1, xlim=xlim, ylim=ylim, log="x", pch=16, col="white", 
       xlab="digital expression", ylab="number of transcripts",
       main=paste(main, "(",nrow(cts), ")"))
  for (sm in colnames(cts)){
    val = sort(cts[useCts[, sm], sm])
    nVal = length(val)
    lines(val, 1:nVal, col=colors[sm])
    points(val[(nVal-9):nVal], (nVal-9):nVal, pch=16, col=colors[sm])
  }
  legend("bottomright", colnames(cts), col=colors[colnames(cts)], cex=1.2, pt.cex=1.5, pch=20, bty="o", pt.bg="white")
  return(NULL)
}


countDensPlot = function(param, cts, colors, main="all transcripts", bw=7){
  
  cts[cts < 0] = 0  ## zeros will be come -Inf and not be part of the area!! Area will be smaller than 1!
  xlim = 0
  ylim = 0
  densList = list()
  for (sm in colnames(cts)){
    densList[[sm]] = density(log2(cts[ ,sm]), bw=bw)
    xlim = range(c(xlim, densList[[sm]]$x))
    ylim = range(c(ylim, densList[[sm]]$y))
  }
  xlim[2] = xlim[2] + 5
  plot(1, 1, xlim=xlim, ylim=ylim, pch=16, col="white", 
       xlab="log2 expression", ylab="density of transcripts",
       main=paste(main, "(",nrow(cts), ")"))
  for (sm in colnames(cts)){
    lines(densList[[sm]]$x, densList[[sm]]$y, col=colors[sm])
  }
  legend("topright", colnames(cts), col=colors[colnames(cts)], cex=1.2, pt.cex=1.5, pch=20, bty="o", pt.bg="white")
  return(NULL)
}

myMdsPlot = function(signal, sampleColors, main){
  require(edgeR)
  y = DGEList(counts=signal,group=colnames(signal))
  #y$counts = cpm(y)
  #y = calcNormFactors(y)
  mds = plotMDS(y)
  plot(mds$x, mds$y, pch=c(15), xlab='Leading logFC dim1', ylab='Leading logFC dim2', main=main,
       xlim=c(1.2*min(mds$x), 1.2*max(mds$x)), ylim=c(1.2*min(mds$y), 1.2*max(mds$y)), col=sampleColors)
  text(mds$x,mds$y,labels = colnames(signal),pos=1,col=c('darkcyan'),cex=0.7)
  par(bg = 'white')
  return(NULL)
}

