###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Plots the count densities
##' @description Plots the count densities.
##' @param cts the signal to use for plotting.
##' @template colors-template
##' @param main a character representing the plot title.
##' @param bw a numeric passed to \code{density()}.
##' @template roxygen-template
countDensPlot = function(cts, colors, main="all transcripts", bw=7){
  
  cts[cts < 0] = 0  ## zeros will become -Inf and not be part of the area!! Area will be smaller than 1!
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

countDensGGPlot <- function(cts, xmin=min(cts, na.rm=TRUE)-5, 
                            xmax=max(cts, na.rm=TRUE)+5, alpha=0.3, colors, main="all transcripts") {
  require(ggplot2, quietly = T)
  cts[cts < 0] = 0
  cts = log2(cts)
  data = data.frame(signal = unlist(cts),sampleName=rep(colnames(cts),each=nrow(cts)), stringsAsFactors = F)
  xvar = 'signal'
  split = 'sampleName'
  p = ggplot(data=data, aes_string(x=xvar, fill=split))
  p = p + xlim(xmin,xmax) + geom_density(alpha=alpha) + scale_fill_manual(values=colors)
  p = p + ylab("density of transcripts") + xlab("log2 expression") + ggtitle(paste(main, "(",nrow(cts), ")"))
  p = p + theme_bw() + theme(plot.background = element_blank(), #Remove default GRID
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank())
  p = p + theme(axis.line.x = element_line(color="black"), #Draw AxisLines
                axis.line.y = element_line(color="black"))
  p = p + theme(legend.position="bottom", legend.title = element_blank()) #Legend
  return(p)
}


##' @title Plots the multi dimensional scaling
##' @description Plots the multi dimensional scaling.
##' @param signal the signal to use for plotting.
##' @param sampleColors a character vector containing colors in hex format.
##' @param main a character representing the plot title.
##' @template roxygen-template
ezMdsPlot = function(signal, sampleColors, main){
  require("edgeR")
  y = DGEList(counts=signal,group=colnames(signal))
  #y$counts = cpm(y)
  #y = calcNormFactors(y)
  
  mds = plotMDS(y, plot=FALSE)
  
  plot(mds$x, mds$y, pch=c(15), xlab='Leading logFC dim1', ylab='Leading logFC dim2', main=main,
       xlim=c(1.2*min(mds$x), 1.2*max(mds$x)), ylim=c(1.2*min(mds$y), 1.2*max(mds$y)), col=sampleColors)
  text(mds$x,mds$y,labels = colnames(signal),pos=1,col=c('darkcyan'),cex=0.7)
  par(bg = 'white')
}

ezMdsGG2 <- function(signal, design, ndim=2, main="MDS plot"){
  require(edgeR)
  require(plotly)
  require(ggrepel)
  if(ndim != 2){
    stop("ggplot2 only produces 2D plot")
  }
  
  y = DGEList(counts=signal, group=colnames(signal))
  mds = plotMDS(y, plot=FALSE, ndim=ndim)
  toPlot <- data.frame(samples=colnames(signal),
                       design,
                       stringsAsFactors = FALSE)
  mdsOut <- mds$cmdscale.out
  colnames(mdsOut) <- c("Leading logFC dim1", "Leading logFC dim2")
  toPlot <- cbind(toPlot, mdsOut)
  if(ncol(design) > 1L){
    p <- ggplot(toPlot, aes(`Leading logFC dim1`, `Leading logFC dim2`)) +
      geom_point(aes_string(colour=colnames(design)[1],
                            shape =colnames(design)[2]),
                 size = 3) + 
      geom_text_repel(aes(label=samples)) +
      theme_bw() + ggtitle(main)
  }else{
    p <- ggplot(toPlot, aes(`Leading logFC dim1`, `Leading logFC dim2`)) +
      geom_point(aes_string(colour=colnames(design)[1]),
                 size = 3) + 
      geom_text_repel(aes(label=samples)) +
      theme_bw() + ggtitle(main)
  }
  p
}

ezMdsPlotly <- function(signal, design, ndim=c(3,2), main, sampleColors=NULL){
  require("edgeR")
  require(plotly)
  y = DGEList(counts=signal, group=colnames(signal))
  mds = plotMDS(y, plot=FALSE, ndim=ndim)
  factorToPlot = colnames(design)[1]
  toPlot <- data.frame(samples=colnames(signal),
                       design,
                       sampleColors,
                       stringsAsFactors = FALSE)
  mdsOut <- mds$cmdscale.out
  
  # names of color vector need to match the selected design factor
  pal = toPlot[,c(factorToPlot,"sampleColors")]
  pal = pal[!duplicated(pal),]
  pal = setNames(pal[['sampleColors']],pal[[factorToPlot]])
    
  if(ndim == 3){
    colnames(mdsOut) <- c("Leading logFC dim1", "Leading logFC dim2", 
                          "Leading logFC dim3")
    toPlot <- cbind(toPlot, mdsOut)
    p <- plot_ly(toPlot, x=~`Leading logFC dim1`, y=~`Leading logFC dim2`, z=~`Leading logFC dim3`,
                 color=formula(paste0("~", factorToPlot)),
                 colors=pal,
                 type='scatter3d',
                 mode='markers+text',
                 text=~samples, textposition = "top right")%>%
      plotly::layout(title=main, 
                     scene=list(xaxis=list(title = 'Leading logFC dim1'),
                                yaxis = list(title = 'Leading logFC dim2'),
                                zaxis = list(title = 'Leading logFC dim3')))
  }else if(ndim ==2){
    colnames(mdsOut) <- c("Leading logFC dim1", "Leading logFC dim2")
    toPlot <- cbind(toPlot, mdsOut)
    if(ncol(design) > 1L){
      p <- plot_ly(toPlot, x=~`Leading logFC dim1`, y=~`Leading logFC dim2`,
              color = formula(paste0("~", factorToPlot)),
              symbol=formula(paste0("~", colnames(design)[2])),
              colors = pal,
              type='scatter',
              mode='markers+text',
              text=~samples, textposition = "top right")
    }else{
      p <- plot_ly(toPlot, x=~`Leading logFC dim1`, y=~`Leading logFC dim2`,
                   color=formula(paste0("~", factorToPlot)),
                   colors = pal,
                   type='scatter',
                   mode='markers+text',
                   text=~samples, textposition = "top right")
    }
    p <- p %>% 
      plotly::layout(title=main,
                     xaxis = list(title = 'Leading logFC dim1'),
                     yaxis = list(title = 'Leading logFC dim2'))
  }else{
    stop("We only support 3D or 2D mds plot.")
  }
  p
}


## REFAC, but function is currently unused. undefined object: sm
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
cumHistPlot = function(cts, png, colors, main="all transcripts"){

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
