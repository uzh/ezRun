###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Chip Stats
##' @template htmlFile-template
##' @seealso \code{\link{EzAppChipStats}}
ezMethodChipStats = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
  setwdNew(basename(output$getColumn("Report")))
  dataset = input$meta
  gff = ezLoadFeatures(param, types = 'gene')
  if (!is.null(gff) && nrow(gff) == 0){
    writeErrorHtml(htmlFile, param = param, dataset = dataset,
                   error = list(error = paste("No features found in given feature file:<br>", 
                                              param$ezRef["refFeatureFile"])))
    return("Error")
  }
  
  ezMclapply(dataset$"BAM [File]", 
             Create_ChIP_QCPlots_ind, param=param, gff=gff, maxX=20, range=c(1,100), 
             mc.cores=param$cores, mc.preschedule =FALSE, mc.set.seed=FALSE)
  
  ###SSD-Coverage:
  files = list.files('.', pattern = 'ssdCoverage.txt')
  data = vector(mode = 'numeric', length = length(files))
  names(data) = sub("_ssdCoverage.txt", "", files)
  for (i in 1:length(files)){
    data[i] = read.table(files[i])$V1
  }
  png('ssdCoveragePlot.png',width=800,height=800)
  par(mar = c(5, 12, 4, 2) + 0.1)
  barplot(data, las=2, horiz=T, col='lightblue', border='lightblue', space=0.1, cex.names=0.9, xlab = 'SSD', main = 'SSD of Coverage')
  par(mar = c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  ###HTML-Report
  html = openHtmlReport(htmlFile, param=param, title=paste("ChIP Statistics:", param$name),
                        dataset=dataset)
  if(param[['paired']]){
    png_images = list.files('.', pattern='fragmentSize')  
    ezWrite("<h1>Density-Plot of Fragment Size</h1>", con=html)
    writeImageRowToHtml(png_images, con=html)
  }
  png_images = list.files('.', pattern='SingleEnrichmentPlot')
  ezWrite("<h1>General Enrichment-Plot</h1>", con=html)
  writeImageRowToHtml(png_images, con=html)
  png_images = list.files('.', pattern='StartPositionPlot')
  ezWrite("<h1>Frequency of identical Reads per Start Position</h1>", con=html)
  writeImageRowToHtml(png_images, con=html)
  png_images = list.files('.', pattern='TSSPlot')
  ezWrite("<h1>TSS-Heatmap</h1>", con=html)
  writeImageRowToHtml(png_images, con=html)
  png_images = list.files('.', pattern='cpgDensityPlot')
  if(length(png_images)>0){
    ezWrite("<h1>CpG-Density-Plot</h1>", con=html)
    writeImageRowToHtml(png_images, con=html)
  }
  
  enrichmentFiles = list.files(path = '.', pattern='EnrichedCoverage.txt')
  enrichmentData = list()
  for(i in 1:length(enrichmentFiles)) {
    enrichmentData[[i]] = read.table(enrichmentFiles[i], header = T, stringsAsFactors = F,quote = '', sep ='\t')
  }
  names(enrichmentData) = sub("_EnrichedCoverage.txt", "", enrichmentFiles)
  cols = rainbow(n = length(enrichmentFiles))
  png('ComprehensiveEnrichmentPlot.png',width=800,height=800)
  plot(enrichmentData[[1]][,1],log2(enrichmentData[[1]][,2]),type='l',xlim=c(0,20),col = cols[1], 
       lwd = 2, main ='Multi Sample Enrichment-Plot',xlab='Normalized Enrichment Level of Reads',ylab='Frequency (log2)')
  for(i in 2:length(enrichmentData)){
    lines(enrichmentData[[i]][,1], log2(enrichmentData[[i]][,2]), col = cols[i], lwd = 2)
  }
  legend('topright', lwd = 2, names(enrichmentData), pch = '', lty = 1, cex = 0.9, col = cols)
  dev.off()
  
  ####CLUSTERING Bin-Data:
  binFiles = list.files(path = '.', pattern = 'BinData.txt')
  tmp = read.table(binFiles[1], header = F, quote = '', sep = '\t')
  clusterData = matrix(0, nrow(tmp), length(binFiles))
  colnames(clusterData) = sub("_BinData.txt", "", binFiles)
  rownames(clusterData) = tmp$V1
  clusterData[,1] = tmp$V2
  for (i in 2:length(binFiles)) {
    clusterData[,i] = read.table(binFiles[i], header = F, quote = '', sep = '\t')$V1
  }
  png('ClusteringAllBins.png',width=800,height=800)
  clusterData(clusterData,title='Clustering of 500 Bp Counts')
  dev.off()
  
  remove(tmp)
  
  ezWrite("<h1>Multi Sample Plots</h1>", con=html)
  png_images = list.files('.', pattern='ComprehensiveEnrichmentPlot')
  png_images = c(png_images,list.files('.', pattern='ssdCoveragePlot'))
  png_images = c(png_images,list.files('.', pattern='ClusteringAllBins'))
  writeImageRowToHtml(png_images, con=html)
  
  closeHTML(html)
  print(warnings())
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodChipStats()
##' @seealso \code{\link{ezMethodChipStats}}
EzAppChipStats <-
  setRefClass("EzAppChipStats",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodChipStats
                  name <<- "EzAppChipStats"
                  appDefaults <<- rbind(estimatedFragmentSize=ezFrame(Type="numeric",  DefaultValue='200',  Description="estimated size of the DNA fragments"),
                                        flank=ezFrame(Type="integer",  DefaultValue="200",  Description="flanking around the transcription start sites used for plotting"))
                }
              )
  )

galp2gal = function(galp){ 
  system('echo Function galp2gal, conversion of paired-end data into single end data... \n')
  betweenPairCigar = paste0(abs(start(right(galp)) - end(left(galp)) + 1), "N")
  galcigar = paste0(cigar(left(galp)), betweenPairCigar, cigar(right(galp)))
  gal = GAlignments(
    seqnames = seqnames(galp),
    pos = start(left(galp)),
    cigar = galcigar,
    strand = strand(left(galp)),
    names = names(left(galp)),
    seqlengths = seqlengths(galp))
  # in case where end of a read exceeds chrom length
  # (which occurs for improper paired reads with large gap)  
  idx = which(end(gal) < seqlengths(gal)[as.character(seqnames(gal))])
  if(length(idx) < length(gal))    
    warning(sprintf("%d read-pairs with end larger than chromosome length are discarded", length(gal) - length(idx) + 1))
  system('echo Function galp2gal done \n')
  return(gal[idx,])    
}

readBam = function(file,isPaired=F){
  requireNamespace("limma")
  requireNamespace("rtracklayer")
  system('echo Function readBam \n')
  if(isPaired){
    system('echo option isPaired=T \n')
    bam = readBamGappedAlignmentPairs(file, use.names=T)
    system('echo bam read, conversion in progress...\n')
    bam = galp2gal(bam)
  } else {
    system('echo option isPaired=F \n')  
    bam = readGAlignments(file)
    #bam = readGAlignmentsFromBam(file, use.names=T)
  }
  sample.name = rev(strsplit2(file,'/'))[1]
  sample.name = gsub('.bam$', '', sample.name)
  sample.list = list(granges(bam))
  names(sample.list) = sample.name
  sample.list = GRangesList(sample.list)
  system('echo Function readBam: done \n')
  return(sample.list)
}

fragmentSize = function(myBam, isPaired = F){
  system('echo Function fragmentSize \n')  
  if(isPaired){
    fragmentSize=width(myBam[[1]])
    medianFS = median(fragmentSize)
    minQuantile=quantile(fragmentSize, 0.1)
    maxQuantile=quantile(fragmentSize, 0.9)
    sample_fragmentSize = fragmentSize[which(fragmentSize>minQuantile & fragmentSize < maxQuantile)]
    system('echo sample_fragmentSize calculated, plotting... \n')
    png(paste0(names(myBam),"_fragmentSize.png"))
    d = density(sample_fragmentSize)
    plot(d, main=names(myBam),xlab=paste0("Median of FragmentSize:",medianFS))
    polygon(d, col="red", border="black") 
    legend("topright",c(names(myBam)))
    dev.off()
  }  else  {
    warning("cannot compute fragmentSize for single-end bam files")      
  }
  system('echo Function fragmentSize: done \n')
}

MakeEnrichmentPlot = function(myBam, isPaired=F, estimatedFragmentSize = 200){
  system('echo Function MakeEnrichmentPlot \n')
  if(isPaired){
    seq.len = median(width(myBam[[1]]))
    system (paste('echo Estimated Median FragmentSize:', seq.len, '\n'))
  } else {
    seq.len = estimatedFragmentSize
    system (paste('echo SE data, defined FragmentSize:', seq.len, '\n'))
  }
  png(paste0(names(myBam), "_SingleEnrichmentPlot.png"))
  data = Repitools::enrichmentPlot(myBam, seq.len, cols = c("brown"), lwd = 4, main = names(myBam))
  dev.off()
  write.table(data[[1]], paste(names(myBam),'_EnrichedCoverage.txt', sep = ''),quote = F, row.names = F, sep = '\t')
  system('echo Function MakeEnrichmentPlot done \n')
}

MakeCpGDensityPlot = function(myBam, isPaired = F, build, estimatedFragmentSize = 200){
  system('echo Function MakeCpGDensity Plot \n')
  bsgenome_table = getBSgenomes()
  if(build %in% bsgenome_table$build){
    bsgenome_table = unlist(bsgenome_table[which(bsgenome_table$build == build),])
    system(paste('echo Found:', bsgenome_table, '\n'))
    ##require(bsgenome_table[1], character.only=T)
    if(isPaired){
      seq.len = median(width(myBam[[1]]))
      system(paste('echo Estimated Median FragmentSize:', seq.len, '\n'))
    } else {
      seq.len = estimatedFragmentSize
      system(paste('echo SE data, defined FragmentSize:', seq.len, '\n'))
    }
    png(paste0(names(myBam), "_cpgDensityPlot.png"))
    Repitools::cpgDensityPlot(myBam, organism = get(bsgenome_table[2]), w.function = "none", seq.len,
                              cols = c("black"), xlim = c(0, 30), lwd = 2, main = names(myBam))
    dev.off()
    system('echo Function MakeCpGDensityPlot done \n')
  } else { 
    system('echo Unsupported organism - cannot generate CpG-Density-Plot \n')
  }
}  

getBSgenomes = function(){
  system('echo Function getBSgenomes \n')
  listall_availablePackages = (.packages(all.available=TRUE))
  bs_packages = listall_availablePackages[grep("^BSgenome.", listall_availablePackages)]
  organism = strsplit2(bs_packages,"\\.")[,2]
  build = strsplit2(bs_packages,"\\.")[,4]
  bsgenome_table = data.frame(bs_packages, organism, build, stringsAsFactors=F)
  system('echo Function getBSgenomes done \n')
  return(bsgenome_table)
}

CoverageVarFunction = function(myBam){
  system('echo Function CoverageVarFunction \n')
  myCov = htSeqTools::ssdCoverage(myBam)
  write.table(myCov, file = paste(names(myBam), "_ssdCoverage.txt",sep= ""), quote = F, row.names = F, col.names = F, sep = "")
  system('echo Function CoverageVarFunction done \n')
  return(myCov)
}

StartPosTableFunction = function(myBam, maxX=20){
  system('echo Function StartPosTableFunction \n')
  startPosTable = table(table(paste(seqnames(myBam[[1]]), start(myBam[[1]]), strand(myBam[[1]]), sep='_')))
  png(paste0(names(myBam), "_StartPositionPlot.png"))
  plot(x = as.numeric(names(startPosTable)),y = log2(startPosTable),type='l', ylab=c('Frequency (log2)'), xlab=c('Reads_per_Start_Position'), xlim=c(1,maxX), main=names(myBam), lwd=3)
  legend("topright", c(names(myBam)))
  dev.off()  
  system('echo Function StartPosTableFunction done \n')
}

MyTable = function(x, range) {
  dat = vector('numeric', 1+(range[2]-range[1]))
  names(dat) = seq(range[1], range[2],1)
  fill = table(x)
  dat[names(fill)] = fill
  return(dat)
}

remove_outliers = function(x, na.rm = TRUE, ...) {
  qnt = quantile(x, probs=c(.25,.9), na.rm = na.rm, ...)
  H = 1.5 * IQR(x, na.rm = na.rm)
  y = x
  y[x < (qnt[1] - H)] = NA
  y[x > (qnt[2] + H)] = NA
  y
}

createTSSPlot = function(myBam, gff, flank, name, range=c(1,100)){
  system('echo Function createTSSPlot \n')
  requireNamespace("rtracklayer")
  cov=coverage(myBam)
  system(paste0('echo ',names(myBam)[1]))
  binMyBam(cov = cov, binLength = 500, sampleName = names(myBam))
  GffSummary = data.frame(
    Start = gff$start, 
    End = gff$end,
    Chr =gff$seqid,
    Strand = gff$strand, stringsAsFactors = F)
  
  positionData = matrix(0, nrow(GffSummary), 2*flank)
  removeList = c()
  for (i in 1:nrow(GffSummary)){
    ezHr = as.character(GffSummary$Chr[i])
    if(GffSummary$Strand[i]=='+'){
      mypos = c(GffSummary$Start[i]-flank, GffSummary$Start[i]+flank)
    }
    else {
      mypos = c(GffSummary$End[i]+flank, GffSummary$End[i]-flank)
    }
    mypos[which(mypos<0)] = 1
    mypos[which(mypos>length(cov[[ezHr]]))] = length(cov[[ezHr]])
    if(abs(diff(mypos))==2*flank){
      positionData[i,] = as.numeric(cov[[ezHr]][mypos[1]:mypos[2]])[-c(1)]
    }  
    else {
      removeList = c(removeList,i)
    }
  }
  if(!is.null(removeList)){
    positionData = positionData[-c(removeList),]
  }
  positionData = positionData[which(rowSums(positionData)>3*flank),]
  positionData = apply(positionData, 2, remove_outliers)
  oldRange = range[2]
  range[2] = max((max(positionData, na.rm=T) + 5), oldRange*0.5)
  posDataCounts = matrix(0, range[2]-range[1]+1, ncol(positionData))
  for(i in 1:ncol(positionData)){
    posDataCounts[,i] = MyTable(shrinkToRange(positionData[!is.na(positionData[,i]),i], range), range)    
  }
  colnames(posDataCounts) = seq(-1*flank, flank-1, 1)
  #require(gplots)
  MyCols = colorRampPalette(c('black','white'))(oldRange)
  file = paste0(names(myBam), "_TSSPlot.png")
  png(file, width=640, height=640)
  image(t(sqrt(posDataCounts[-c(1:2),]+1)), col=MyCols, axes=F)
  title(main=names(myBam), col.main="black",
        sub="", col.sub="blue",
        xlab="Relative Genomic Position around TSS", ylab="Normalized Coverage",
        col.lab="black", cex.lab=1) 
  legend("topright",c(names(myBam)))
  axis(2, at=c(0,1)) 
  axis(1, at=seq(0,1,0.25), labels=c(paste0('-', flank), paste0('-',flank/2), 0, flank/2,flank)) 
  dev.off()
  system('echo Function createTSSPlot done \n')
} 

binMyBam = function(cov, binLength, sampleName){
  MyBins = c()
  if (binLength>1){
    for (i in 1:length(cov)){
      y = rep(runValue(cov[[i]]), runLength(cov[[i]]))
      x = c(1:length(y))
      bx = seq(0, rev(x)[1], binLength)
      result = matrixStats::binMeans(y, x=x, bx=bx)
      cat(i, '\n')
      names(result) = paste0(names(cov)[i], '_', bx[-length(bx)])
      MyBins = c(MyBins, result)
    }
  } else {
    for (i in 1:length(cov)){
      MyBins = c(MyBins,rep(runValue(cov[[i]]),runLength(cov[[i]])))
    }
  }
  write.table(MyBins,paste0(sampleName,'_',binLength,'_BinData.txt'),col.names=F,quote=F,sep='\t')
}

clusterData = function(data, distmethod='pearson', clustermethod='ward', title = title){
  n = round(nrow(data)*0.2)
  keep = names(sort(apply(data,1,var),decreasing=T,na.last=NA))[1:n]
  data = data[keep,]
  c = cor(data, method = distmethod, use='complete.obs')
  d = as.dist(1-c)
  hr = hclust(d, method = clustermethod, members=NULL) 
  par(mfrow = c(1, 1))
  plot(hr, hang = -1,main=title)
}


Create_ChIP_QCPlots_ind = function(file, param, maxX=20, gff=gff, name='', range=c(1,100)){
  isPaired = as.logical(param[['paired']])
  file = paste(param[['dataRoot']],file,sep='/')
  estimatedFragmentSize = param[['estimatedFragmentSize']]
  build = param[['refBuild']]
  flank = param[['flank']]
  myBam = readBam(file, isPaired)
  #createBigWig(aligns=myBam, name=basename(file))  ## TODO: refactor
  fragmentSize(myBam, isPaired)
  MakeEnrichmentPlot(myBam,isPaired, estimatedFragmentSize)
  #MakeCpGDensityPlot(myBam,isPaired, build, estimatedFragmentSize)
  StartPosTableFunction(myBam, maxX)
  CoverageVarFunction(myBam)
  createTSSPlot(myBam, gff, flank, name='', range)
}
