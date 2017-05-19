###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFastQCRmd = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  require(rmarkdown)
  setwdNew(basename(output$getColumn("Report")))
  dataset = input$meta
  samples = rownames(dataset)
  files = c()
  for (sm in samples){
    files[paste0(sm, "_R1")] = input$getFullPaths("Read1")[sm]
    if (!is.null(dataset$Read2)){
      files[paste0(sm, "_R2")] = input$getFullPaths("Read2")[sm]    
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

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "FastQC.Rmd", "FastQC_overview.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  
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
  for(i in 1:length(plots)){
    plotPage <- plotPages[i]
    pngs <- file.path(reportDirs, "Images", plots[i])
    render(input="FastQC_overview.Rmd", output_dir=".", output_file=plotPage)
  }
  
  ## establish the main report
  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset
  
  if (!is.null(dataset$"Read Count")){
    readCount = signif(dataset$"Read Count" / 1e6, digits=3)
    names(readCount) = rownames(dataset)
  } else {
    readCount = integer()
    for (i in 1:nFiles){
      x = ezRead.table(file.path(reportDirs[i], "fastqc_data.txt"), 
                       header=FALSE, nrows=7, fill=TRUE)
      readCount[names(files)[i]] = signif(as.integer(x["Total Sequences", 1]) / 1e6, digits=3)
    }
  }
  ans4Report[["Read Counts"]] <- readCount
  
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
  ans4Report[["Fastqc quality measures"]] <- tbl
  
  qualMatrixList = ezMclapply(files, getQualityMatrix, mc.cores=ezThreads())
  ans4Report[["Per Base Read Quality"]] <- qualMatrixList
  
  ## debug
  ## save(ans4Report, file="ans4Report.rda")
  
  ## generate the main reports
  render(input="FastQC.Rmd", output_dir=".", output_file=htmlFile)
  
  ezSystem(paste("rm -rf ", paste0(reportDirs, ".zip", collapse=" ")))
  return("Success")
}

plotReadCountToLibConcRmd = function(dataset, colname){
  if(colname %in% colnames(dataset) && nrow(dataset) > 1){
    if(!all(dataset[[colname]]==0)){
      dataset = dataset[order(dataset$'Read Count',decreasing = T),]
      dataset$'Read Count' = dataset$'Read Count'/10^6
      corResult = cor.test(dataset$'Read Count',dataset[[colname]],
                           method = 'spearman')
      regressionResult = lm(dataset[[colname]]~dataset$'Read Count')
      label = sub(' \\[.*','',colname)
        plot(dataset$'Read Count', dataset[[colname]], pch=c(18), 
             cex=1.5, main=label,
             xlab='ReadCount in Mio', ylab=sub('\\[.*','', colname), 
             xlim=c(0, max(dataset$'Read Count', na.rm=TRUE)*1.2), #min(dataset$'Read Count')*0.8
             ylim=c(min(dataset[[colname]], na.rm = TRUE) * 0.8, 
                    max(dataset[[colname]], na.rm=TRUE) * 1.2))
        legend("topright", paste('r=', round(corResult$estimate, 2)), bty="n")
        abline(regressionResult, col='red',lty=c(2))
        text(dataset$'Read Count', dataset[[colname]], pos=2,
             labels=rownames(dataset), cex=0.7, col='darkcyan')
    }
  }
  return(NULL)
}

plotQualityMatrixAsHeatmapGG2 = function(qualMatrixList, isR2=FALSE, 
                                         xScale=1, yScale=1){
  colorsGray = gray((30:256)/256)
  minPercent = 0
  maxPercent = sqrt(40)
  minDiff = -5
  maxDiff = 5
  plotList <- list() # to store all the ggplot objects
  
  ## test if R2 exists
  index=list("R1"=which(!isR2))
  
  if(any(isR2)){
    stopifnot(sum(isR2) == sum(!isR2))
    index[["R2"]]= which(isR2)
  }
  for(nm in names(index)){
    plotList[[nm]] <- list()
    idx = index[[nm]]
    ## Plot the color key for the average quality heatmap R1_1
    # colorKeyFile = paste0("averageReadsQuality-Key_", nm, ".png")
    by.label = 1
    at=seq(from=minPercent, to=maxPercent, by=by.label)
    p = ezColorLegendGG2(colorRange=c(minPercent, maxPercent), 
                    colors=colorsGray, vertical=FALSE, by.label=by.label,
                    at=at, labels=as.character(at^2))
    plotList[[nm]][["Avg Qual Colors"]] <- p
    
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
    p = plotQualityHeatmapGG2(result=sqrt(avgQual), 
                              colorRange=c(minPercent, maxPercent), 
                              colors=colorsGray, 
                              main=paste("averageReadsQuality", nm, sep="_"), 
                              xScale=xScale, yScale=yScale)
    plotList[[nm]][["Average"]] <- p
    
    ## plot the difference quality heatmap for R1_1
    at=seq(from=minDiff, to=maxDiff, by=by.label)
    p = ezColorLegendGG2(colorRange=c(minDiff, maxDiff), 
                         colors=getBlueRedScale(), vertical=FALSE, 
                         by.label=by.label,
                         at=at, labels=as.character(at))
    plotList[[nm]][["Diff Qual Colors"]] <- p
    
    for(sampleName in names(qualMatrixList[idx])){
      qm = qualMatrixList[[sampleName]]
      diffResult = signif(prop.table(qm,2)*100, digits=3) - avgQual[1:nrow(qm), 
                                                                    1:ncol(qm)]
      p = plotQualityHeatmapGG2(diffResult, colorRange=c(minDiff, maxDiff),
                                colors=getBlueRedScale(), 
                                main=paste("diffReadsQuality", sampleName, 
                                           sep="_"),
                                xScale=xScale, yScale=yScale)
      plotList[[nm]][[sampleName]] <- p
    }
  }
  return(plotList)
}

plotQualityHeatmapGG2 = function(result, name=NULL, colorRange=c(0,sqrt(40)), 
                                 colors=gray((1:256)/256), main=NULL, 
                                 xScale=1, yScale=1){
  require(ggplot2)
  require(reshape2)
  ## some ugly controls of labels
  labCol = seq(0, ncol(result), by=10)
  labCol[1] = 1
  labRow = seq(0, nrow(result)-1, by=5)
  
  result[result > colorRange[2]] = colorRange[2]
  result[result < colorRange[1]] = colorRange[1]
  toPlot <- melt(result)
  p <- ggplot(data=toPlot, aes(x=Var2, y=Var1)) + geom_raster(aes(fill=value)) +
    scale_fill_gradientn(colours = colors, limits=colorRange) + theme_bw() +
    scale_y_continuous(breaks=seq(1, nrow(result), by=5), labels=labRow,
                       expand = c(0, 0)) +
    scale_x_continuous(breaks=labCol, labels=labCol, expand = c(0, 0)) +
    theme(panel.border=element_blank(), panel.grid=element_blank(),
          legend.position="none", 
          axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Read Position") + ylab("Read Quality") + ggtitle(main)
  return(p)
}
