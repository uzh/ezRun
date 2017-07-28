###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFastQC = function(input=NA, output=NA, param=NA, 
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
    cmd = paste("fastqc", "--extract -o . -t", min(ezThreads(), 8), 
                param$cmdOptions,
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
      readCount[names(files)[i]] = signif(as.integer(x["Total Sequences", 1]) / 
                                            1e6, digits=3)
    }
  }
  ans4Report[["Read Counts"]] <- readCount
  
  for (i in 1:nFiles){
    smy = ezRead.table(file.path(reportDirs[i], "summary.txt"), row.names=NULL, header=FALSE)
    if (i == 1){
      rowNames = paste0("<a href=", reportDirs, "/fastqc_report.html>", 
                        names(files), "</a>")
      colNames = ifelse(smy[[2]] %in% names(plotPages),
                        paste0("<a href=", plotPages[smy[[2]]], ">", smy[[2]], 
                               "</a>"),
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
  #save(ans4Report, file="ans4Report.rda")
  
  ## generate the main reports
  render(input="FastQC.Rmd", output_dir=".", output_file=htmlFile)
  
  ezSystem(paste("rm -rf ", paste0(reportDirs, ".zip", collapse=" ")))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }
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

plotReadCountToLibConc = function(dataset, colname){
  if(colname %in% colnames(dataset) && nrow(dataset) > 1){
    if(!all(dataset[[colname]]==0)){
      dataset = dataset[order(dataset$'Read Count',decreasing = T),]
      dataset$'Read Count' = dataset$'Read Count'/10^6
      corResult = cor.test(dataset$'Read Count',dataset[[colname]],
                           method = 'spearman')
      regressionResult = lm(dataset[[colname]]~dataset$'Read Count')
      label = sub(' \\[.*','',colname)
        # plot(dataset$'Read Count', dataset[[colname]], pch=c(18), 
        #      cex=1.5, main=label,
        #      xlab='ReadCount in Mio', ylab=sub('\\[.*','', colname), 
        #      xlim=c(0, max(dataset$'Read Count', na.rm=TRUE)*1.2), #min(dataset$'Read Count')*0.8
        #      ylim=c(min(dataset[[colname]], na.rm = TRUE) * 0.8, 
        #             max(dataset[[colname]], na.rm=TRUE) * 1.2))
        #legend("topright", paste('r=', round(corResult$estimate, 2)), bty="n")
        #abline(regressionResult, col='red',lty=c(2))
        #text(dataset$'Read Count', dataset[[colname]], pos=2,
        #     labels=rownames(dataset), cex=0.7, col='darkcyan')
        
        ## plotly
      require(plotly)
        #a function to calculate your abline
        xmin <- min(dataset$'Read Count') - 5
        xmax <- max(dataset$'Read Count') + 5
        intercept <- regressionResult$coefficients[1]
        slope <- regressionResult$coefficients[2]
        p_abline <- function(x, a, b){
          y <- a * x + b
          return(y)
        }
        
        p <- plot_ly(x=dataset$'Read Count', y=dataset[[colname]], 
                 text=rownames(dataset)) %>%
           add_markers() %>%
           add_text(textposition = "top right") %>%
           layout(showlegend = FALSE)
        a <- list(
          x = max(dataset$'Read Count'),
          y = max(dataset[[colname]]),
          text = paste0("r=", round(corResult$estimate, 2)),
          xref = "x",
          yref = "y",
          showarrow = FALSE
        )
        p <- p %>% layout(
          shapes=list(type='line', line=list(dash="dash"),
                      x0=xmin, x1=xmax,
                      y0=p_abline(xmin, slope, intercept), 
                      y1=p_abline(xmax, slope, intercept)),
          annotations = a,
          title=label, yaxis = list(title = label),
          xaxis = list(title = "Counts [Mio]")
        )
        return(p)
        ## ggplot2
        # toPlot <- data.frame(x=dataset$`Read Count`,
        #                      y=dataset[[colname]],
        #                      label=rownames(dataset))
        # p <- ggplot(toPlot, aes(x=x, y=y, label=label)) +
        #   geom_point() + geom_text(hjust = 0, nudge_x = 0.05) +
        #   theme_bw()
    }
  }
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

getQualityMatrix = function(inputFile){
  ## files could be fastq, or gziped fastq.gz
  require("ShortRead", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)  
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

plateStatistics <- function(dataset,
                            colname=c("Read Count", 
                                      "LibConc_qPCR [Characteristic]",
                                      "LibConc_100_800bp [Characteristic]")){
  colsExist <- colname %in% colnames(dataset)
  if(any(!colsExist)){
    warning("No column ", colname[!colsExist], " in dataset!")
    colname <- colname[colsExist]
  }
  colsZero <- colSums(dataset[, colname, drop=FALSE]) == 0
  if(any(colsZero)){
    message("The column ", colname[colsZero], " is empty!")
    colname <- colname[!colsZero]
  }
  
  if(!is.null(dataset$`PlatePosition [Characteristic]`)){
    plateChar <- dataset$`PlatePosition [Characteristic]`
    ## Plate position should be in the format of *_C4;
    ## * is the plate number
    ## C is the column; 4 is the row
    require(stringr, quietly = TRUE)
    plateNumber <- sub("_[[:alpha:]]\\d+$", "", plateChar)
    datasetByPlate <- split(dataset, plateNumber)
    ans <- list()
    for(i in seq_len(length(datasetByPlate))){
      platePos <- str_extract(datasetByPlate[[i]]$`PlatePosition [Characteristic]`,
                              "_[[:alpha:]]\\d+$")
      if(any(is.na(platePos))){
        warning("The PlatePosition format is not supported!")
        return(NA)
      }
      plateRow <- str_extract(platePos, "[[:alpha:]]")
      plateCol <- str_extract(platePos, "\\d+")
      ans[[names(datasetByPlate)[i]]] <- list()
      for(oneCol in colname){
        counts <- datasetByPlate[[i]][[oneCol]]
        countMatrix <- ezMatrix(0, rows=LETTERS[1:which(LETTERS==max(plateRow))],
                                cols=seq_len(max(as.integer(plateCol))))
        for(j in seq_len(length(counts))){
          countMatrix[plateRow[j], plateCol[j]] <- counts[j]
        }
        ans[[names(datasetByPlate)[i]]][[oneCol]] <- countMatrix
      }
    }
    return(ans)
  }else{
    return(NA)
  }
}
