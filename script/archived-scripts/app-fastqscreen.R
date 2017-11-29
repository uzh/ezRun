fastqscreenReport = function(dataset, param, htmlFile="00index.html", fastqData, fastqDataAdapters, speciesPercentageTop, speciesPercentageTopVirus=NULL){
  titles = list()
  titles[["FastQ Screen"]] = paste("FastQ Screen:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  addDataset(doc, dataset, param)
  
  settings = character()
  settings["Configuration File:"] = param$confFile
  settings["RefSeq mRNA Reference:"] = REFSEQ_mRNA_REF
  settings["FastqScreen Version:"] = basename(dirname("$FastQScreen"))
  settings["Bowtie2 Version:"] = basename(Sys.getenv("Bowtie2"))
  settings["Bowtie2 Parameters:"] = param$cmdOptions
  settings["Minimum AlignmentScore:"] = param$minAlignmentScore
  settings["TopSpecies:"] = param$nTopSpecies
  addFlexTable(doc, ezFlexTable(as.data.frame(settings), add.rownames=TRUE))
  
  titles[["rRNA-Check"]] = "FastqScreen Mapping Rates"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  titles[["Per Dataset"]] = "Overview"
  addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    bplt = barplot(fastqData$MappingRate, las=2, ylim=c(0,100), ylab="MappedReads in %", main="Overall MappingRate", col="royalblue3",
                   names.arg=rep('',length(ezSplitLongLabels(names(fastqData$MappingRate)))))
    if(min(fastqData$MappingRate) < 8){
      text(y=fastqData$MappingRate+2, font=2, x=bplt, labels=as.character(fastqData$MappingRate), cex= 1, xpd=TRUE)
    } else {
      text(y=fastqData$MappingRate-5, font=2, x=bplt, labels=as.character(fastqData$MappingRate), cex= 1.1, col='white', xpd=TRUE)
    }
    text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = ezSplitLongLabels(names(fastqData$MappingRate)), xpd = TRUE)
  })
  
  mappingRateLink = ezImageFileLink(plotCmd, file="MappingRate.png", width=min(600 + (nrow(dataset)-10)* 30, 2000)) # nSamples dependent width
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    bplt = barplot(fastqDataAdapters$MappingRate, las=2, ylim=c(0,100), ylab="MappedReads in %", main="MappingRate to Adapters without trimming", col="royalblue3",
                   names.arg=rep('',length(ezSplitLongLabels(names(fastqDataAdapters$MappingRate)))))
    if(min(fastqDataAdapters$MappingRate) < 8){
      text(y=fastqDataAdapters$MappingRate+2, font=2, x=bplt, labels=as.character(fastqDataAdapters$MappingRate), cex= 1, xpd=TRUE)
    } else {
      text(y=fastqDataAdapters$MappingRate-5, font=2, x=bplt, labels=as.character(fastqDataAdapters$MappingRate), cex= 1.1, col='white', xpd=TRUE)
    }
    text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = ezSplitLongLabels(names(fastqDataAdapters$MappingRate)), xpd = TRUE)
  })
  
  mappingRateAdaptersLink = ezImageFileLink(plotCmd, file="MappingRateAdapters.png", width=min(600 + (nrow(dataset)-10)* 30, 2000)) # nSamples dependent width
  
  plotCmd = expression({
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    bplt = barplot(fastqData$Reads/1000, las=2, ylab="#Reads in K", main="ProcessedReads", col="lightblue",
                   names.arg=rep('',length(ezSplitLongLabels(names(fastqData$MappingRate)))))
    text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = ezSplitLongLabels(names(fastqData$MappingRate)), xpd = TRUE)
  })
  readsLink = ezImageFileLink(plotCmd, file="Reads.png", width=min(600 + (nrow(dataset)-10)* 30, 2000)) # nSamples dependent width
  addFlexTable(doc, ezGrid(cbind(mappingRateLink, mappingRateAdaptersLink, readsLink)))
  
  screenLinks = list()
  detectedSpeciesLinks = list()
  for (nm in rownames(dataset)){
    plotCmd = expression({
      par(mar=c(10.1, 4.1, 4.1, 2.1))
      x = fastqData$CommonResults[[nm]]
      if (nrow(x) > 0){
        bplt = barplot(t(x), las=2, ylim=c(0,100), 
                       legend.text=T, ylab="Mapped Reads in %", main=nm, names.arg=rep('', nrow(x)))
        text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = rownames(x), xpd = TRUE)
      } else {
        plot(1,1, type="n", axes=FALSE, ann=FALSE)
        text(1,1, "no hits found")
      }
    })
    screenLinks[[nm]] = ezImageFileLink(plotCmd, name = nm, plotType = "-rRNA-countsBySpecies-barplot", width=400, height=400)
    
    plotCmd = expression({
      par(mar=c(10.1, 4.1, 4.1, 2.1))
      x = speciesPercentageTop[[nm]]
      if (is.null(x)) x = matrix(0, 2, 1, dimnames=list(c('UniqueSpeciesHits','MultipleSpeciesHits'),'Misc'))
      bplot = barplot(t(x), col=c("royalblue3", "lightblue"), las=2, ylim=c(0,100),
                      legend.text=T, ylab="Mapped Reads in %", main=nm, names.arg=rep('',nrow(x)) )
      text(y=t(x)[ 1,] + 5, x=bplot, font = 2, labels=t(x)[ 1, ], cex=1.1, col='black')
      text(x = bplot, y = par("usr")[3] - 2, srt = 45, adj = 1, 
           labels = rownames(x), xpd = TRUE)
    })
    detectedSpeciesLinks[[nm]] = ezImageFileLink(plotCmd, name=nm, plotType="-mRNA-countsBySpecies-barplot", width=400, height=400)
  }
  IMAGESperROW = 4
  if (ezIsSpecified(screenLinks)){
    titles[["Per Sample"]] = "Per Sample"
    addTitle(doc, titles[[length(titles)]], 3, id=titles[[length(titles)]])
    if(length(screenLinks) <= IMAGESperROW){
      addFlexTable(doc, ezGrid(rbind(screenLinks)))
    } else {
      addFlexTable(doc, ezGrid(rbind(screenLinks[1:IMAGESperROW])))
      for (i in 1:(ceiling(length(screenLinks)/IMAGESperROW)-1)){
        addFlexTable(doc, ezGrid(rbind(screenLinks[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(screenLinks))])))
      }
    }
  }
  if (ezIsSpecified(detectedSpeciesLinks)){
    titles[["Mapping to RefSeq mRNA"]] = "Mapping to RefSeq mRNA"
    addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
    if(length(detectedSpeciesLinks) <= IMAGESperROW){
      addFlexTable(doc, ezGrid(rbind(detectedSpeciesLinks)))
    } else {
      addFlexTable(doc, ezGrid(rbind(detectedSpeciesLinks[1:IMAGESperROW])))
      for (i in 1:(ceiling(length(detectedSpeciesLinks)/IMAGESperROW)-1)){
        addFlexTable(doc, ezGrid(rbind(detectedSpeciesLinks[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(detectedSpeciesLinks))])))
      }
    }
  }
  
  if(param[['virusCheck']]){
    detectedVirusLinks = list()
    for (nm in rownames(dataset)){
      plotCmd = expression({
        par(mar=c(18.1, 7.1, 2.1, 2.1))
        x = speciesPercentageTopVirus[[nm]]
        if (is.null(x)) x = matrix(0, 2, 1, dimnames=list(c('UniqueSpeciesHits','MultipleSpeciesHits'),'Misc'))
        bplot = barplot(t(x), col=c("royalblue3", "lightblue"), las = 2, ylim = c(0,100),
                        legend.text=T, ylab="Mapped Reads in %", main=nm, names.arg=rep('',nrow(x)) )
        text(y=t(x)[ 1,] + 5, x=bplot, font = 2, labels=t(x)[ 1, ], cex = 1.1, col = 'black')
        text(x = bplot, y = par("usr")[3] - 2, srt = 60, adj = 1, 
             labels = rownames(x), xpd = TRUE)
      })
      detectedVirusLinks[[nm]] = ezImageFileLink(plotCmd, name=nm, plotType="-refSeq-countsByVirus-barplot", width=400, height=400)
    }
    if (ezIsSpecified(detectedVirusLinks)){
      titles[["Mapping to RefSeq human pathogenic Viruses"]] = "Mapping to RefSeq human pathogenic Viruses"
      addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
      if(length(detectedVirusLinks) <= IMAGESperROW){
        addFlexTable(doc, ezGrid(rbind(detectedVirusLinks)))
      } else {
        addFlexTable(doc, ezGrid(rbind(detectedVirusLinks[1:IMAGESperROW])))
        for (i in 1:(ceiling(length(detectedVirusLinks)/IMAGESperROW)-1)){
          addFlexTable(doc, ezGrid(rbind(detectedVirusLinks[(i*IMAGESperROW+1):min((i+1)*IMAGESperROW,length(detectedVirusLinks))])))
        }
      }
    }
  }
  
  closeBsdocReport(doc=doc, file=htmlFile, titles)
}

