###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFastqScreen = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  dataset = input$meta
  # fastqscreen part
  files = input$getFullPaths("Read1")
  resultFiles = executeFastqscreenCMD(param, files)
  fastqData = collectFastqscreenOutput(dataset, files, resultFiles)
  # bowtie2 reference part
  if ("Adapter1" %in% input$colNames && all(!is.na(input$getColumn("Adapter1")))){
    input = ezMethodTrim(input = input, param = param)
  } else {
    param$trimAdapter = FALSE
    input = ezMethodTrim(input = input, param = param)
  }
  countFiles = executeBowtie2CMD(param, input)
  speciesPercentageTop = collectBowtie2Output(param, input$meta, countFiles)
  setwdNew(basename(output$getColumn("Report")))
  fastqscreenReport(dataset, param, htmlFile, fastqData, speciesPercentageTop)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastqScreen(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
##' @section Functions:
##' \itemize{
##'   \item{\code{executeFastqscreenCMD(param, files): }}
##'   {Executes the fastqscreen command with given parameters on the input files.}
##'   \item{\code{collectFastqscreenOutput(dataset, files, resultFiles): }}
##'   {Collects the fastqscreen output after the result files have been obtained by \code{executeFastqscreenCMD()}.}
##'   \item{\code{executeBowtie2CMD(param, input): }}
##'   {Executes the bowtie2 command with given parameters on the input files.}
##'   \item{\code{collectBowtie2Output(param, dataset, countFiles): }}
##'   {Collects the bowtie2 output after the count files have been obtained by \code{executeBowtie2CMD()}.}
##'   \item{\code{fastqscreenReport(dataset, param, htmlFile="00index.html", fastqData, speciesPercentageTop): }}
##'   {Generates a report with plots and other information about the outcome of the run.}
##' }
EzAppFastqScreen <-
  setRefClass("EzAppFastqScreen",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFastqScreen
                  name <<- "EzAppFastqScreen"
                  appDefaults <<- rbind(nTopSpecies=ezFrame(Type="integer",  DefaultValue=10,  Description="number of species to show in the plots"),
                                        confFile=ezFrame(Type="character",  DefaultValue="",  Description="the configuration file for fastq screen"),
                                        minAlignmentScore=ezFrame(Type="integer",  DefaultValue="-20",  Description="the min alignment score for bowtie2"),
                                        trimAdapter=ezFrame(Type="logical",  DefaultValue=TRUE,  Description="whether to search for the adapters and trim them"))
                }
              )
  )

executeFastqscreenCMD = function(param, files){
  confFile = paste0(FASTQSCREEN_CONF_DIR, param$confFile)
  opt = ""
  if (param$nReads > 0){
    opt = paste(opt, "--subset", ezIntString(param$nReads))
  }
  cmd = paste(FASTQSCREEN, opt, " --threads", param$cores, " --conf ", confFile,
              paste(files, collapse=" "), "--outdir . --aligner bowtie2",
              "> fastqscreen.out", "2> fastqscreen.err")
  ezSystem(cmd)
  resultFiles = paste0(sub(".fastq$", "", sub(".gz$", "", basename(files))), "_screen.txt") ## remove the suffix .fastq[.gz] with _screen.txt
  names(resultFiles) = names(files)
  return(resultFiles)
}

collectFastqscreenOutput = function(dataset, files, resultFiles){
  fastqData = list()
  fastqData$MappingRate = numeric()
  fastqData$Reads = integer()
  fastqData$CommonResults = list()
  for(nm in rownames(dataset)){
    cat('Process ',files[nm], " - ", resultFiles[nm], ' :')
    x = ezRead.table(resultFiles[nm], skip=1, stringsAsFactors=F, blank.lines.skip=T, fill=T, row.names=NULL)
    fastqData$Reads[nm] = x$"#Reads_processed"[1]
    UnmappedReads = as.numeric(unlist(strsplit(x$Genome[nrow(x)], split = " "))[2])
    fastqData$MappingRate[nm] = round((100 - UnmappedReads), digits = 2)
    rownames(x) = x$Genome
    x = x[-nrow(x), grep('%.*hit',colnames(x))]
    fastqData$CommonResults[[nm]] = x
  }
  return(fastqData)
}

executeBowtie2CMD = function(param, input){
  r1Files = input$getFullPaths("Read1")
  if (param$paired){
    r2Files = input$getFullPaths("Read2")
  }
  countFiles = character()
  for (nm in names(r1Files)){
    countFiles[nm] = paste0(nm, "-counts.txt")
    bowtie2options = param$cmdOptions
    writeLines("ReadID\tRefSeqID\tAlignmentScore", countFiles[nm])
    if(!param$paired){
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                  " -U ", r1Files[nm], bowtie2options ,"-p",param$cores,
                  "--no-unal --no-hd", "2> ", paste0(nm, "_bowtie2.err"),
                  "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm])
    } else {
      cmd = paste(file.path(BOWTIE2_DIR,'bowtie2'),"-x",REFSEQ_mRNA_REF, 
                  " -1 ", r1Files[nm]," -2 ", r2Files[nm], bowtie2options, "-p",param$cores,
                  "--no-discordant --no-mixed --no-unal --no-hd",
                  "2> ", paste0(nm, "_bowtie2.err"),
                  "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm])
    }
    ezSystem(cmd)
  }
  return(countFiles)
}

collectBowtie2Output = function(param, dataset, countFiles){
  tax2name = read.table('/srv/GT/reference/RefSeq/mRNA/20150301/Annotation/tax2name.txt',header=F,stringsAsFactors=F,sep='|', 
                        colClasses="character",quote='', comment.char="")
  colnames(tax2name) = c('TAX_ID','Name')
  tax2name$TAX_ID = sub("\t", "", tax2name$TAX_ID)
  tax2name$Name = gsub('\t','',tax2name$Name)
  rownames(tax2name) = tax2name$TAX_ID
  speciesPercentageTop = list()
  for (nm in names(countFiles)){
    countData = ezRead.table(countFiles[nm], row.names = NULL)
    bestScores = tapply(countData$AlignmentScore, countData$ReadID, max)
    countData = countData[countData$AlignmentScore == bestScores[countData$ReadID], , drop=FALSE]
    countData = countData[countData$AlignmentScore >= param$minAlignmentScore, , drop=FALSE]
    if (nrow(countData) > 0){
      countData$species = sub("_.*", "", countData$RefSeqID)
      speciesHitsPerRead = tapply(countData$species, countData$ReadID, unique)
      uniqSpeciesHitsPerRead = names(speciesHitsPerRead)[sapply(speciesHitsPerRead, length) == 1]
      ###Result UniqHits:
      uniqSpeciesHits = sort(table(unlist(speciesHitsPerRead[uniqSpeciesHitsPerRead])), decreasing = T)
      ###Results MultipleHits:
      multipleSpeciesHitsPerRead = countData[!(countData$ReadID %in% uniqSpeciesHitsPerRead), ]
      by = paste(multipleSpeciesHitsPerRead$ReadID, multipleSpeciesHitsPerRead$species,sep='_')
      ##multipleSpeciesHits = sort(table(multipleSpeciesHitsPerRead$species[!duplicated(by)]), decreasing=T) # is equivalent to the row below
      multipleSpeciesHits = sort(table(tapply(multipleSpeciesHitsPerRead$species,by,unique)),decreasing = T)
      
      if (length(uniqSpeciesHits) > param$nTopSpecies){
        topSpeciesUniq = uniqSpeciesHits[1:param$nTopSpecies]
      } else {
        topSpeciesUniq = uniqSpeciesHits
      }
      ## Special case where all hits are multi hits --- in that case we sort according to the multi-hits
      if (length(uniqSpeciesHits) == 0){
        topSpeciesUniq = rep(0, min(param$nTopSpecies, length(multipleSpeciesHits)))
        names(topSpeciesUniq) = names(multipleSpeciesHits)[1:min(param$nTopSpecies, length(multipleSpeciesHits))]
      }

      multipleSpeciesHits[setdiff(names(topSpeciesUniq), names(multipleSpeciesHits))] = 0
      topSpeciesMultiple = multipleSpeciesHits[names(topSpeciesUniq)]

      if(param$nReads > 0){
        totalCount = param$nReads
      } else {
        totalCount = dataset[nm, 'Read Count']
      }
      taxIds = names(topSpeciesUniq)
      taxNames = tax2name[taxIds, 'Name']
      hasNoName = is.na(taxNames)
      taxNames[hasNoName] = taxIds[hasNoName]
      x = ezFrame(UniqueSpeciesHits=signif(100 * as.matrix(topSpeciesUniq)/totalCount, digits=4),
                  MultipleSpeciesHits=signif(100 * as.matrix(topSpeciesMultiple)/totalCount, digits=4),
                  row.names = taxNames)
      speciesPercentageTop[[nm]] = x
    } else {
      speciesPercentageTop[[nm]] = NULL
    }
  }
  return(speciesPercentageTop)
}

fastqscreenReport = function(dataset, param, htmlFile="00index.html", fastqData, speciesPercentageTop){
  titles = list()
  titles[["FastQ Screen"]] = paste("FastQ Screen:", param$name)
  doc = openBsdocReport(title=titles[[length(titles)]])
  
  addDataset(doc, dataset, param)
  
  settings = character()
  settings["Configuration File:"] = param$confFile
  settings["RefSeq mRNA Reference:"] = REFSEQ_mRNA_REF
  settings["FastqScreen Version:"] = basename(dirname(FASTQSCREEN))
  settings["Bowtie2 Version:"] = basename(BOWTIE2_DIR)
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
    bplt = barplot(fastqData$MappingRate, las=2, ylim=c(0,100), ylab="MappedReads in %", main="MappingRate", col="royalblue3",
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
    bplt = barplot(fastqData$Reads/1000, las=2, ylab="#Reads in K", main="ProcessedReads", col="lightblue",
            names.arg=rep('',length(ezSplitLongLabels(names(fastqData$MappingRate)))))
    text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = ezSplitLongLabels(names(fastqData$MappingRate)), xpd = TRUE)
  })
  readsLink = ezImageFileLink(plotCmd, file="Reads.png", width=min(600 + (nrow(dataset)-10)* 30, 2000)) # nSamples dependent width
  addFlexTable(doc, ezGrid(cbind(mappingRateLink, readsLink)))
  
  screenLinks = list()
  detectedSpeciesLinks = list()
  for (nm in rownames(dataset)){
    plotCmd = expression({
      par(mar=c(10.1, 4.1, 4.1, 2.1))
      bplt = barplot(t(fastqData$CommonResults[[nm]]), las=2, ylim=c(0,100), 
                     legend.text=T, ylab="Mapped Reads in %", main=nm, names.arg=rep('', nrow(fastqData$CommonResults[[nm]])))
      text(x = bplt, y = par("usr")[3] - 2, srt = 45, adj = 1, labels = rownames(fastqData$CommonResults[[nm]]), xpd = TRUE)
      
      
    })
    screenLinks[[nm]] = ezImageFileLink(plotCmd, name = nm, plotType = "-rRNA-countsBySpecies-barplot", width=400, height=400)
    
    plotCmd = expression({
      par(mar=c(10.1, 4.1, 4.1, 2.1))
      x = speciesPercentageTop[[nm]]
      if (is.null(x)) x = 0
      bplot = barplot(t(x), col=c("royalblue3", "lightblue"), las=2, ylim=c(0,100),
                      legend.text=T, ylab="Mapped Reads in %", main=nm, names.arg=rep('',nrow(x)) )
      text(y=t(x)[ 1,] + 2, x=bplot, font = 2, labels=t(x)[ 1, ], cex=1.1,col='black', xpd=TRUE)
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
  closeBsdocReport(doc=doc, file=htmlFile, titles)
}
