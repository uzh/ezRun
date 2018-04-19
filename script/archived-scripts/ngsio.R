loadCountDataset = function(input, param){
  require(tools)
  files = input$getFullPaths("Count")
  suffix = unique(toupper(file_ext(files)))
  if (length(suffix) > 1){
    return(list(error=paste("different file suffixes not supported: <br>",
                            paste(files, collapse="<br>"))))
  }
  x = ezRead.table(files[1])
  if (ezIsSpecified(param$expressionName)){
    columnName = param$expressionName
  } else {
    columnName = intersect(param$knownExpressionNames, colnames(x))[1]
  }
  if (!columnName %in% colnames(x)){
    return(list(error=paste0("Specified column name not found in data!<br>columnName: '", columnName, "'\n",
                             "<br>Available column names:<br>\n",
                             paste0("'", colnames(x), "'", collapse="<br>"),
                             "<br>Set the option columnName to one of the names above!")))
  }
  dataFeatureLevel = unique(input$getColumn("featureLevel"))
  stopifnot(length(dataFeatureLevel) == 1)
  if (ezIsSpecified(param$ezRef@refBuild)){
    seqAnno = ezFeatureAnnotation(param, rownames(x), dataFeatureLevel)
  } else {
    seqAnno = x[ , intersect(c("type", "gene_name", "gene_id", "transcript_id", "Description", "GO BP", "GO MF", "GO CC", "gc", "width"), colnames(x)), drop=FALSE]
  }
  signal = ezMatrix(0, rows=rownames(seqAnno), cols=names(files))
  columnNameStart = grep(paste(columnName, "[first"), colnames(x), fixed=TRUE, value=TRUE)
  if (length(columnNameStart) == 1){
    signalStart = signal
  } else {
    signalStart = NULL
  }
  columnNameEnd = grep(paste(columnName, "[last"), colnames(x), fixed=TRUE, value=TRUE)
  if (length(columnNameEnd) == 1){
    signalEnd = signal
  } else {
    signalEnd = NULL
  }
  for (i in 1:length(files)){
    message("loading file: ", files[i])
    x = ezRead.table(files[i], strip.white = FALSE)
    if(!setequal(rownames(x), rownames(seqAnno))){
      if (all( rownames(seqAnno) %in% rownames(x))){
        warning("inconsistent ID set")
      } else {
        stop("later arrays have IDs not present in the first array")
      }
    }
    y = x[rownames(seqAnno), columnName]
    y[is.na(y)] = 0
    signal[ , i] = y
    if (!is.null(signalStart)){
      y = x[rownames(seqAnno), columnNameStart]
      y[is.na(y)] = 0
      signalStart[ , i] = y
    }
    if (!is.null(signalEnd)){
      y = x[rownames(seqAnno), columnNameEnd]
      y[is.na(y)] = 0
      signalEnd[ , i] = y
    }
  }
  
  if (ezIsSpecified(param$useTranscriptType)){
    use = seqAnno$type == param$useTranscriptType
  } else {
    use = TRUE
  }
  
  if (param$useSigThresh){
    sigThresh = param$sigThresh
  } else {
    sigThresh = 0
  }
  
  rawData = list(counts=signal[use, ,drop=FALSE], countsStart=signalStart[use, ,drop=FALSE], 
                 countsEnd=signalEnd[use, ,drop=FALSE], isLog=FALSE,
                 presentFlag=signal[use, ,drop=FALSE] > sigThresh, 
                 seqAnno=seqAnno[use, , drop=FALSE], featureLevel=dataFeatureLevel,
                 type="Counts", countName=columnName, dataset=input$meta)
  if (dataFeatureLevel == "isoform" && param$featureLevel == "gene"){
    rawData = aggregateCountsByGene(param, rawData)
  }
  rawData$rpkm = getRpkm(rawData)  
  rawData$tpm = getTpm(rawData)
  return(rawData)
}
