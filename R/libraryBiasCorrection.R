


ezCorrectBias = function(counts, gc, width,
                       gcCore = c(0.45, 0.50),
                       widthCore = c(10, 10.5),
                       widthOffset = 200,
                       minCount = 3,
                       minGenesInBin = 20,
                       logWidthBreaks = c(9.5, 10.5, 11.5),
                       gcBreaks = c(0.4, 0.43, 0.46, 0.49, 0.53, 0.57, 0.61, 0.64, 0.67),
                       quantileValue = 0.5,
                       minPresFraction=0.5){
  
  nGenes = nrow(counts)
  nSamples = ncol(counts)
  if (is.null(colnames(counts))){
    colnames(counts) = paste("sample", 1:ncol(counts), sep="_")
  }
  
  logWidth = log2(width + widthOffset)
  
  ###### normalize with respect to sequencing depth
  useForNorm = gc >= gcCore[1] & gc <= gcCore[2] & logWidth >= widthCore[1] & logWidth <= widthCore[2]
  useForNorm = useForNorm & apply(counts > minCount, 1, mean) > minPresFraction
  
  
  require(DESeq2)
  dds = DESeqDataSetFromMatrix(countData=round(counts), colData=data.frame(names=colnames(counts)), ~ names)
  dds = estimateSizeFactors(dds, controlGenes=useForNorm)
  scalingFactors = 1/dds@colData$sizeFactor

  logCountNorm = log(sweep(counts, 2, scalingFactors, FUN="*"))
  
  medLogCount = apply(logCountNorm, 1, median)
  logRatios = sweep(logCountNorm, 1, medLogCount)
  
  ## divide the genes in width - gc bins
  useForFit = is.finite(medLogCount)
  lwClasses = ezCut(logWidth, breaks = logWidthBreaks)
  gcClasses = ezCut(gc, gcBreaks)
  
  binMatrix = tapply((1:nGenes)[useForFit], list(logWidth=lwClasses[useForFit], gc=gcClasses[useForFit]), identity)
  binMatrixTemplate = matrix(0, nrow = nrow(binMatrix), ncol=ncol(binMatrix),
                             dimnames = list(rownames(binMatrix), colnames(binMatrix)))
  binCounts = binMatrixTemplate
  for (i in 1:length(binMatrix)){
    binCounts[i] = length(binMatrix[[i]])
  }
  ## we recompute the bin centers; this is needed for the extreme bins because they cover half-open intervals; 
  ## and is done for consistency reasons for all bins
  gcMids = tapply(gc, gcClasses, mean)
  gcmDef = (c(gcBreaks[1], gcBreaks) + c(gcBreaks, gcBreaks[length(gcBreaks)])) /2
  gcMids = ifelse(!is.na(gcMids), gcMids, gcmDef)
  binGcMids = binMatrixTemplate
  for (nm in colnames(binMatrixTemplate)){
    binGcMids[ ,nm] = gcMids[nm]
  }
  logWidthMids = tapply(logWidth, lwClasses, mean)
  lwmDef = (c(logWidthBreaks[1], logWidthBreaks) + c(logWidthBreaks, logWidthBreaks[length(logWidthBreaks)])) /2
  logWidthMids = ifelse(!is.na(logWidthMids), logWidthMids, lwmDef)
  binWidthMids = binMatrixTemplate
  for (nm in rownames(binMatrixTemplate)){
    binWidthMids[nm, ] = logWidthMids[nm]
  }
  
  ## shrink the values to the centers so that akima::interpp can get an interpolated value for all genes
  ## this is equivalent to extrapolate with a constant value
  gcShrinked = gc
  gcShrinked[gcShrinked > max(gcMids)] = max(gcMids)
  gcShrinked[gcShrinked < min(gcMids)] = min(gcMids)
  logWidthShrinked = logWidth
  logWidthShrinked[logWidthShrinked > max(logWidthMids)] = max(logWidthMids)
  logWidthShrinked[logWidthShrinked < min(logWidthMids)] = min(logWidthMids)
  
  ## now we fit 
  logOffset = counts
  logOffset[] = NA
  for (sm in colnames(counts)){
    binScores = binMatrixTemplate
    ## get the estimates for the bins
    for (i in 1:length(binMatrix)){
      values = logRatios[binMatrix[[i]], sm]
      binScores[i] = quantile(values, quantileValue)
      if (is.infinite(binScores[i])){
        if (sum(is.finite(values)) > 3){
          if (binScores[i] == Inf){
            binScores[i] = max(values[is.finite(values)])
          } else {
            binScores[i] = min(values[is.finite(values)])
          }
        } else {
          binScores[i] = NA
        }
      }
    }
    binScores[binCounts < minGenesInBin] = NA
    binScoresFilled = fillNaWithClosestNeighbor(binScores)
    
    ## 2d interpolation --------
    require(akima)
    fit2d = interpp(x=binGcMids, y=binWidthMids, z=binScoresFilled,
                    xo=gcShrinked, yo=logWidthShrinked, extrap = FALSE, linear=TRUE)
    logOffset[ ,sm] = - log(scalingFactors[sm]) + fit2d$z
  }
  
  correctedCounts = counts / exp(logOffset)
  hasNoCorr = !is.finite(correctedCounts)
  correctedCounts[hasNoCorr] = counts[hasNoCorr]
  return(list(correctedCounts=correctedCounts, logOffset=logOffset, nGenesForNorm=sum(useForNorm), nGenesForFit=sum(useForFit),
              gcShrinked=gcShrinked, logWidthShrinked=logWidthShrinked, binScoresFilled=binScoresFilled))
}



fillNaWithClosestNeighbor = function(x){
  rowPosMatrix = matrix(rep(1:nrow(x), ncol(x)), nrow=nrow(x), ncol=ncol(x))
  colPosMatrix = matrix(rep(1:ncol(x), each=nrow(x)), nrow=nrow(x), ncol=ncol(x))
  isKnown = !is.na(x)
  for (i in 1:nrow(x)){
    for (j in 1:ncol(x)){
      if (!isKnown[i,j]){
        d = (i - rowPosMatrix)^2 + (j - colPosMatrix)^2
        d[!isKnown] = Inf
        isNearest = d == min(d)
        x[i,j] = mean(x[isNearest])
      }
    }
  }
  return(x)
}
