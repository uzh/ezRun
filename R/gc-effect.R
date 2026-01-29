ezGcEffect = function(rawData, param) {
  if (is.null(rawData$signal)) {
    rawData$signal = ezNorm(
      rawData$counts,
      presentFlag = rawData$presentFlag,
      method = param$normMethod
    )
  }
  isPresent = rowMeans(rawData$presentFlag) >= 0.5
  logSig = log2(rawData$signal + param$backgroundExpression)
  logRatio = logSig - rowMeans(logSig)
  isHighGc = rawData$seqAnno$gc > 0.6 & isPresent
  isLowGc = rawData$seqAnno$gc < 0.4 & isPresent
  gcEffect = apply(logRatio, 2, function(x) {
    mean(x[isHighGc], na.rm = TRUE) - mean(x[isLowGc], na.rm = TRUE)
  })
  return(gcEffect)
}
