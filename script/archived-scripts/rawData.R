getSignal = function(rawData){
  if (rawData$isLog){
    return(2^rawData$signal)
  } else {
    return(rawData$signal)
  }
}
