


cleanupTwoGroupsInput = function(input, param){
  dataset = input$meta
  if (param$useFactorsAsSampleName){
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(ezDesignFromDataset(dataset, param), 1, paste, collapse="_"))
  }
  if (!is.null(param$removeOutliers) && param$removeOutliers && !is.null(dataset$Outlier)){
    dataset = dataset[toupper(dataset$Outlier) %in% c("", "NO", '""', "FALSE") == TRUE, ]
  }
  inputMod = EzDataset(meta=dataset)
  if (!is.null(param$markOutliers) && param$markOutliers){
    stopifnot(!is.null(dataset$Outlier))
    grouping = inputMod$getColumn(param$grouping)
    isOut = dataset$Outlier %in% c("", "NO", '""', "FALSE") == FALSE
    grouping[isOut] = paste(grouping[isOut], "OUTLIER", sep="_")
    inputMode$setColumn(grouping)
  }
  retrurn(inputMod)
}
