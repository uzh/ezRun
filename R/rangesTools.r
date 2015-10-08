###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezSequenceFromRanges = function(granges, chromSeqs){
  stopifnot(unique(seqnames(granges)) %in% names(chromSeqs))
  
  ## TODO implement according to the example in p1432/funcs-v30.r
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
expandGRanges = function(x, width=2000){
  
  sl = seqlengths(x)
  sl[is.na(sl)] = Inf
  start(x) = pmax(start(x) - width, 1)
  end(x) = pmin(end(x) + width, sl[as.character(seqnames(x))], na.rm=TRUE)
  return(x)
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
ezCountOverlaps = function(query, subject, ignoreStrand=FALSE, ...){
  #qSeqs = levels(seqnames(query))
  #sSeqs = levels(seqnames(subjects))
  #if (length(setdiff(qSeqs, sSeqs)) > 0 & length(setdiff(sSeqs, qSeqs)) > 0){
  #  warning("incompatible seqnames!\nqSeqs: ", paste(qSeqs, collapse=" "),
  #       "\nsSeqs: ", paste(sSeqs, collapse=" "))
  #}
  if (ignoreStrand){
    if ( class(query) == "GRanges"){
      strand(query) = "*"  
    } else {
      sr = unlist(query)
      strand(sr) = "*"
      nms = sub("\\..*$", "", names(sr))
      query = split(sr, nms)      
    }
    if ( class(subject) == "GRanges"){
      strand(subject) = "*"  
    } else {
      sr = unlist(subject)
      strand(sr) = "*"
      nms = sub("\\..*$", "", names(sr))
      subject = split(sr, nms)      
      #foo = lapply(elementLengths(subject), function(x){Rle("*", x)})
      #s = RleList(foo)
      #strand(subject) = s
      #for (nm in names(subject)){
      #  message(nm)
      #  strand(subject[[nm]]) = "*"
      #}
    }
  }
  return(countOverlaps(query=query, subject=subject, ...))
}


##' @describeIn ezCountOverlaps
ezHasOverlaps = function(query, subject, ignore.strand=FALSE, ...){  
  return(countOverlaps(query=query, subject=subject, ignore.strand=ignore.strand, ...) > 0)
}


# TODOP: purpose? seems identical to as.numeric()
ezViewValues = function(x){
  return(as.numeric(x))
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
getRangeValues = function(cov, targetRanges, doRev=as.character(strand(targetRanges)) == "-", asMatrix=TRUE){
  names(targetRanges) = 1:length(targetRanges)
  rgs = ranges(RangedData(targetRanges))
  gc()
  #stopifnot(setequal(names(cov), names(rgs)))
  targetViews = RleViewsList(rleList=cov[names(rgs)], rangesList=rgs)  ###
  targetVal = viewApply(targetViews, function(x){as.vector(x)}, simplify=FALSE)
  idx = match(names(targetRanges), names(unlist(targetViews, use.names=FALSE)))
  values = unlist(targetVal)@listData[idx]
  if (any(doRev)){
    values[doRev] = lapply(values[doRev], rev)
  }
  if (asMatrix){
    profMatrix  = t(as.matrix(as.data.frame(values)))
    rownames(profMatrix) = names(targetRanges)
    gc()
    return(profMatrix)
  } else {
    gc()
    return(values)
  }
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
computeRangeStats = function(cov, targetRanges, FUN=mean){
  rgs = ranges(RangedData(targetRanges))
  gc()
  targetViews = RleViewsList(rleList=cov[names(rgs)], rangesList=rgs)  ###
  targetVal = viewApply(targetViews, FUN)
  values = unlist(targetVal)
  names(values) = names(targetRanges)
  return(values)
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
getRleFromRanges = function(x, r){
  gc()
  message("getRleFromRanges")
  rle = mapply(function(start, end, x) window(x, start=start, end=end),
               start(r), end(r), MoreArgs=list(x=x), SIMPLIFY = FALSE )
  names(rle) = names(r)
  isNeg = as.character(strand(r)) == "-"
  rle[isNeg] = lapply(rle[isNeg], rev)
  return(rle)
}


##' @title 1
##' @description 1
##' @param x
##' @template roxygen-template
##' @return Returns
##' @examples
##' 1
subSampleRle = function(x, idx){
  return(as.vector(x[idx]))
}





.getAvgProfileByOrientation = function(readRanges, targetRanges){
  
  covPos = coverage(readRanges[as.character(strand(readRanges)) == "+"])
  covNeg = coverage(readRanges[as.character(strand(readRanges)) == "-"])
  targetProfilesPos = getRangeValues(covPos, targetRanges)
  targetProfilesNeg = getRangeValues(covNeg, targetRanges)
  targetStrand = as.character(strand(targetRanges))
  stopifnot(targetStrand %in% c("+", "-"))
  isPosTarget =  targetStrand == "+"
  forw = (colSums(targetProfilesPos[isPosTarget, ]) + colSums(targetProfilesNeg[!isPosTarget, ])) /length(isPosTarget)
  rev = (colSums(targetProfilesPos[!isPosTarget, ]) + colSums(targetProfilesNeg[isPosTarget, ])) /length(isPosTarget)
  return(list(forw=forw, rev=rev))
}

.getProfileMatrixForRanges = function(targetRanges, cov, xPos){
  targetRangesList = split(targetRanges, seqnames(targetRanges))
  profiles = unlist(mapply(getRleFromRanges, cov[names(targetRangesList)],
                           targetRangesList, SIMPLIFY=FALSE))
  names(profiles) = sub(".*?\\.", "", names(profiles))
  if (!is.null(idx)){
    profMatrix = t(sapply(profiles, subSampleRle, idx))
  } else {
    ## xPos have to be quantiles
    profLength = sapply(profiles, length)
    idxFrame = as.data.frame(round(xPos %o% profLength))
    profMatrix = t(mapply(subSampleRle, profiles, idxFrame))
    
  }
  colnames(profMatrix) = as.character(xPos)
  gc()
  return(profMatrix)
}
