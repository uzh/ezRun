###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Expands genomic ranges
##' @description Expands genomic ranges of a GRanges object. All ranges will be expanded at both ends by \code{width}.
##' @param x an object of the class GRanges.
##' @param width an integer specifying how much to expand the ranges.
##' @template roxygen-template
##' @seealso \code{\link[GenomicRanges]{GRanges}}
##' @return Returns an expanded GRanges object.
##' @examples
##' x = GRanges(c("chr1", "chr1", "chr2"), IRanges(5000:5002, 8000:8002), strand=c("+", "-", "+"))
##' expandGRanges(x)
expandGRanges = function(x, width = 2000) {
  sl = seqlengths(x)
  sl[is.na(sl)] = Inf
  start(x) = pmax(start(x) - width, 1)
  end(x) = pmin(end(x) + width, sl[as.character(seqnames(x))], na.rm = TRUE)
  return(x)
}

##' @title Gets the values of genomic ranges
##' @description Gets the values of genomic ranges
##' @param cov an object of the class RleList.
##' @param targetRanges an object of the class GRanges.
##' @param doRev a logical vector indicating whether to reverse ranges. Usually derived from \code{targetRanges}.
##' @param asMatrix a logical indicating whether to return the values as a matrix.
##' @template roxygen-template
##' @seealso \code{\link[IRanges]{RleViewsList}}
##' @seealso \code{\link[IRanges]{viewApply}}
##' @return Returns the values or a matrix containing them.
getRangeValues = function(
  cov,
  targetRanges,
  doRev = as.character(strand(targetRanges)) == "-",
  asMatrix = TRUE
) {
  names(targetRanges) = 1:length(targetRanges)
  rgs = ranges(RangedData(targetRanges))
  gc()
  #stopifnot(setequal(names(cov), names(rgs)))
  targetViews = RleViewsList(rleList = cov[names(rgs)], rangesList = rgs) ###
  targetVal = viewApply(
    targetViews,
    function(x) {
      as.vector(x)
    },
    simplify = FALSE
  )
  idx = match(
    names(targetRanges),
    names(unlist(targetViews, use.names = FALSE))
  )
  values = unlist(targetVal)@listData[idx]
  if (any(doRev)) {
    values[doRev] = lapply(values[doRev], rev)
  }
  if (asMatrix) {
    profMatrix = t(as.matrix(as.data.frame(values)))
    rownames(profMatrix) = names(targetRanges)
    gc()
    return(profMatrix)
  } else {
    gc()
    return(values)
  }
}

##' @title Computes stats of genomic ranges
##' @description Computes stats of genomic ranges.
##' @param cov an object of the class RleList.
##' @param targetRanges an object of the class GRanges.
##' @param FUN a function to apply to the views.
##' @template roxygen-template
##' @seealso \code{\link[IRanges]{RleViewsList}}
##' @seealso \code{\link[IRanges]{viewApply}}
##' @return Returns the computed values belonging to their genomic ranges.
computeRangeStats = function(cov, targetRanges, FUN = mean) {
  rgs = ranges(RangedData(targetRanges))
  gc()
  targetViews = RleViewsList(rleList = cov[names(rgs)], rangesList = rgs) ###
  targetVal = viewApply(targetViews, FUN)
  values = unlist(targetVal)
  names(values) = names(targetRanges)
  return(values)
}

##' @title Subsets an Rle object
##' @description Subsets an Rle object.
##' @param x an object of the class Rle.
##' @param idx an integer vector specifying which indices to keep.
##' @template roxygen-template
##' @return Returns a subset of an Rle object as an integer or numeric vector.
##' @examples
##' rleobj = Rle(c(1, 1, 1, 2, 2, 3, 2, 2))
##' subSampleRle(rleobj, 1:8)
subSampleRle = function(x, idx) {
  return(as.vector(x[idx]))
}


.ezSequenceFromRanges = function(granges, chromSeqs) {
  stopifnot(unique(seqnames(granges)) %in% names(chromSeqs))

  ## TODO implement according to the example in p1432/funcs-v30.r
}

.getRleFromRanges = function(x, r) {
  gc()
  message("getRleFromRanges")
  rle = mapply(
    function(start, end, x) window(x, start = start, end = end),
    start(r),
    end(r),
    MoreArgs = list(x = x),
    SIMPLIFY = FALSE
  )
  names(rle) = names(r)
  isNeg = as.character(strand(r)) == "-"
  rle[isNeg] = lapply(rle[isNeg], rev)
  return(rle)
}

.getAvgProfileByOrientation = function(readRanges, targetRanges) {
  covPos = coverage(readRanges[as.character(strand(readRanges)) == "+"])
  covNeg = coverage(readRanges[as.character(strand(readRanges)) == "-"])
  targetProfilesPos = getRangeValues(covPos, targetRanges)
  targetProfilesNeg = getRangeValues(covNeg, targetRanges)
  targetStrand = as.character(strand(targetRanges))
  stopifnot(targetStrand %in% c("+", "-"))
  isPosTarget = targetStrand == "+"
  forw = (colSums(targetProfilesPos[isPosTarget, ]) +
    colSums(targetProfilesNeg[!isPosTarget, ])) /
    length(isPosTarget)
  rev = (colSums(targetProfilesPos[!isPosTarget, ]) +
    colSums(targetProfilesNeg[isPosTarget, ])) /
    length(isPosTarget)
  return(list(forw = forw, rev = rev))
}

.getProfileMatrixForRanges = function(targetRanges, cov, xPos, idx = NULL) {
  targetRangesList = split(targetRanges, seqnames(targetRanges))
  profiles = unlist(mapply(
    .getRleFromRanges,
    cov[names(targetRangesList)],
    targetRangesList,
    SIMPLIFY = FALSE
  ))
  names(profiles) = sub(".*?\\.", "", names(profiles))
  if (!is.null(idx)) {
    profMatrix = t(sapply(profiles, subSampleRle, idx)) ## TODO: idx should be defined before the if or be passed as an argument.
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

### Get the consensus peaks from a list of peaks
consensusPeaks <- function(x) {
  if (!is(x, "GRangesList")) {
    stop("x must be a GRangesList object.")
  }
  if (length(x) == 1L) {
    return(x[[1]])
  }
  if (length(x) > 2) {
    consensusPeaks(c(
      GRangesList(consensusPeaks(x[1:(length(x) - 1)])),
      x[length(x)]
    ))
  }
  gr1 <- x[[1]]
  gr2 <- x[[2]]
  hits <- findOverlaps(gr1, gr2)
  grPairs <- Pairs(first = gr1, second = gr2, hits = hits)
  intersectedGR <- pintersect(first(grPairs), second(grPairs))
  toKeep <- width(intersectedGR) / width(first(grPairs)) > 0.5 |
    width(intersectedGR) / width(second(grPairs)) > 0.5
  return(first(grPairs)[toKeep])
}
