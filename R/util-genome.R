###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Make a coordinate
##' @description Creates a chromosome coordinate with the provided input.
##' @param chrom the chromosome name.
##' @param strand the strand of the chromosome.
##' @param start where the coordinate starts.
##' @param stop where the coordinate stops.
##' @return Returns a chromosome coordinate.
##' @template roxygen-template
##' @examples
##' makeCoordinate("chrm", 3, 150, 300)
makeCoordinate = function(chrom, strand, start, stop) {
  paste0(chrom, "(", strand, ")", ":", start, "-", stop)
}

##' @title Splits fragment coordinates
##' @description Use this function to split chromosome range information in the format of chromosome:start-stop.
##' @param pos a character vector of chromosome ranges to split.
##' @return Returns a data.frame containing the chromosome name, the start and the stop position.
##' @template roxygen-template
##' @examples
##' splitCoordinate("chrm:1-50")
splitCoordinate = function(pos) {
  chrom = sub(":.*", "", pos)
  start = as.numeric(sub("-.*", "", sub(".*:", "", pos)))
  stop = as.numeric(sub(".*-", "", sub(".*:", "", pos)))
  result = data.frame(
    Chrom = chrom,
    Start = start,
    Stop = stop,
    row.names = pos,
    stringsAsFactors = FALSE
  )
  return(result)
}

##' @describeIn splitCoordinate Does the same as \code{splitCoordinate}, but returns the result as a list and is more robust.
splitRegion = function(pos) {
  if (grepl(":", pos)) {
    chrom = sub(":.*", "", pos)
    stopifnot(grepl("-", pos))
    se = sub(".*:", "", pos)
    start = as.integer(sub("-.*", "", se))
    end = as.integer(sub(".*-", "", se))
  } else {
    chrom = pos
    start = NA
    end = NA
  }
  return(list(seq = chrom, start = start, end = end))
}
