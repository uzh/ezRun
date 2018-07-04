###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


doEnrichr = function(param){
  grepl("^Homo_|^Rattus_|^Sus_|^Canis_|^Mus_", getOrganism(param$ezRef)) & param$featureLevel == "gene"
}

##' @title Run Enrichr
##' @description Runs Enrichr on the specified genes against all libraries
##'
##' If you find that this function is not flexible enough, you can run Enrichr as follows.
##'
##' \itemize{
##' \item Update the internal list of the libraries (if necessary) by running \code{\link{parseEnrichrLibNames}}
##' \item Register gene list with \code{\link{enrichrAddList}}
##' \item Run enrichment analysis with \code{\link{enrichrEnrich}}
##' \item Check the length \code{failure} element of the response to see if there were any failures
##' \item Filter the \code{success} element of the response by using \code{\link{filterEnrichrResults}}
##' }
##'
##' @param genes A list or vector containing gene names
##' @param minScore Combined score (\eqn{log(AdjP) * z}) threshold to filter the results.
##' @param maxAdjP Adjusted p (Benjamini-Hochberg method) threshold to filter the results.
##' @param minOverlapGenes Minimum number of genes in the overlap
##' @param softFilter select libraries with at least one significant result but report all results
##' @param maxResult maximum number of results to return
##' @param connectionN Maximum number of connections to open simultaneously for the enrichment analysis.
##' @return A named \code{list} with library names as names and \code{data.frame} objects as values.
##' An entire result set is retained if there is at least one significant gene based on the supplied
##' criteria.
##'
##' @examples
##' \dontrun{
##' genes <- c("HPS1", "HPS3", "HPS5", "PSMD11", "CUL7", "SNAPC4", "GNB2")
##' resList <- runEnrichr(genes, minScore = 15, maxAdjP = 1e-4, minOverlapGenes = 3, connectionN = 20)
##' resList = lapply(names(resList), function(nm){return(cbind("Gene-set library"=nm, resList[[nm]]))}) ## add the name as a first column
##' resMerged = do.call("rbind", resList)
##' }
##'
##' # Alternative usage
##' \dontrun{
##' addResp <- enrichrAddList(genes)
##' names(addResp)
##' # [1] "shortId"    "userListId"
##' # This will check all libraries specified in the internal list (extdata/enrichr_libnames.txt)
##' res <- enrichrEnrich(addResp$userListId)
##' # Check if any requests failed
##' if (length(res$failure) > 0) {
##'   cat(length(res$failure), " requests failed")
##'   allLibs <- getEnrichrLibNames()
##'   # Libraries that failed
##'   setdiff(allLibs, names(res$success))
##' }
##' # Filter the data sets
##' resFilt <- filterEnrichrResults(res$success, maxPAdj = 0.01, maxZ = -3, minCombinedScore = 12, minOverlapGenes = 3)
##' }
##' @author Roman Briskine
runEnrichr <- function(genes, minScore = 12, maxAdjP = 0.01, minOverlapGenes = 3, softFilter = F, 
                       maxResult = NA, connectionN = 10) {
  geneListResp <- enrichrAddList(genes)
  res <- enrichrEnrich(geneListResp$userListId, connectionN = connectionN)
  failureN <- length(res$failure)
  if (failureN > 0) {
    message("There were ", failureN, " failures while running Enrichr. First is displayed below.")
    message(res$failure[[1]])
  }
  resFilt <- list()
  if (length(res$success) > 0) {
    resFilt <- filterEnrichrResults(res$success, minCombinedScore = minScore, maxPAdj = maxAdjP, 
                                    minOverlapGenes = minOverlapGenes, softFilter = softFilter,
                                    maxResult = maxResult)
  }
  resFilt
}


##' @title Add gene list
##' @description Registers a gene list with the Enrichr server
##' @param genes A list or vector containing gene names
##' @template roxygen-template
##' @return Returns a list with userListId and shortId
##' @seealso \code{\link{runEnrichr}}
##' @author Roman Briskine
enrichrAddList <- function(genes) {
  require(httr, quietly = T, warn.conflicts = WARN_CONFLICTS)

  reqUrl <- paste(ENRICHR_BASE_URL, "addList", sep = "/")

  if (is.null(genes) || length(genes) < 1) {
    stop("The gene list should contain at least one gene")
  }

  geneStr <- paste(genes, collapse = "\n")
  resp <- POST(reqUrl, body = list(list = geneStr))
  if (http_error(resp)) {
    stop_for_status(resp, "register the gene list with Enrichr")
  }
  respParsed <- jsonlite::fromJSON(httr::content(resp, as = "text"))

  if (is.null(respParsed$userListId) || is.null(respParsed$shortId)) {
    stop("Enrichr server returned an invalid response: userListId or shortId is missing")
  }

  return(respParsed)
}

##' @title Query string generator
##' @description Helper function that generates a query string out of named parameters. The main
##' purpose of this function is to generate query strings for curl requests. All special characters
##' in the values are escaped.
##' @param ...  Named parameters that will form the query
##' @template roxygen-template
##' @return A string in the following format: \code{param1=value1&param2=value2}.
##' @author Roman Briskine
mkCurlQryString <- function(...) {
  x <- c(...)
  paste(paste(names(x), curl::curl_escape(x), sep = "="), collapse = "&")
}


##' @title Enrichment analysis
##' @description Runs enrichment analysis against the specified databases. By default, the function
##' will use the internal list (\code{\link{getEnrichrLibNames}}) of the library names.
##' @param userListId userListId returned by the Enrichr server via \code{\link{enrichrAddList}}
##' @param libNames vector of library names that should be a subset of those specified at the
##'   \href{http://amp.pharm.mssm.edu/Enrichr/#stats}{Enrichr} website. This function does not
##'    validate the list. If the value is NULL, the function will use the internal list.
##' @param connectionN maximum number of concurrent connections to make
##' @param retryN total number of times to retry a query if it fails
##' @template roxygen-template
##' @return \code{list} with two elements: success and failure. The first element is a \code{list}
##' of \code{data.frame} objects that contain enrichment results per library. The names of the
##' success list match the library names. The failure \code{list} contains error messages for the
##' libraries that failed enrichment tests. The name of this \code{list} may not match the library
##' names. You can determine which libraries failed by comparing the success list with the supplied
##' library names.
##' @seealso \code{\link{runEnrichr}, \link{enrichrAddList}, \link{getEnrichrLibNames}}
##' @author Roman Briskine
enrichrEnrich <- function(userListId, libNames = getEnrichrLibNames(), connectionN = 10, retryN = 3) {
  # While httr is easier to deal with, it does not support asynchrous requests. Neither does Rcurl.
  # So, we have to resort to the use of the curl package.
  require(curl, quietly = T, warn.conflicts = WARN_CONFLICTS)

  reqUrl <- paste(ENRICHR_BASE_URL, "enrich", sep = "/")

  # The functions below are internal helper functions to make the code clear. There is no reason to
  # expose them.

  # Concatenates gene list into a string. Otherwise, the unlist function would expand them. Note
  # that the recursive flag will not help you in this case.
  concatGenes <- function(lst, sep = ";") {
      lapply(lst, function(x) {
        x[[10]] = length(x[[6]])
        x[[6]] <- paste(x[[6]], collapse = sep);
        x})
  }

  # The jsonlite parser returns all values as character() and you cannot call as.numeric on a subset
  # of columns.
  respToNumeric <- function(x, idx = c(3,4,5,7,8,9)) {
    for (i in idx) {
      x[,i] <- as.numeric(x[,i])
    }
    x
  }

  # Parses the response
  parseResp <- function(resp) {
    respParsed <- jsonlite::fromJSON( rawToChar(resp$content) )
    if ("expired" %in% names(respParsed) && respParsed$expired == T) {
      stop("The server response indicates that the gene set has expired")
    }
    respParsed <- respParsed[[1]]
    fieldNames <- c("Rank", "Term", "p_value", "z_score", "Combined.Score", "Overlapping.Genes",
                    "Adjusted.p_value", "Old.p_value", "Old.Adjusted.p_value", "nOverlapping.Genes")
    if (length(respParsed) > 0) {
      ds <- as.data.frame(
        do.call("rbind", lapply(concatGenes(respParsed), unlist)),
        stringsAsFactors = F
      )
      ds <- respToNumeric(ds)
    } else {
      ds <- as.data.frame(matrix(nrow = 0, ncol = length(fieldNames)))
    }
    names(ds) <- fieldNames
    ds
  }

  success <<- list()
  failure <<- list()
  tryId <- 0

  while ((tryId == 0 || (length(failure) > 0) && tryId < retryN + 1)) {
    tryId <- tryId + 1
    message(sprintf("Try: %d; failures: %d", tryId, length(failure)))
    failure <<- list()
    libNamesToQuery <- setdiff(libNames, names(success))
    # All connections are to the same host, so both parameters should have the same value
    pool <- new_pool(total_con = connectionN, host_con = connectionN, multiplex = T)
    for (libName in libNamesToQuery) {
      qryString <- mkCurlQryString(userListId = userListId, backgroundType = libName)
      
      curl_fetch_multi(
          paste(reqUrl, qryString, sep="?"),
          # failonerror=T will cause curl to fail for any response status >= 400
          handle = new_handle(failonerror = T),
          pool = pool,
          done = function(x) {
            # libName may be bound after the loop is done rather than immediately. In that case, all
            # results will be saved under the same libName. To avoid that, we have to extract the
            # libName here.
            saveRDS(x, file=paste0("enrichr-done-", ezTime(), "-debug.rds"))
            libName <- sub('^.+backgroundType=([^&]+).*$', '\\1', x$url)
            success[[libName]] <<- tryCatch(
              parseResp(x),
              error = function(e) {
                saveRDS(e, file=paste0("enrichr-error-", ezTime(), "-debug.rds"))
                failure[[libName]] <<- paste("Response parsing failure with", e) }
            )
          },
          # In case of a failure, curl returns only the error message that does not contain the URL,
          # so we have no way to determine which request failed. So, we'll just append the message to
          # the end of the list.
          fail = function(x) {
            saveRDS(x, file=paste0("enrichr-fail-", ezTime(), "-debug.rds"))
            k <- length(failure) + 1; failure[[k]] <<- x }
        )
    }
    multi_run(pool = pool, timeout=300)
  }
  ezSystem("rm -f enrich-*-debug.rds")
  list(success = success, failure = failure)
}


##' @title Result filtering
##' @description Filters the results returned by Enrichr. Libraries that fail to yield any
##' significant results are removed. However, libraries that contain significant results, retain all
##' records whether significant or not)
##' @param resList list of results, i.e. \code{success} field in the list returned by \code{\link{enrichrEnrich}}
##' @param maxP p-value threshold (maximum)
##' @param maxPAdj adjusted p-value threshold (maximum)
##' @param maxZ z-value threshold (maximum)
##' @param minCombinedScore combined score threshold (minimum)
##' @param minOverlapGenes minimum number of genes in the overlap
##' @param softFilter select libraries with at least one significant result but report all results
##' @param maxResult maximum number of results to return
##' @template roxygen-template
##' @return \code{list} that contains the results for libraries with at least one significant gene
##' that satisfies the criteria.
##' @seealso \code{\link{runEnrichr}, \link{enrichrEnrich}}
##' @author Roman Briskine
filterEnrichrResults <- function(resList, maxP = 1, maxPAdj = 1, maxZ = 0, minCombinedScore = 0, 
                                 minOverlapGenes=3, softFilter = F, maxResult = NA) {
  resF <- lapply(resList, function(x) {
    mask <- x$p_value <= maxP &
      x$Adjusted.p_value <= maxPAdj &
      x$z_score <= maxZ &
      x$Combined.Score >= minCombinedScore &
      x$nOverlapping.Genes >= minOverlapGenes
    maskIdx = which(mask)
    if (!is.na(maxResult) && length(maskIdx) > maxResult){
      maskIdx = maskIdx[1:maxResult]
    }
    x[maskIdx, ]
  })
  listMask <- sapply(resF, function(x) { nrow(x) > 0 }, simplify = T)
  if (softFilter) {
    resList[listMask]
  } else {
    resF[listMask]
  }
}


##' @title Get Enrichr library names
##' @description reads Enrichr library names from the internal file
##'   (\code{extdata/enrichr_libnames.txt})
##' @template roxygen-template
##' @return a vector containing library names
##' @seealso \code{\link{runEnrichr}}
##' @author Roman Briskine
getEnrichrLibNames <- function() {
  file <- system.file(file.path("extdata", "enrichr_libnames.txt"), package = "ezRun", mustWork = T)
  scan(file, character(), quiet = T)
}


##' @title Retrieve Enrichr library names
##' @description Retrieves Enrichr library names from the specified HTML file and saves it to the
##'   internal package file (\code{extdata/enrichr_libnames.txt}). Ideally, this method
##'   would fetch the data from \url{http://amp.pharm.mssm.edu/Enrichr/#stats}. However, the page
##'   uses an Ajax query to populate the table, so it is empty when you access it with libcurl.
##'   Note: this does not modify the file in your source directory but the file in your installed package
##' @param file location of the HTML file saved from \url{http://amp.pharm.mssm.edu/Enrichr/#stats}
##' @template roxygen-template
##' @return (invisibly) a vector containing library names
##' @seealso \code{\link{runEnrichr}}
##' @author Roman Briskine
parseEnrichrLibNames <- function(file) {
  require(XML, quietly = T, warn.conflicts = WARN_CONFLICTS)
  mainPage <- htmlParse(file)
  tblNode <- getNodeSet(mainPage, "//table[@id='stats']")
  tbl <- readHTMLTable(tblNode[[1]], as.data.table = T, stringsAsFactors=F)
  fileOut <- system.file(file.path("extdata", "enrichr_libnames.txt"), package = "ezRun", mustWork = T)
  write(tbl[, 1], file = fileOut, ncolumns = 1)
  invisible(tbl[, 1])
}
