###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Writes an html report
##' @description Writes an html report and returns the connection to it.
##' @template htmlFile-template
##' @param param a list of parameters to extract the \code{projectId} from.
##' @param title a character specifying the title of the html report.
##' @template dataset-template
##' @template roxygen-template
##' @return Returns a connection to the written html report.
##' @seealso \code{\link{writeJavaScriptIgvStarter}}
##' @seealso \code{\link{writeTxtLinksToHtml}}
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' openHtmlReport(htmlFile,param)
openHtmlReport = function(htmlFile, param=param, title="", dataset=NULL){
	html = file(htmlFile, open="wt")
	
	ezWrite("<html><head>", con=html)
	ezWrite("<!-- ", getwd(), " -->", con=html)
	writeLines(readLines(ezCSSFile()), con=html)
	file.copy(ezBannerFile(), ".")
  flush(html)
  ezWrite("<title>", title, "</title>", con=html)
	
	writeJavaScriptIgvStarter(htmlFile, param$projectId, html)
	
	ezWrite("</head><body>", con=html)
  
  ## formatting
  ezWrite('<div id="margins">', con=html)

	ezWrite('<table id="screenhead" bgcolor="#000000">', con=html)
  ezWrite('<tr>', con=html)
  ezWrite('<td align="left" class="identity-image" valign="bottom" width="800"><img src="', basename(ezBannerFile()), '" alt="" usemap="#logo" height="70px" border="0"></td><td align="left" valign="bottom" width="500">', con=html)
  ezWrite('<td align="left">&#160;</td>', con=html)
  ezWrite('</tr>', con=html)
  ezWrite('</table>', con=html)
  ezWrite('<div class="pathline">', con=html)
  ezWrite('<a href="http://www.fgcz.ethz.ch" title="FGCZ Home">Functional Genomics Center Zurich</a>', con=html)
  ezWrite('</div>', con=html)
  ezWrite('', con=html)
	ezWrite("Started on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " -- <a href='", DOC_URL, "'>Documentation</a>", con=html);
	ezWrite('<table>', con=html)
  ezWrite('<tr>', con=html)
  ezWrite('<td valign="top" align="left" id="contentblock">', con=html)
  ezWrite('<div class="listing">', con=html)
  
	ezWrite("<h1>", title, "</h1>", con=html)
	flush(html)
	if (!is.null(dataset)){
    ezWrite.table(dataset, file="dataset.tsv", head="Name")
    writeTxtLinksToHtml(txtNames = "dataset.tsv", con=html)
	}
	flush(html)
	return(html)
}

##' @describeIn openHtmlReport Returns a link to the report css file if it exists.
ezCSSFile = function(){
  if (REPORT_CSS_FILE == ""){
    return(system.file("extdata/style.css", package="ezRun", mustWork = TRUE))
  } else {
    if (file.exists(REPORT_CSS_FILE)){
      return(REPORT_CSS_FILE)
    } else {
      stop("file does not exist: ", REPORT_CSS_FILE)
    }
  }
}

##' @describeIn openHtmlReport Returns a link to the banner png file if it exists.
ezBannerFile = function(){
  if (REPORT_BANNER_FILE == ""){
    return(system.file("extdata/banner.png", package="ezRun", mustWork = TRUE))
  } else {
    if (file.exists(REPORT_BANNER_FILE)){
      return(REPORT_BANNER_FILE)
    } else {
      stop("file does not exist: ", REPORT_BANNER_FILE)
    }
  }
}

##' @title Closes and wraps up an html connection
##' @description Closes and wraps up an html connection adding the closing tags and a "Finished..." line.
##' @param html a connection to an html file.
##' @template roxygen-template
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' html =  openHtmlReport(htmlFile,param)
##' closeHTML(html)
closeHTML = function(html){
 
 ## formatting
 ezWrite("</div></td></tr></table>", con=html)
 ezWrite("Finished ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), con=html);
 ezWrite("</div></body></html>", con=html)
 flush(html)
 close(html)
}

##' @title Writes an error message
##' @description Writes an error message to an html file. Also creates the file and closes it.
##' @template htmlFile-template
##' @param param a list of parameters to extract the \code{projectId} and \code{name} from.
##' @template dataset-template
##' @param error a character vector representing the error message(s).
##' @template roxygen-template
##' @seealso \code{\link{openHtmlReport}}
##' @seealso \code{\link{closeHTML}}
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' writeErrorHtml(htmlFile, param)
writeErrorHtml = function(htmlFile, param=param, dataset=NULL, error="Unknown Error"){
	html = openHtmlReport(htmlFile, param=param,
	       title=paste("Error:", param$name), dataset=dataset)
	ezWrite("<h2>Error message</h2>", con=html)
	ezWrite(paste("<p>", error, "</p>", collapse="\n"), con=html)
	flush(html)
  closeHTML(html)
}

##' @title Writes an image table
##' @description Writes an image table to an html file, usually containing plot windows.
##' @param pngMatrix a matrix containing picture files.
##' @template connection-template
##' @template roxygen-template
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' html =  openHtmlReport(htmlFile,param)
##' writeImageTableToHtml(as.matrix("example"),html)
##' writeImageRowToHtml("example",html)
##' writeImageColumnToHtml("example",html)
writeImageTableToHtml = function(pngMatrix, con=stdout()){
  ezWrite("<table>", con=con)
  if (nrow(pngMatrix) > 0 && ncol(pngMatrix) > 0){
    for (i in 1:nrow(pngMatrix)){
      ezWrite("<tr>", con=con)
      for (j in 1:ncol(pngMatrix)){
        ezWrite(paste0("<td><img src='", pngMatrix[i, j], "'></td>"), con=con)
      }
      ezWrite("</tr>", con=con)
    }
  }
  ezWrite("</table>", con=con)
}

##' @describeIn writeImageTableToHtml Write the pictures in a row instead.
writeImageRowToHtml = function(pngNames, con=stdout()){
  writeImageTableToHtml(t(as.matrix(pngNames)), con=con)
}

##' @describeIn writeImageTableToHtml Write the pictures in a column instead.
writeImageColumnToHtml = function(pngNames, con=stdout()){
  writeImageTableToHtml(as.matrix(pngNames), con=con)
}

##' @title Writes a text link
##' @description Writes a text link to an html file.
##' @param txtNames a character representing the link name.
##' @param mime a character representing the type of the link.
##' @template connection-template
##' @template roxygen-template
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' html =  openHtmlReport(htmlFile,param)
##' writeTxtLinksToHtml("dataset.tsv", con=html)
writeTxtLinksToHtml = function(txtNames, mime="text/plain", con=stdout()){
  ezWrite("<table>", con=con)
  ezWrite(paste0("<tr><td><a href='", txtNames, "' type='", mime, "'>", txtNames, "</a></td></tr>", collapse="\n"), con=con)
  ezWrite("</table>", con=con)
}

##' @title Writes a table
##' @description Writes a table to an html file.
##' @param x a matrix or data.frame to paste a table from.
##' @template connection-template
##' @param bgcolors a matrix specifying the background colors in html format.
##' @param valign a character specifying with an html tag where to align the table vertically.
##' @param border an integer or numeric specifying the border width.
##' @param head a character specifying the header of the table.
##' @template roxygen-template
##' @examples
##' x = matrix(1:25,5)
##' rownames(x) = letters[1:5]
##' colnames(x) = LETTERS[1:5]
##' param = ezParam()
##' htmlFile = "example_html"
##' html =  openHtmlReport(htmlFile,param)
##' writeTableToHtml(x,html)
writeTableToHtml = function(x, con=stdout(), bgcolors=matrix("#ffffff", nrow=nrow(x), ncol=ncol(x)),
   valign="middle", border=1, head=""){

  ezWrite("<table border='", border, "'><tr>", con=con)
  ezWrite(paste0("<th>", c(head, colnames(x)), "</th>", collapse="\n"), con=con)
  ezWrite("</tr>", con=con)
  if (nrow(x) > 0){
    for (i in 1:nrow(x)){
      ezWrite("<tr><th>", rownames(x)[i], "</th>", con=con)
      ezWrite(paste0("<td valign='", valign, "' bgcolor='", bgcolors[i,], "'>", x[i,], "</td>", collapse="\n"), con=con)
      ezWrite("</tr>", con=con)
    }
  }
  ezWrite("</table>", con=con)
}


##' @title Writes a count result summary
##' @description Writes a count result summary to an html file.
##' @param html a connection to an html file.
##' @param param a list of parameters to influence the output:
##' \itemize{
##'  \item{batch}{ a logical indicating whether the second factor was used.}
##'  \item{comparison}{ which comparison was used.}
##'  \item{normMethod}{ the normalization method.}
##'  \item{sigThresh}{ the threshold...}
##'  \item{useSigThresh}{ ...and whether it should be used.}
##' }
##' @template result-template
##' @template roxygen-template
##' @examples
##' 1
writeCountResultSummary = function(html, param, result){

  ezWrite("<h2>Result Summary</h2>", con=html)

	ezWrite("<table border='0'>", con=html)
  ezWrite("<tr><td>Analysis:</td><td>", result$analysis, "</td></tr>", con=html)
  ezWrite("<tr><td>Feature level:</td><td>", result$featureLevel, "</td></tr>", con=html)
  ezWrite("<tr><td>Data Column Used:</td><td>", result$countName, "</td></tr>", con=html)
  ezWrite("<tr><td>Method:</td><td>", result$method, "</td></tr>", con=html)
  if (ezIsSpecified(param$batch)){
    ezWrite("<tr><td>Statistical Model:</td><td>", "used provided second factor", "</td></tr>", con=html)    
  }
  ezWrite("<tr><td>Comparison:</td><td>", param$comparison, "</td></tr>", con=html)
  if (!is.null(param$normMethod)){
    ezWrite("<tr><td>Normalization:</td><td>", param$normMethod, "</td></tr>", con=html)
  }
  ezWrite("<tr><td>Number of features:</td><td>", length(result$pValue), "</td></tr>", con=html)
  if (!is.null(result$isPresentProbe)){
    ezWrite("<tr><td>Number of features with counts above threshold:</td><td>", sum(result$isPresentProbe), "</td></tr>", con=html)
  }

  if (param$useSigThresh){
     ezWrite("<tr><td>Log2 signal threshold:</td><td>", signif(log2(param$sigThresh), digits=4), "</td></tr>", con=html)
     ezWrite("<tr><td>Linear signal threshold:</td><td>", signif(param$sigThresh, digits=4), "</td></tr>", con=html)
  }
  ezWrite("</table>", con=html)
}

## will get deprecated eventually
writeResultCounts = function(html, param, result, geneIds=NULL, pThresh=c(0.1, 0.05, 1/10^(2:5))){

	sigTable = getSignificantCountsTable(result, pThresh=pThresh)
	sigFcTable = getSignificantFoldChangeCountsTable(result, pThresh=pThresh)
	
  ezWrite("<p>Feature counts by significance and fold-change (fc)</p>", con=html)
  ezWrite("<table border=0><tr><td>", con=html)
  writeTableToHtml(sigTable, con=html)
  ezWrite("</td><td>", con=html)
  writeTableToHtml(sigFcTable, con=html)
  ezWrite("</td></tr></table>", con=html)
}

##' @title Writes a result file
##' @description Writes a result file in text format or zipped.
##' @param html a connection to an html file.
##' @param param a list of parameters that pastes the \code{comparison} into the file name and does a zip file if \code{doZip} is true.
##' @template result-template
##' @template rawData-template
##' @param useInOutput a logical specifying whether to use most of the result information.
##' @param file a character representing the name of the result file.
##' @template roxygen-template
##' @return Returns the name of the result file.
##' @seealso \code{\link{writeTxtLinksToHtml}}
writeResultFile = function(html, param, result, rawData, useInOutput=TRUE,
  file=paste0("result--", param$comparison, ".txt")){

  seqAnno = rawData$seqAnno
  probes = names(result$pValue)[useInOutput]
  y = data.frame(row.names=probes, stringsAsFactors=FALSE, check.names=FALSE)
  y[ , colnames(seqAnno)] = sapply(seqAnno[match(probes, rownames(seqAnno)), ], as.character)
  y$"log2 Signal" = result$log2Expr[useInOutput]
	y$"isPresent" = result$isPresentProbe[useInOutput]
  y$"log2 Ratio" = result$log2Ratio[useInOutput]
  y$"gfold (log2 Change)" = result$gfold[useInOutput]
  y$"log2 Effect" = result$log2Effect[useInOutput]
  y$"probesetCount" = result$nProbes[useInOutput]
  y$"presentProbesetCount" = result$nPresentProbes[useInOutput]
  y$ratio = result$ratio[useInOutput]
  y$pValue = result$pValue[useInOutput]
  y$fdr = result$fdr[useInOutput]
	for (nm in grep("Tukey pValue", names(result), value=TRUE)){
		y[[nm]] = result[[nm]][useInOutput]
	}
  if (!is.null(result$groupMeans)){
    groupMeans = result$groupMeans[useInOutput, ]
    colnames(groupMeans) = paste("log2 Avg of", colnames(groupMeans))
    y = data.frame(y, groupMeans, check.names=FALSE, stringsAsFactors=FALSE)
  }
  # 	if (!is.null(result$resVar)){
  # 		resVar = result$resVar[useInOutput, ]
  # 		colnames(resVar) = paste("Res. Var. of", colnames(resVar))
  # 	  y = data.frame(y, resVar, check.names=FALSE, stringsAsFactors=FALSE)
  # 		y$"Variance Outlier" = result$varianceOutlier
  # 	}
  
  if (!is.null(result$xNorm)){
    yy = result$xNorm[useInOutput, ]
    colnames(yy) = paste(colnames(yy), "[normalized count]")
    y = cbind(y, yy)
  }
  yy = getRpkm(rawData)[useInOutput, ]
  if (!is.null(yy)){
    colnames(yy) = paste(colnames(yy), "[FPKM]")
    y = cbind(y, yy)
  }
  y = y[order(y$fdr, y$pValue), ]
  ezWrite.table(y, 
     file=file, head="Identifier", digits=4)
  if (param$doZip){
    zipLink = zipFile(file)
    if (!is.null(html)){
      writeTxtLinksToHtml(zipLink, mime="application/zip", con=html)
      flush(html)
    }
  } else {
    if (!is.null(html)){
      writeTxtLinksToHtml(file, mime="application/text", con=html)
      flush(html)
    }
  }
  return(list(resultFile=file))
}


##' @title Writes QC scatter plots
##' @description Writes QC scatter plots to an html file.
##' @param html a connection to an html file.
##' @param param a list of parameters. If \code{writeScatterPlots} is false, the function returns NULL.
##' @param design a data.frame containing the factorial design.
##' @param conds a named character vector containing the conditions of the factorial design.
##' @template rawData-template
##' @param signalCond a set of values containing signals averaged by the conditions.
##' @param isPresentCond Either NULL or a data.frame containing coloring information.
##' @template colors-template
##' @template types-template
##' @template roxygen-template
writeQcScatterPlots = function(html, param, design, conds,
															 rawData, signalCond, isPresentCond, seqAnno,
															 colors=getBlueRedScale(), types=NULL){
  samples = rownames(design)
	nConds = length(unique(conds))
	signal = getSignal(rawData)
	signal[signal <= 0] = NA
  isPresent = ezPresentFlags(signal, presentFlag=rawData$presentFlag, param=param, isLog=rawData$isLog)
	signalRange = range(signal, na.rm=TRUE)

  ezWrite("<h2>Scatter Plots by Conditions</h2>", con=html)
  if (!is.null(rawData$seqAnno$gc)){
    gcTypes = data.frame("GC < 0.4"=as.numeric(rawData$seqAnno$gc) < 0.4,
                         "GC > 0.6"=as.numeric(rawData$seqAnno$gc) > 0.6,
                         check.names=FALSE)
  } else {
    gcTypes = NULL
  }
  if (!is.null(rawData$seqAnno$featWidth)){
    widthTypes = data.frame("width < 500nt"=as.numeric(rawData$seqAnno$featWidth) < 500, 
                            "width > 5000nt"=as.numeric(rawData$seqAnno$featWidth) > 5000,
                            check.names=FALSE)
  } else {
    widthTypes = NULL
  }
  
  if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
    pngNames = c()
    pngNames["def"] = "allPairs-scatter.png"
    ezAllPairScatter(signalCond, file=pngNames["def"], isPresent=isPresentCond, types=types)
    if (!is.null(gcTypes)){
      pngNames["gc"] = "allPairs-scatter-byGc.png"
      ezAllPairScatter(signalCond, file=pngNames["gc"], main="color by GC", isPresent=isPresentCond, types=gcTypes)      
    }
    if (!is.null(widthTypes)){
      pngNames["width"] = "allPairs-scatter-byWidth.png"
      ezAllPairScatter(signalCond, file=pngNames["width"], main="color by width", isPresent=isPresentCond, types=widthTypes)
    }
    writeImageRowToHtml(pngNames, con=html)
    flush(html)
  }

  for (i in 1:min(4, ncol(design))){
    for (cond in unique(design[,i])){
      idx = which(cond == design[,i])
      if (length(idx) > 1){
        idx = idx[order(samples[idx])] ## order alphabetically
        condName = paste(colnames(design)[i], cond)
        ezWrite("<h3>", condName, "</h3>", con=html)
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        ezScatter(NULL, signal[ ,idx], file=pngName, isPresent=isPresent[ ,idx], types=types,
                     lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
        ezWrite("<img src='", pngName, "'><br>", con=html)
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          ezScatter(NULL, signal[ ,idx], file=pngName, isPresent=isPresent[ ,idx], types=gcTypes,
                       lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          ezWrite("<img src='", pngName, "'><br>", con=html)          
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          ezScatter(NULL, signal[ ,idx], file=pngName, isPresent=isPresent[ ,idx], types=widthTypes,
                       lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          ezWrite("<img src='", pngName, "'><br>", con=html)          
        }
      }
    }
  }
  flush(html)
}


##' @title Writes test scatter plots
##' @description Writes test scatter plots to an html file.
##' @param html a connection to an html file.
##' @param param a list of parameters .....
##' @param x a dataset .....
##' @template result-template
##' @param seqAnno the sequence annotation. This is used if types is NULL.
##' @template types-template
##' @template roxygen-template
##' @return Returns
writeTestScatterPlots = function(html, param, x, result, seqAnno, colorRange=c(-3, 3),
   colors=getBlueRedScale(), types=NULL){
  if (is.null(types)){
      types = data.frame(row.names=rownames(x))
			if ("IsControl" %in% colnames(seqAnno)){
				types$Controls = seqAnno[ rownames(x), "IsControl"]
			}
  }
	msg = "Highlighting significants with: "
  if (!is.null(param$pValueHighlightThresh)){
    significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    types$Significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
		msg = paste(msg, "p <= ", param$pValueHighlightThresh)
		if (!is.null(param$log2RatioHighlightThresh)){
			msg = paste(msg, "and log ratio >= ", param$log2RatioHighlightThresh)
			if (!is.null(result$log2Ratio)){
				types$Significants = types$Significants & abs(result$log2Ratio) >= param$log2RatioHighlightThresh
			} else {
				types$Significants = types$Significants & result$log2Effect >= param$log2RatioHighlightThresh
			}
		}
  }
  ezWrite("<h2>Scatter Plots</h2>", con=html)
  ezWrite("<p>", msg, "</p>", con=html)
  ezWrite("<h3>Between-group Comparison</h3>", con=html)
  pngNames = character()
  pngNames["scatter"] = paste0(param$comparison, "-scatter.png")
  if (ncol(result$groupMeans) == 2 & !is.null(param$sampleGroup) & !is.null(param$refGroup)){
    sampleValues = 2^result$groupMeans[ , param$sampleGroup]
    refValues = 2^result$groupMeans[ , param$refGroup]
    ezScatter(refValues, sampleValues, file=pngNames["scatter"],
          isPresent=result$usedInTest,
          types=types,
          xlab=param$refGroup, ylab=param$sampleGroup)
		pngNames["volcano"] = paste0(param$comparison, "-volcano.png")
    ezVolcano(result$log2Ratio, result$pValue, file=pngNames["volcano"],
		   isPresent=result$usedInTest, types=types, main=param$comparison)
		pngNames["fdr-volcano"] = paste0(param$comparison, "-FDR-volcano.png")
    ezVolcano(result$log2Ratio, result$fdr, file=pngNames["fdr-volcano"],
		   isPresent=result$usedInTest, types=types, main=param$comparison, yType="FDR")
  } else {
    pngNames["allpair"] = paste0(param$comparison, "-scatter.png")
    ezAllPairScatter(2^result$groupMeans, file=pngNames["allpair"],
       isPresent=result$usedInTest, types=types)
  }
  myBreaks = seq(0, 1, by=0.002)
  pngNames["pValueHist"] = paste0(param$comparison, "-pValueHist.png")
  png(pngNames["pValueHist"], height=400, width=800)
  histUsed = hist(result$pValue[result$usedInTest], breaks=myBreaks, plot=FALSE)
  histAbs = hist(result$pValue[!result$usedInTest], breaks=myBreaks, plot=FALSE)
  xx =rbind(used=histUsed$counts, absent=histAbs$counts)
  xx = shrinkToRange(xx, c(0, max(xx["used", ])))
  #colnames(x) = histUsed$mids
  barplot(xx, space=0, border=NA, col=c("blue", "darkorange"), 
          xlab="p-value", ylab="counts", ylim=c(0, max(xx["used", ])),
          main="p-value histogram")
  abline(h=sum(result$usedInTest)/ncol(xx))
  at = c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
  axis(1, at=at*ncol(xx), labels = at)
  legend("top", c("used", "not expressed"), col=c("blue", "darkorange"), pch=20, cex=1)
  dev.off()
  writeImageRowToHtml(pngNames, con=html)
  
  if (is.null(x)){
    flush(html)
    return()
  }
  
  theRange = 2^(range(x, na.rm=TRUE))
  if (!ezIsSpecified(param$batch)){ ## TODO: we no longer use pairing, we now use batch which is more general; however these plots only work if batch is a real pairing
    for (group in unique(c(param$refGroup, colnames(result$groupMeans)))){
      idx = which(group == param$grouping)
      if (length(idx) > 1){
         #ezWrite("<h3>Intra-group Comparison: ", group, "</h3>", con=html)
         pngName = paste0(group, "-scatter.png")
         xlab = paste("Avg of", group)
         refValue = result$groupMeans[ , group]
         ezScatter(2^refValue, 2^x[, idx, drop=FALSE],
                file=pngName,
                isPresent=result$isPresent[, idx, drop=FALSE],
                types=types,
                lim=theRange, xlab=xlab)
         #ezWrite("<img src='", pngName, "'>", con=html)
         #ezWrite("<br>", con=html)
         if (ncol(result$groupMeans) == 2){
           otherGroup = setdiff(colnames(result$groupMeans), group)
           pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
           xlab = paste("Avg of", otherGroup)
           refValue = result$groupMeans[ , otherGroup]
           ezScatter(2^refValue, 2^x[, idx, drop=FALSE],
					        file=pngName,
                  isPresent=result$isPresent[, idx, drop=FALSE],
                  types=types,
                  lim=theRange, xlab=xlab)
           #ezWrite("<img src='", pngName, "'>", con=html)
           #ezWrite("<br>", con=html)
         }
      }
    }
  } else {
    
    #ezWrite("<h3>Pairs: ", param$sampleGroup, " over ", param$refGroup, "</h3>", con=html)
    use = param$grouping %in% c(param$sampleGroup, param$refGroup)
    if (all(table(param$batch[use], param$grouping[use]) == 1)){
      groups = paste(param$grouping, param$batch, sep="--")
      sampleGroups = sort(unique(groups[param$grouping == param$sampleGroup]))
      refGroups = sort(unique(groups[param$grouping == param$refGroup]))
      avgValues = averageColumns(x[ ,use], groups[use], mean)
      avgPresent= averageColumns(x[ ,use], groups[use], function(x){mean(x) > 0.5})
      sampleValues = avgValues[ , sampleGroups, drop=FALSE]
      refValues = avgValues[ , refGroups, drop=FALSE]
      samplePresent = avgPresent[ ,sampleGroups, drop=FALSE]
      refPresent = avgPresent[ , refGroups, drop=FALSE]
      pngName = paste0(param$sampleGroup, "-over-", param$refGroup, "-pairs.png")
      ezScatter(2^refValues, 2^sampleValues, file=pngName,
                isPresent=samplePresent | refPresent,
                types=types, lim=theRange, xlab=colnames(refValues))
      #ezWrite("<img src='", pngName, "'>", con=html)
      #ezWrite("<br>", con=html)
    }
  }
  flush(html)
}


##' @title Gets significant ratios
##' @description Gets significant ratios and counts them in a table.
##' @template result-template
##' @param pThresh a numerical indicating the p-value threshold.
##' @param ratioThresh a ratio threshold .....
##' @param genes if not NULL, the function counts the number of different genes that are significant.
##' @template roxygen-template
##' @return Returns a matrix containing the counts of significant ratios.
##' @examples
##' 1
getSignificantRatioCountsTable = function(result, pThresh=1/10^(1:5), ratioThresh = c(1, 1.5, 2, 3, 4, 8, 10), genes=NULL){
  
  ## counts the significant entries
  ## if genes is given counts the number of different genes that are significant
  if (!is.null(result$log2Ratio)){
    ratio = 2^result$log2Ratio
  } else {
    stop("not supported")
  }
  
  sigFcTable = ezMatrix(NA, rows=paste("p <", pThresh),
                        cols=c(paste("ratio >=", ratioThresh), paste("ratio <=", signif(1/ratioThresh, digits=3))))
  for (i in 1:length(pThresh)){
    for (j in 1:length(ratioThresh)){
      isSig = result$pValue < pThresh[i] & result$usedInTest == 1 & ratio >= ratioThresh[j]
      if (is.null(genes)){
        sigFcTable[i, j] = sum(isSig, na.rm=TRUE)
      } else {
        sigFcTable[i, j] = length(unique(na.omit(genes[isSig])))
      }
    }
    for (j in (length(ratioThresh)+1):(2*length(ratioThresh))){
      isSig = result$pValue < pThresh[i] & result$usedInTest == 1 & ratio <= 1/ratioThresh[j-length(ratioThresh)]
      if (is.null(genes)){
        sigFcTable[i, j] = sum(isSig, na.rm=TRUE)
      } else {
        sigFcTable[i, j] = length(unique(na.omit(genes[isSig])))
      }
    }
  }
  sigFcTable
}
