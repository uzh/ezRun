\name{ezVolcano}
\alias{ezVolcano}
\title{
  Volcano plot
}
\description{
  Plot the interactive volvano plot or volcano plot with labelled genes.
}
\usage{
  ezVolcano(log2Ratio, pValue, yType = c("p-value", "FDR"), 
              xlim = NULL, ylim = NULL, isPresent = NULL, 
              names = NULL, types = NULL, 
              colors = rainbow(ncol(types)), main = NULL, labelGenes = NULL,
              mode=c("plotly", "ggplot2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{log2Ratio}{
    \code{numeric}(n): the log2 fold change of each gene
  }
  \item{pValue}{
    \code{numeric}(n): the p-value or FDR from 
    differential gene expression analysis.
  }
  \item{yType}{
    \code{character}(1): p-value or FDR in \code{pValue}.
  }
  \item{xlim}{
    \code{numeric}(2): the range of xaxis.
  }
  \item{ylim}{
    \code{numeric}(2): the range of yaxis.
  }
  \item{isPresent}{
    \code{boolean}(n): the genes are expressed or not.
  }
  \item{names}{
    \code{character}(n): the gene names.
  }
  \item{types}{
    \code{data.frame} of \code{boolean}: indices for each type per column.
  }
  \item{colors}{
    colors for the types.
  }
  \item{main}{
    \code{character}(1): the title of volcano plot.
  }
  \item{labelGenes}{
    \code{character}(n): subset of \code{names} to label.
  }
  \item{mode}{
    generate interactive plot with \sQuote{plotly} or 
    static plot with \sQuote{ggplot2}.
    To label the genes, only \sQuote{ggplot2} is supported.
  }
}
\value{
  A \code{plotly} or \code{ggplot2} object.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{ezXYScatter}}
}
\examples{
  log2Ratio <- -3:3
  pValue <- 10^(c(-3:0, -1:-3))
  yType <- "p-value"
  isPresent <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
  types <- data.frame(Significants=c(TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE),
                      ZincFinger=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  names <- c("SALL1", "PAX6", "MEIS1", "IRX3", "ARX", "TBX2", "SHH")
  main <- "A--over--B"
  labelGenes <- c("MEIS1", "ARX")
  p_plotly <- ezVolcano(log2Ratio=log2Ratio, pValue=pValue,
                        yType="p-value", isPresent=isPresent,
                        names=names,
                        types=types, main=main, labelGenes=labelGenes,
                        mode="plotly")
  export(p_plotly, file="volcano_plotly.pdf")
  library(htmlwidgets)
  saveWidget(as_widget(p_volcano), "volcano_plotly.html")
  
  p_ggplot2 <- ezVolcano(log2Ratio=log2Ratio, pValue=pValue,
                         yType="p-value", isPresent=isPresent,
                         names=names,
                         types=types, main=main, labelGenes=labelGenes,
                         mode="ggplot2")
  p_ggplot2 + theme(title=element_text(size=20),
                    axis.title.x=element_text(size=30),
                    axis.title.y=element_text(size=30))
}
