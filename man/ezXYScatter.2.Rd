\name{ezXYScatter.2}
\alias{ezXYScatter.2}
\title{
  XYScatter plot
}
\description{
  Plot the interactive XYScatter plot or XYScatter plot with labelled genes.
}
\usage{
  ezXYScatter.2(xVec, yVec, absentColor = "gray", shrink = FALSE, 
                xlim = range(xVec, yVec, na.rm = TRUE), 
                ylim = xlim, isPresent = NULL, names = NULL, types = NULL, 
                colors = rainbow(ncol(types)), main = NULL, 
                xlab = NULL, ylab = NULL, labelGenes = NULL,
                mode = c("plotly", "ggplot2"))
}
\arguments{
  \item{xVec}{
    \code{numeric}(n): the group means of one group.
  }
  \item{yVec}{
    \code{numeric}(n): the group means of another group.
  }
  \item{absentColor}{
    \code{character}(1): the colour for absent genes.
  }
  \item{shrink}{
    \code{boolean}(1): shrink the group means or not.
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
     \code{character}(1): the title.
  }
  \item{xlab, ylab}{
    \code{character}(1): the xlab and ylab.
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
  \code{\link{ezVolcano}}
}
\examples{
  sampleValues <- c(0.28, 0.87, 31.70, 361.18)
  refValues <- c(0.27, 0.60, 28.81, 364.13)
  isPresent <- c(TRUE, TRUE, FALSE, TRUE)
  types <- data.frame(Significants=c(TRUE, FALSE, FALSE, TRUE),
	                    ZincFingers=c(TRUE, FALSE, FALSE, FALSE))
  names <- c("SALL1", "PAX6", "IRX3", "MEIS1")
  xlab <- "A"
  ylab <- "B"
  main="Comparison of average expression"
  
  p_plotly <- ezXYScatter.2(xVec=refValues, yVec=sampleValues,
                            isPresent=isPresent, 
                            types=types, names=names,
                            xlab=xlab, ylab=ylab,
                            main="Comparison of average expression",
                            labelGenes=c("SALL1", "IRX3"),
                            mode="plotly")

  p_ggplot2 <- ezXYScatter.2(xVec=refValues, yVec=sampleValues,
                             isPresent=isPresent, 
                             types=types, names=names,
                             xlab=xlab, ylab=ylab,
                             main="Comparison of average expression",
                             labelGenes=c("SALL1", "IRX3"),
                             mode="ggplot2")
}
