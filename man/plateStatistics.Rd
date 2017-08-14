\name{plateStatistics}
\alias{plateStatistics}
\alias{heatmapPlate}

\title{
  Bias inspection per plate
}
\description{
  Calculate counts/concentration in a plate layout and plot in heatmaps.
  
}
\usage{
  plateStatistics(dataset, colname = c("Read Count",
                                       "LibConc_qPCR [Characteristic]",
                                       "LibConc_100_800bp [Characteristic]"))
  heatmapPlate(x, title="unnamed", center=TRUE, log10=TRUE)                                       
}

\arguments{
  \item{dataset}{
    A \code{data.frame} read from \dQuote{input_dataset.tsv}.
  }
  \item{title}{
    A \code{character(1)}: title of the plate heatmap.
  }
  \item{x}{
    A \code{matrix} of counts/concentration of a plate.
  }
  \item{center}{
    A \code{boolean}. When TRUE, plot the heatmaps value from half of median to
    twice of median.
  }
  \item{log10}{
    A \code{boolean}. When TRUE, plot the \code{log10} value.
  }
}
\details{
  \code{plateStatistics} fetch the plate position information from the column
  \dQuote{PlatePosition [Characteristic]}.
  The current supported format of \dQuote{PlatePosition [Characteristic]} is
  \dQuote{PlateNumber_[A-Z][Integer]}, e.g. \dQuote{1_A2}.
  
  \code{heatmapPlate} plot the heatmap with the results 
  from \code{plateStatistics}.
}
\value{
  \code{plateStatistics}: a \code{list} of \code{list} of \code{matrix} with 
  rows in alphabets and columns in numbers of the plate.
  
  \code{heatmapPlate}: \code{plotly} object of heatmap.
}
\references{
  Plate layout: \url{https://www.dnaseq.co.uk/images/samples-on-well-plate.png}
}
\author{
  Ge Tan
}

\examples{
  \dontrun{
    dataset <- ezRead.table("/srv/gstore/projects/p2438/Fastqc_18564_2017-06-07--14-06-33/input_dataset.tsv")
    ans <- plateStatistics(dataset)
    
    p1 <- heatmapPlate(ans[[1]][["Read Count"]], center=TRUE, log10=TRUE)
    p2 <- heatmapPlate(ans[[1]][["Read Count"]], center=TRUE, log10=FALSE)
    p3 <- heatmapPlate(ans[[1]][["Read Count"]], center=FALSE, log10=TRUE)
    p4 <- heatmapPlate(ans[[1]][["Read Count"]], center=FALSE, log10=FALSE)
    subplot(p1, p3,p2,p4, nrows=2)
  }
}
