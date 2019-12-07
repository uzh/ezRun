\name{writeErrorReport}
\alias{writeErrorReport}
\title{Writes an error report}
\usage{
writeErrorReport(htmlFile, param = param, error = "Unknown Error")
}
\arguments{
\item{htmlFile}{a character, ending with .html, representing the file path of an html file.}

\item{param}{a list of parameters to extract the \code{name} from.}

\item{error}{a character vector representing the error message(s).}
}
\description{
Writes an error report to an html file. Also creates the file and closes it.
}
\examples{
param = ezParam()
htmlFile = "example.html"
writeErrorReport(htmlFile, param)
}
\seealso{
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
