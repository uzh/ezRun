###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Format an integer into a character using only digits
##' @description convenience function that prevents scientific notation and surrounding white space. Such characters are need for system calls to commandline tools.
##' @param x an integer.
##' @return Returns the input integer as a character.
##' @template roxygen-template
##' @examples 
##' ezIntString(4.5e7)
ezIntString = function(x){
  stopifnot(as.integer(x) == x) ## must be an integer number; can be in numeric (floating point) format.
  format(x, scientific=FALSE, trim=TRUE)
}

##' @title Converts a number into a character representing millions
##' @description used for labels etc. where typically millions of reads are indicated.
##' @param x an integer.
##' @return Returns a character representing the input integer in millions.
##' @template roxygen-template
##' @examples 
##' mioString(4e7)
mioString = function(x){
  format(x/1e6, scientific=FALSE, digits=3, trim=TRUE)
}

##' @title Archives files
##' @description Archives one or several files to a .zip extension. Optionally, a separate name can be given to the archive.
##' @param zipName a character naming the working directory to be archived.
##' @param inputs one or several input files to archive.
##' @param zipped optional character to name the archive.
##' @return Returns the file name of the zipped archive.
##' @template roxygen-template
##' @examples 
##' write("a file","exampleFile")
##' zipFile("exampleFile")
zipFile = function(inputs, zipped=NULL){
  
  if (is.null(zipped)){
    x = paste(sub("\\.[[:alpha:]]*$", "", inputs), collapse="_")
    zipped = paste0(x, ".zip")
  }
  inputs = paste0("'", inputs, "'") ## surround with quotes to support files that contain white space
  inputs = paste(inputs, collapse=" ") ## concat the files if multiple files are given;
  cmd = paste("zip -q -r", zipped, inputs, sep=" ")
  ezSystem(cmd, echo=FALSE)
  return(zipped)
}

##' @describeIn zipFile Archives the whole working directory.
zipWorkingDir = function(zipName){
  
  zipFile(".", zipped=zipName)
  file.remove(list.files(pattern=".*.RData"))
  ezSystem(paste("zip -q -r", zipName, ".", sep=" "))
  return(zipName)
}

##' @title Modified version of read.table
##' @description Modified version of \code{read.table()} with a different default.
##' @param file the name of the file to read the data from.
##' @return Returns a data.frame.
##' @template roxygen-template
##' @template addargs-template
##' @templateVar fun read.table()
##' @seealso \code{\link[utils]{read.table}}
##' @examples 
##' m1 = ezMatrix(1:20, rows=1:5, cols=1:4)
##' ezWrite.table(m1, "exampleTable")
##' ezRead.table("exampleTable")
## simple wrapper to read.table with useful defaults
ezRead.table = function(file, header=TRUE, sep="\t", as.is=TRUE, row.names=1, quote="", skip=0, comment.char="", check.names=FALSE, ...){
  data = read.table(file, header=header, sep=sep, as.is=as.is, row.names=row.names, quote=quote, skip=skip,
                    comment.char=comment.char, check.names=check.names, ...)
  return(data)
}

##' @title Modified version of write.table
##' @description Modified version of \code{write.table()} with a different default. Best used together with \code{ezRead.table()}.
##' @param values a vector, matrix of data.frame to write a table from.
##' @param file the name of the output file.
##' @param head the names of the header.
##' @param digits the number of digits to round to, if rounding is desired.
##' @return Returns a table written into a separate file.
##' @template roxygen-template
##' @seealso \code{\link[utils]{write.table}}
##' @seealso \code{\link{ezRead.table}}
##' @examples 
##' m1 = matrix(seq(0.01, 1, 0.01), 10)
##' colnames(m1) = letters[1:10]
##' ezWrite.table(m1, "exampleTable", digits=1)
## the example produces an interesting behaviour: signif(0.15,1) == 0.2 == signif(0.15,1). R seems to prefer to round to even digits...
## convenience wrapper that writes data.frames and matrices with row and column names such that they can be read again; also supports vectors
## NOTE: if x has row and column names and you use write.table(x, file="foo.txt", row.names=TRUE, col.names=TRUE) you can not read in that table with read.table("foo.txt")
ezWrite.table = function(values, file=file, head="Identifier", row.names=TRUE, col.names=TRUE,
                         append=FALSE, quote=FALSE, sep="\t", na="NA", digits=NA){
  
  if (is.vector(values)){
    values = as.data.frame(values)
  }
  if (!is.na(digits)){
    if (is.data.frame(values)){
      for (i in 1:length(values)){
        if (is.numeric(values[[i]])){
          hasDecimals = suppressWarnings({as.integer(values[[i]]) != values[[i]]})
          hasDecimals[is.na(hasDecimals)] = FALSE
          if (any(hasDecimals)){
            values[[i]] = signif(values[[i]], digits=digits)
          }
        }
      }
    } else {
      if (is.numeric(values)){
        values = signif(values, digits=digits)
      }
    }
  }
  
  if (row.names){
    if (col.names){
      write.table(matrix(c(head, colnames(values)), nrow=1), file=file, sep=sep, quote=quote, col.names=FALSE, row.names=FALSE, append=append, na=na)
      write.table(values, file=file, sep=sep, quote=quote, col.names=FALSE, row.names=TRUE, append=TRUE, na=na)
    } else {
      write.table(values, file=file, sep=sep, quote=quote, col.names=FALSE, row.names=TRUE, append=append, na=na)
    }
  } else {
    write.table(values, file=file, sep=sep, quote=quote, col.names=col.names, row.names=FALSE, append=append, na=na)
  }
}

##' @title Saves an interactive table
##' @description Saves an interactive table accessible with the provided \code{tableLink}.
##' @param values a data.frame or table to create an interactive table from.
##' @param tableLink a character ending with .html representing the link to the interactive table
##' @param digits the number of digits to round to, if rounding is desired.
##' @param colNames a character vector specifying the column names of the interactive table.
##' @param title a character representing the title of the interactive table.
##' @param format formatting options passed as an expression. The table argument in formatting functions must be named \code{interactiveTable}.
##' @param envir the environment to evaluate \code{format} in.
##' @template roxygen-template
##' @seealso \code{\link[DT]{datatable}}
##' @seealso \code{\link[DT]{saveWidget}}
##' @examples 
##' tableLink = "exampleTable.html"
##' table = data.frame(a=c(1.11, 2:100), b=201:300)
##' ezInteractiveTable(table, tableLink)
ezInteractiveTable = function(values, tableLink, digits=NULL, colNames=colnames(values), title="", format=NULL, envir=parent.frame()){
  require(DT, quietly=TRUE)
  require(htmltools, quietly=TRUE)
  if (!is.null(digits)){
    for (i in 1:ncol(values)) {
      if(typeof(values[ ,i]) == "double"){
        values[ ,i] = signif(values[ ,i], digits=digits)
      }
    }
    captionText = paste("Numeric values are rounded to", digits, "digits.")
    caption = htmltools::tags$caption(htmltools::h1(title), 
                                      htmltools::p(captionText))
  } else {
    caption = htmltools::tags$caption(htmltools::h1(title))
  }
  interactiveTable = datatable(values, 
                                   extensions=c("Buttons"), filter="top", caption=caption, colnames=colNames,
                                   options=list(dom = 'Bfrtip', buttons = c('colvis','copy', 'csv', 'excel', 'pdf', 'print'), pageLength=25, autoWidth=TRUE)
                                   )
  if (!is.null(format)){
    currEnv = environment()
    interactiveTable = eval(format, envir=c(envir, currEnv))
  }
  saveWidget(interactiveTable, tableLink)
}

##' @title Generates an interactive table
##' @description Generates an interactive table embeded in rmarkdown document.
##' @param values a data.frame or table to create an interactive table from.
##' @param digits the number of digits to round to, if rounding is desired.
##' @param colNames a character vector specifying the column names of the interactive table.
##' @param title a character representing the title of the interactive table.
##' @param format formatting options passed as an expression. The table argument in formatting functions must be named \code{interactiveTable}.
##' @param envir the environment to evaluate \code{format} in.
##' @template roxygen-template
##' @seealso \code{\link[DT]{datatable}}
##' @seealso \code{\link[DT]{saveWidget}}
##' @examples 
##' table = data.frame(a=c(1.11, 2:100), b=201:300)
##' ezInteractiveTableRmd(table)
ezInteractiveTableRmd = function(values, digits=NULL, 
                                 colNames=colnames(values), title="", 
                                 format=NULL, envir=parent.frame()){
  require(DT, quietly=TRUE)
  require(htmltools, quietly=TRUE)
  if (!is.null(digits)){
    for (i in 1:ncol(values)) {
      if(typeof(values[ ,i]) == "double"){
        values[ ,i] = signif(values[ ,i], digits=digits)
      }
    }
    captionText = paste("Numeric values are rounded to", digits, "digits.")
    caption = htmltools::tags$caption(htmltools::h1(title), 
                                      htmltools::p(captionText))
  } else {
    caption = htmltools::tags$caption(htmltools::h1(title))
  }
  interactiveTable <- datatable(values, 
                                extensions=c("Buttons"), filter="top", 
                                caption=caption, colnames=colNames,
                                options=list(dom = 'Bfrtip', 
                                             buttons = c('colvis','copy', 'csv',
                                                         'excel', 'pdf', 'print'), 
                                             pageLength=25, autoWidth=TRUE)
                                )
  if (!is.null(format)){
    currEnv = environment()
    interactiveTable = eval(format, envir=c(envir, currEnv))
  }
  return(interactiveTable)
}

##' @title Write in a single line
##' @description Concatenates its arguments and writes it as a single line.
##' @param ... the arguments to concatenate.
##' @param sep a character specifying how to separate the arguments from each other.
##' @param collapse a character specifying how to separate entries of each argument from each other.
##' @return Returns a single line written into a separate file.
##' @template roxygen-template
##' @template connection-template
##' @examples 
##' con = file("example.txt", "w")
##' ezWrite("concatenates","these",4,"arguments",con=con,sep=" ")
##' close(con)
## concatenates its arguments and writes it as a single line
## usage:
## con = file("result.txt", "w")
## ezWrite(a, b, c, con=con)
ezWrite = function(..., sep="", collapse=" ", con=stdout()){
  args = list(...)  ## see function message
  #args = sapply(args, print)# as.character)
  text = paste(sapply(args, paste, collapse=collapse), collapse=sep)
  writeLines(text, con=con)
}

##' @title Enforces a valid file name
##' @description Replaces invalid characters and returns a file name with these replaced
##' @param name a character representing a file name.
##' @param replace a character defining the replacement character(s).
##' @return Returns a character where invalid characters for file names are replaced.
##' @template roxygen-template
##' @examples
##' ezValidFilename("example:filename")
##' ezValidFilename("or?to/remove(them",replace="")
## replace characters that are not allowed in filenames
ezValidFilename = function(name, replace="_"){
  gsub("[$%#!?/:;() '=]", replace, name)
}

##' @describeIn getSuffix Gets and removes the suffix of a file name.
removeSuffix = function(filename){
  sapply(filename, function(x){sub("\\.$", "", sub(getSuffix(x), "", x))})
}

##' @title Gets the suffix of a file name
##' @description Gets the suffix of a file name or the file name without the suffix.
##' @param filename a filename to extract the suffix from.
##' @return Returns a character with either the suffix of the file name or the file name without suffix.
##' @template roxygen-template
##' @examples
##' getSuffix("example.file")
##' removeSuffix("example.file")
getSuffix = function(filename){
  return( sub(".*\\.", "", filename))
}

##' @title Is x an absolute file path?
##' @description Checks whether \code{x} is an absolute file path.
##' @param x a character representing a file path to check.
##' @return Returns TRUE if \code{x} starts with a /.
##' @template roxygen-template
##' @examples
##' ezIsAbsolutePath("/absolutepath")
ezIsAbsolutePath = function(x){
  !is.null(x) & grepl("^/", x)
}


ezRandomString = function(length){
  paste(sample(c(0:9, letters, LETTERS),
               length, replace=TRUE), collapse="")
}
