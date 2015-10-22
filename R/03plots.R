###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title The R5 class representing the basic plotter
##' @description This class is usually not used, but all other plotters inherit from it.
##' @field name a character specifying the name of the plotter instance.
##' @field data a data.frame, matrix, table or another sensible data object to create plots from.
##' @field helpText a character specifying the help text of the plotter.
##' @field mouseOverText a character specifying the mouseOver text of images created with the plotter.
##' @template roxygen-template
##' @examples
##' EzPlotter$new()
EzPlotter =
  setRefClass("EzPlotter",
              fields = c("name", "data", "helpText", "mouseOverText"),
              methods = list(
                plot = function(...)
                {
                  "Plots \\code{data} with the default plot function from the graphics package."
                  graphics::plot(data, ...)
                },
                plotPng = function(file=NULL, width=480, height=480, ...)
                {
                  "Creates a .png file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".png", sep="")
                  }
                  png(filename = filename, width, height)
                  .self$plot(...)
                  dev.off()
                  return(filename)
                },
                plotPdf = function(file=NULL, width=480, height=480, ...)
                {
                  "Creates a .pdf file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".pdf", sep="")
                  }
                  width = round(width/72, digits=2)
                  height = round(height/72, digits=2)
                  pdf(file = filename, width, height)
                  .self$plot(...)
                  dev.off()
                  return(filename)
                },
                writeData = function()
                {
                  
                }
              )
  )

# A simple iris plotter as an example:
EzPlotterIris =
  setRefClass("EzPlotterIris",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL)
                {
                  "Initializes the fields and uses the iris dataset from R."
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterIris"
                  }
                  data <<- iris
                  helpText <<- "Iris is a flower dataset."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot=function(...)
                {
                  "Plots the provided iris data with optional arguments passed with ..."
                  graphics::plot(data$Sepal.Length, data$Sepal.Width, ...)
                }
              )
  )
