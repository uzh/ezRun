

## all plots would be generated from a plotter class

EzPlotter =
  setRefClass("EzPlotter",
              fields=c("name", "data", "helpText", "mouseOverText"),
              methods=list(
                plot=function(){},
                plotPng=function(){
                  ## create a png file and call then plot
                },
                plotPdf=function(){
                  ## create a pdf file call then plot
                },
                writeData=function(){}
              )
  )
              
EzPlotterIris =
  setRefClass("EzPlotterIris",
              contains="EzPlotter",
              methods=list(
                plot=function(){
                  plot(data$Sepal.length, Sepal.Width)
                }
              )
  )


## in the report generating scripts I want to write

theDoc = bsdoc( title = 'My document' )

irisData = iris 

theDoc = addEzImage(theDoc, ezPlotIris$new(data=irisData, param=NULL,
                                           name="myIrisPlot"),
                    mouseOverText="myMouseOverText", helpText="myHelpText") ## if mouseOverText and helpText is null the default help text will be used


addEzImage = function(theDoc, ezPlotter, mouseOverText=NULL, helpText=NULL, addPdfLink=TRUE) {
  ## create the png plot and add it to the doc add the mouse over text
  # put it into a 2x1 table
  # put in the second row the pdf link and the help text
}

#plot(iris$Sepal.Length, iris$Sepal.Width)
