

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

#plot(iris$Sepal.Length, iris$Sepal.Width)
