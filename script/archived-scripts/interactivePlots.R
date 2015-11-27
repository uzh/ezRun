


# volcano scatter plots: red dots mouse over description
# use addPlot and add.plot.interactivity from ReporteRs.
library(ReporteRs)

myfoo = function(){
  data(mtcars)
  plot(mtcars$mpg, mtcars$wt, type="n", xlab="mpg", ylab="wt")
  dbl_click_actions = paste0("window.open('https://www.google.fr/#q=", row.names( mtcars ),"');")
  popup.labels = paste(row.names( mtcars ), "<br>", "double click to google the car")
  add.plot.interactivity(fun=points, x=mtcars$mpg, y=mtcars$wt,
                          col=mtcars$gear, pch=16,
                          dblclick.actions=dbl_click_actions,
                          popup.labels=popup.labels)
}


## own testing
doc = openBsdocReport("example")
doc = addPlot(doc, fun = myfoo, fontname = "serif")
closeBsdocReport(doc, "example.html")

refValues = 1:100
names(refValues) = letters[1:20]
sampleValues = runif(100)
myPlotFoo = function(){
  ezScatter(x=refValues, y=sampleValues)
  add.plot.interactivity(fun=points, x=refValues, y=sampleValues, pch=16, popup.labels = names(refValues))
}
doc = addPlot(doc, fun=myPlotFoo, fontname="serif")


