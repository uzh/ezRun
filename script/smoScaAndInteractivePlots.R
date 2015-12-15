

doc = openBsdocReport("testSmoScaAndAddPlot")
myplotFoo = function(){
  bla = data.frame(a=runif(1000), b=runif(1000))
  smoothScatter(bla$a, bla$b)
  add.plot.interactivity(points, x=rep(0.5,10), y=seq(0.1, 0.9, length.out = 10))
}
addPlot(doc, myplotFoo, font="")
closeBsdocReport(doc, "test.html")
