

doc = openBsdocReport("testSmoScaAndAddPlot")
bla = data.frame(a=runif(100), b=runif(100))
bla = data.frame(a=runif(10000), b=runif(10000))

myplotFoo = function(){
  print(system.time({
    smoothScatter(bla$a, bla$b, nbin=32)
  }))
  print(system.time({
    add.plot.interactivity(points, x=rep(0.5,10), y=seq(0.1, 0.9, length.out = 10))
  }))
}
myplotFoo2 = function(){
  print(system.time({
    ezSmoothScatter(bla$a, bla$b, nbin=32)
  }))
  print(system.time({
    add.plot.interactivity(points, x=rep(-5,10), y=seq(-1, -9, length.out = 10))
  }))
}
myplotFoo3 = function(){
  print(system.time({
    ezScatter(bla$a, bla$b)
  }))
  print(system.time({
    add.plot.interactivity(points, x=rep(0.01,10), y=seq(0.01, 0.5, length.out = 10))
  }))
}

print(system.time({
  addPlot(doc, myplotFoo, font="")
}))
print(system.time({
  addPlot(doc, myplotFoo2, font="")
}))
print(system.time({
  addPlot(doc, myplotFoo3, font="")
}))

closeBsdocReport(doc, "test.html")
