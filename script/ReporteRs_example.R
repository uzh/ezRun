require( ggplot2 )
doc = docx( title = 'My document' )

doc = addTitle( doc , 'First 5 lines of iris', level = 1)
doc = addFlexTable( doc , light.table(iris[1:5, ]), layout.properties= get.light.tableProperties())

doc = addTitle( doc , 'ggplot2 example', level = 1)
myggplot = qplot(Sepal.Length, Petal.Length, data = iris, color = Species, size = Petal.Width )
doc = addPlot( doc = doc , fun = print, x = myggplot )

doc = addTitle( doc , 'Text example', level = 1)
doc = addParagraph( doc, 'My tailor is rich.', stylename = 'Normal' )

writeDoc( doc, 'my_first_doc.docx' )
