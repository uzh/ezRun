# tables within a table within a table works flawlessly with this code, I never experienced the <br></br> problem
# perhaps tables in the functions goClusterTable() and goUpDownTables() should be done as a list with each flex table being passed as.html()
param = ezParam(userParam = list('refBuild' = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'))
clusterResult = list("GO"=list("A"=letters[1:3], "B"=letters[4:6], "C"=letters[7:9]), nClusters=3, clusterColors=rainbow(3))
names(clusterResult$GO) = c("CC", "BP", "MF")
doc = openBsdocReport("My title")
# goLink = as.html(ezGrid(c("Background color corresponds to the row colors in the heatmap plot.",
#                           as.html(goClusterTable(param, clusterResult)))))
#tbl = ezGrid(ezFrame("Cluster Plot"=imgLinks("Rplot001.png"), "GO categories of feature clusters"=goLink), header.columns = TRUE)
innerTables = list()
for (i in seq(2,24,2)){
  innerTables[[i]] = as.html(ezFlexTable(matrix(rep(i, i), ifelse(i %% 4 == 0, 4, 2)), header.columns = TRUE))
}
table1 = ezMatrix(unlist(innerTables), rows=1:3, cols=letters[1:4])
tableLink = as.html(ezFlexTable(table1, header.columns = TRUE))
file.remove("Rplot001.png")
png()
plot(10:1)
dev.off()
tbl = ezGrid(ezFrame("Cluster Plot"=imgLinks("Rplot001.png"), "tables"=tableLink), header.columns = TRUE)
doc = addFlexTable(doc, tbl)
closeBsdocReport(doc, "example.html")
