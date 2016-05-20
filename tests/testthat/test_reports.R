context("Tests the functions in reports.R")

cwd = getwd()
setwdNew("./ReportTest")




doc = openBsdocReport("Testreport")
dat = data.frame(a=1:5, b=11:15)
param = ezParam()
input = EzDataset$new(file=system.file("extdata/yeast_10k/dataset.tsv", package = "ezRun", mustWork = TRUE))

test_that("Tests ezFlexTable() and ezGrid()", {
  flexTable = ezFlexTable(dat)
  flexGrid = ezGrid(dat)
  expect_is(flexTable, "FlexTable")
  expect_is(flexGrid, "FlexTable")
  addParagraph(doc, "Testing ezFlexTable()")
  addFlexTable(doc, flexTable)
  addParagraph(doc, "Testing ezGrid()")
  addFlexTable(doc, flexGrid)
})


test_that("tests adding links",{
  addParagraph(doc, ezPot("normal link", hyperlink = "http://www.google.com"))
  addParagraph(doc, ezPot("link with explicit type", hyperlink = "https://github.com/jgm/pandoc/archive/1.15.1.tar.gz", linkType="application/image"))
  addParagraph(doc, ezPot("link opening in a new window", hyperlink = "http://www.google.com", linkTarget = "_blank"))
})



test_that("Tests ezImageFileLink()", {
  x = 1:10
  plotCmd = expression({
    plot(x)
    text(2,1, "my Text")
  })
  fileLink = ezImageFileLink(plotCmd)
  expect_is(fileLink, "character")
  expect_true(grepl("span><img src", fileLink))
  addParagraph(doc, "Testing ezImageFileLink()")
  addParagraph(doc, fileLink)
})

test_that("Tests addDataset()", {
  addParagraph(doc, "Testing addDataset()")
  addDataset(doc, input$meta, param)
})

test_that("Tests addTxtLinksToReport()", {
  addParagraph(doc, "Testing addTxtLinksToReport()")
  addTxtLinksToReport(doc, "dataset.tsv")
})

# test_that("Tests newWindowLink()", {
#   newWindow = newWindowLink("testLink")
#   expect_is(newWindow, "pot")
#   expect_true(grepl("javascript", newWindow[[1]]$value))
# })

test_that("Tests openBsdocReport() and closeBsdocReport()", {
  expect_is(doc, "bsdoc")
  closeBsdocReport(doc, "TestReport.html")
})

setwd(cwd)
ezSystem("rm -fr ./ReportTest")
