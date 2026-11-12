## Cell Ontology names contain commas ("CD4-positive, alpha-beta T cell"), so a list with
## one name per line must be split on newlines only; a single line is split on commas.

test_that("one name per line keeps names that contain commas", {
  x <- "CD4-positive, alpha-beta T cell\nCD8-positive, alpha-beta T cell\nB cell\n"
  expect_equal(parseReferenceCellTypes(x),
               c("CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "B cell"))
})

test_that("Windows line endings and blank lines are dropped", {
  expect_equal(parseReferenceCellTypes("B cell\r\n\r\n natural killer cell \r\n"),
               c("B cell", "natural killer cell"))
})

test_that("a single line is split on commas", {
  expect_equal(parseReferenceCellTypes("B cells, NK cells,plasma cells"),
               c("B cells", "NK cells", "plasma cells"))
})

test_that("empty input gives no cell types", {
  expect_length(parseReferenceCellTypes(""), 0)
  expect_length(parseReferenceCellTypes(NULL), 0)
})

test_that("the app default keeps the CD4 and CD8 T cell names whole", {
  default <- EzAppSCEVANApp$new()$appDefaults["referenceCellTypes", "DefaultValue"]
  expect_true("CD4-positive, alpha-beta T cell" %in% parseReferenceCellTypes(default))
  expect_true("CD8-positive, alpha-beta T cell" %in% parseReferenceCellTypes(default))
  expect_length(parseReferenceCellTypes(default), 5)
})
