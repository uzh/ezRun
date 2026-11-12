## Every app must build: an appDefaults entry with an extra column (EzAppSCEVANApp had
## `Values = c("human", "mouse")`) makes rbind() fail in initialize(), so every job of
## the app stops before any work and its citation() cannot be called for the Methods.

test_that("every EzApp can be instantiated", {
  ns <- asNamespace("ezRun")
  apps <- Filter(function(n) grepl("^EzApp", n) && inherits(get(n, ns), "refObjectGenerator"),
                 ls(ns, all.names = TRUE))
  expect_gt(length(apps), 100)
  failed <- Filter(Negate(is.null), lapply(stats::setNames(apps, apps), function(n) {
    tryCatch({ get(n, ns)$new(); NULL }, error = function(e) conditionMessage(e))
  }))
  expect(length(failed) == 0,
         paste0("apps that fail to build:\n", paste(names(failed), unlist(failed), sep = ": ", collapse = "\n")))
})
