## The writer's prompt (methods_identity / methods_task in R/01classes.R): the rules adopted
## from methods-task-no-stage-prescription, and the rules of the audited prompt it dropped.

flatPrompt <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))

test_that("the prompt has the plain-text, no-path, one-step and job-setting rules", {
  p <- flatPrompt(c(EzApp$new()$methods_identity(), EzApp$new()$methods_task()))
  for (rule in c(
    "the authoritative record of this analysis, but often not a complete one",
    "Plain text only: no Markdown and no code formatting, including backticks.",
    "Describe settings in words, not as command-line options: write 'a minimum read length of 25 bases', not '--length_required 25'.",
    "Never include a path in any form",
    "no folder or file name that shows where something is stored",
    "describe a supplied file by what it contains",
    "Describe the input as samples, not files.",
    "one step of a larger analysis and is combined with the text of the other steps, so do not introduce the analysis as a whole",
    "Leave out settings that only affect how the job ran: compute resources, threads, memory, output and temporary locations, and logging",
    "whether the run completed",
    "Where the record does not show a setting, write [not recorded]; do not describe it as a default unless app_facts.txt states it.",
    "Some things are not worth describing: the framework that launched the job (SUSHI, ezRun, and R itself unless the version of an R package matters), the layout of the report and the names of output files."))
    expect_true(grepl(rule, p, fixed = TRUE), label = rule)
  expect_false(grepl("complete and authoritative", p, fixed = TRUE))
})

test_that("the prompt keeps the sessionInfo, vendor-pipeline and completeness rules", {
  p <- flatPrompt(c(EzApp$new()$methods_identity(), EzApp$new()$methods_task()))
  for (rule in c(
    "A module that was loaded but never invoked is not part of the analysis.",
    "a package that appears only in a sessionInfo listing",
    "A vendor pipeline often runs steps such as clustering or differential expression inside itself",
    "Unless the record shows the step was skipped, say the step was performed and write [not recorded] for the settings",
    "every tool that was actually invoked must appear somewhere in the text",
    "Before finishing, check your draft against the tools named in the scripts and logs.",
    "Never assert that something was NOT done.",
    "Do not complete a partial identifier."))
    expect_true(grepl(rule, p, fixed = TRUE), label = rule)
})

test_that("[not recorded] is for one specific setting, and vendor steps follow run_summary.txt", {
  p <- flatPrompt(EzApp$new()$methods_task())
  for (rule in c(
    "[not recorded] names one specific setting a reader would need; never write it for 'any other' or 'additional' settings",
    "never for a setting listed in the parameters: an empty value there means the app's own rule applied",
    "Where run_summary.txt lists the steps a vendor pipeline ran and did not run, describe only the ones it ran."))
    expect_true(grepl(rule, p, fixed = TRUE), label = rule)
})
