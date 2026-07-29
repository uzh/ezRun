context("CellRanger annotatable-reference alias")

## CellRanger >=10.1.0 runs its local Pan-Human Azimuth model by default, but
## only when reference.json declares a genome matching hg19/GRCh38/GRCh39 as a
## substring. FGCZ references declare the index directory name, so annotation is
## skipped SILENTLY (count succeeds, no outs/cell_types/). These tests pin the
## alias that works around it without touching the shared production references.
##
## Uses a synthetic reference directory, so no 31 GB index is needed.

makeFakeRef <- function(genome) {
  d <- file.path(tempfile("ref"))
  dir.create(file.path(d, "fasta"), recursive = TRUE)
  dir.create(file.path(d, "genes"), recursive = TRUE)
  dir.create(file.path(d, "star"), recursive = TRUE)
  writeLines("x", file.path(d, "star", "SA"))
  writeLines(
    jsonlite::toJSON(
      list(genomes = genome, fasta_hash = "abc", mem_gb = 60),
      auto_unbox = TRUE
    ),
    file.path(d, "reference.json")
  )
  d
}

FGCZ_GENOME <- "genes_10XGEX_SC_Mt_rRNA-Mt_tRNA-protein_coding-rRNA-tRNA_Index"
HUMAN_BUILD <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
MOUSE_BUILD <- "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"

test_that("an FGCZ human reference gets an alias declaring GRCh38", {
  ref <- makeFakeRef(FGCZ_GENOME)
  out <- ezRun:::cellRangerAnnotatableRef(
    ref,
    list(refBuild = HUMAN_BUILD),
    aliasParent = tempfile("alias")
  )

  expect_false(identical(out, ref))
  expect_equal(
    as.character(jsonlite::fromJSON(file.path(out, "reference.json"))$genomes),
    "GRCh38"
  )

  ## The heavy directories must be SYMLINKS, never copies: the real index is
  ## 31 GB and copying it per job would be ruinous.
  for (sub in c("fasta", "genes", "star")) {
    expect_true(nzchar(Sys.readlink(file.path(out, sub))))
  }
  ## and they must still resolve to the real content
  expect_true(file.exists(file.path(out, "star", "SA")))

  ## the original reference is untouched
  expect_equal(
    as.character(jsonlite::fromJSON(file.path(ref, "reference.json"))$genomes),
    FGCZ_GENOME
  )
})

test_that("non-human references are left alone", {
  ## The model is pan-HUMAN. Aliasing a mouse reference to GRCh38 would be
  ## actively harmful: it would coax CellRanger into annotating mouse cells
  ## against human symbols.
  ref <- makeFakeRef("genes_10XGEX_SC_protein_coding_Index")
  out <- ezRun:::cellRangerAnnotatableRef(
    ref,
    list(refBuild = MOUSE_BUILD),
    aliasParent = tempfile("alias")
  )
  expect_identical(out, ref)
})

test_that("a reference that is already annotatable is not aliased", {
  ## 10x-supplied references declare GRCh38 already; building an alias would be
  ## pointless churn.
  for (g in c("GRCh38", "hg19", "GRCh39", "GRCh38-2020-A")) {
    ref <- makeFakeRef(g)
    out <- ezRun:::cellRangerAnnotatableRef(
      ref,
      list(refBuild = HUMAN_BUILD),
      aliasParent = tempfile("alias")
    )
    expect_identical(out, ref)
  }
})

test_that("multi-genome references are left alone", {
  ## Cell annotation rejects multi-genome references outright, so there is
  ## nothing to gain and a barnyard reference must not be relabelled human.
  ref <- makeFakeRef(c("GRCh38_x", "GRCm39_y"))
  out <- ezRun:::cellRangerAnnotatableRef(
    ref,
    list(refBuild = HUMAN_BUILD),
    aliasParent = tempfile("alias")
  )
  expect_identical(out, ref)
})

test_that("a broken reference never fails the job", {
  ## An optional annotation must never cost us the count. Anything unexpected
  ## returns the original path so cellranger runs exactly as before.
  expect_identical(
    ezRun:::cellRangerAnnotatableRef(
      "/no/such/reference",
      list(refBuild = HUMAN_BUILD),
      aliasParent = tempfile("alias")
    ),
    "/no/such/reference"
  )

  noJson <- tempfile("refNoJson")
  dir.create(noJson)
  expect_identical(
    ezRun:::cellRangerAnnotatableRef(
      noJson,
      list(refBuild = HUMAN_BUILD),
      aliasParent = tempfile("alias")
    ),
    noJson
  )
})
