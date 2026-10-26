context("getSpecies: refBuild -> species")

## getSpecies gates the whole CellMarker2/AUCell annotation section: anything it
## calls "other" gets no annotation at all, silently. It was a bare
## startsWith(refBuild, "Homo_sapiens") / "Mus_musculus", which three refBuild
## spellings in live gStore results do not satisfy even though they ARE human or
## mouse. Measured over 1,407 ScSeurat-family runs in /srv/gstore (2024-09 .. 2026-09):
##
##   13 mouse runs (p40923, p40924, Jan-Mar 2026) passed refBuild as an ABSOLUTE
##      path, "/srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"
##    3 runs passed a bare assembly name, "GRCm38"
##
## Those 16 runs lost mouse annotation to string parsing, not to biology.

test_that("canonical genus-first refBuilds still resolve", {
  expect_equal(getSpecies("Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-03-01"), "Human")
  expect_equal(getSpecies("Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"), "Mouse")
})

test_that("an absolute reference path resolves (regression: 13 lost mouse runs)", {
  expect_equal(
    getSpecies("/srv/GT/reference/Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"),
    "Mouse")
  expect_equal(
    getSpecies("/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-03-01"),
    "Human")
})

test_that("a trailing slash or doubled separator does not break the match", {
  expect_equal(getSpecies("/srv/GT/reference//Mus_musculus/GENCODE/GRCm39"), "Mouse")
  expect_equal(getSpecies("Homo_sapiens/"), "Human")
})

test_that("a bare assembly name resolves (regression: 3 lost runs)", {
  expect_equal(getSpecies("GRCm38"), "Mouse")
  expect_equal(getSpecies("GRCm39"), "Mouse")
  expect_equal(getSpecies("GRCh38.p13"), "Human")
  expect_equal(getSpecies("GRCh37"), "Human")
})

test_that("genuinely unsupported species are still 'other' - the gate must keep working", {
  ## This is the negative control. A fix that makes everything resolve to Human
  ## or Mouse would pass every test above and be worse than the bug.
  expect_equal(getSpecies("Danio_rerio/GENCODE/GRCz11/Annotation/Release_110-2023-08-01"), "other")
  expect_equal(getSpecies("Canis_familiaris/Ensembl/ROS_Cfam_1.0/Annotation/Release_110"), "other")
  expect_equal(getSpecies("Felis_catus/Ensembl/Felis_catus_9.0"), "other")
  expect_equal(getSpecies("Equus_caballus/Ensembl/EquCab3.0"), "other")
  expect_equal(getSpecies("/srv/GT/reference/Danio_rerio/GENCODE/GRCz11"), "other")
})

test_that("a human-mouse chimera is 'other', deliberately", {
  ## Neither species' marker sets are valid on a mixed-genome matrix, and the
  ## gene symbols carry a species suffix, so scoring it as Human would be wrong.
  ## 10 such runs exist; they SHOULD skip annotation.
  expect_equal(getSpecies("Chimera_GRCh38.p13_GRCm39"), "other")
  expect_equal(getSpecies("Chimera_GRCm39_mRatBN7.2"), "other")
  expect_equal(getSpecies("Chimera_human_mouse"), "other")
})

test_that("empty and missing refBuilds do not error", {
  expect_equal(getSpecies(""), "other")
  expect_equal(getSpecies(NA_character_), "other")
  expect_equal(getSpecies(NULL), "other")
  expect_equal(getSpecies(character(0)), "other")
})

test_that("a species whose name merely CONTAINS a supported genus is not misread", {
  ## Guard against a lazy grepl() fix: Mus_minutoides has a reference installed
  ## at FGCZ and is not Mus musculus.
  expect_equal(getSpecies("Mus_minutoides/Ensembl/MMinutoides1.0"), "other")
})
