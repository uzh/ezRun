## STAR (like most aligners) gives multimapping reads MAPQ 0-3; with the SUSHI defaults
## minMapQuality 10 and keepMultiHits true, featureCounts drops them before keepMultiHits applies.
test_that("FeatureCounts warns when minMapQuality removes the multimapping reads keepMultiHits would count", {
  expect_warning(warn_multi_hits_dropped(list(keepMultiHits = TRUE, minMapQuality = 10)), "minMapQuality 10")
  expect_warning(warn_multi_hits_dropped(list(keepMultiHits = TRUE, minMapQuality = 1)), "minMapQuality 1")
  expect_no_warning(warn_multi_hits_dropped(list(keepMultiHits = TRUE, minMapQuality = 0)))
  expect_no_warning(warn_multi_hits_dropped(list(keepMultiHits = FALSE, minMapQuality = 10)))
})
