context("Testing conversion functions")

## Note convert.mask() and convert.capt() are tested in test-fits.r.

test_that("Converting/creating capture histories", {
    dets <- lightfooti$dets
    traps <- lightfooti$traps
    capt.test <- convert.pamguard(dets = dets, mics = traps)
    expect_that(capt.test, is_identical_to(lightfooti$capt))
    mask.test <- create.mask(traps = traps, buffer = 40)
    expect_that(mask.test, is_identical_to(lightfooti$mask))
})
