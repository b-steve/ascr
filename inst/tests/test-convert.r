context("Testing conversion functions")

## Note convert.mask() and convert.capt.to.secr() are tested in test-fits.r.

test_that("Converting/creating capture histories", {
    dets <- lightfooti$dets
    traps <- lightfooti$traps
    ## Testing old allocation method.
    capt.test <- convert.pamguard(dets = dets, mics = traps, new.allocation = FALSE)
    ## New convert.pamguard() retains original time data.
    capt.test$toa[capt.test$toa > 0] <- capt.test$toa[capt.test$toa > 0] - (capt.test$toa[1, 2] - 1)
    expect_that(capt.test, is_identical_to(lightfooti$capt))
    ## Testing new allocation method.
    capt.test.new <- convert.pamguard(dets = dets, mics = traps)
    expect_that(capt.test.new, is_identical_to(lightfooti$capt.new))
    mask.test <- create.mask(traps = traps, buffer = 40)
    expect_that(mask.test, is_identical_to(lightfooti$mask))
})
