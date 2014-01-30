context("Testing admbsecr() errors")

test_that("error for missing bincapt", {
    test.capt <- example.capt["toa"]
    expect_that(admbsecr(capt = test.capt, traps = example.traps,
                         mask = example.mask),
                throws_error("The binary capture history must be provided as a component of 'capt'."))
})
