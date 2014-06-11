context("Testing admbsecr() errors")

test_that("error for missing bincapt", {
    test.capt <- example$capt["toa"]
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask),
                throws_error("The binary capture history must be provided as a component of 'capt'."))
})

test_that("error for non-matrix in capt object", {
    test.capt <- example$capt["bincapt"]
    test.capt$bearing <- 1:10
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask),
                throws_error("At least one component of 'capt' is not a matrix."))
})

test_that("error for mismatch in number of individuals or traps", {
    ## Testing error checking for equality in number of rows.
    test.capt <- example$capt[c("bincapt", "bearing", "dist")]
    test.capt$bearing <- test.capt$bearing[-1, ]
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask),
                throws_error("Components of 'capt' object have different dimensions."))
    ## Testing error checking for equality in number of columns.
    test.capt <- example$capt[c("bincapt", "bearing", "dist")]
    test.capt$bearing <- test.capt$bincapt[, -1]
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask),
                throws_error("Components of 'capt' object have different dimensions."))
    ## Testing error checking for matching in number of traps.
    test.capt <- example$capt["bincapt"]
    test.capt$bincapt <- test.capt$bincapt[, -1]
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask),
                throws_error("There must be a trap location for each column in the components of 'capt'."))
})

test_that("arguments are of the right type", {
    test.capt <- example$capt["bincapt"]
    ## Testing error checking for 'sv' type.
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask, sv = c(D = 1000, g0 = 1)),
                throws_error("The 'sv' argument must be 'NULL' or a list."))
    ## Testing error checking for 'bounds' type.
    bounds <- matrix(1:6, ncol = 2)
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask, bounds = bounds),
                throws_error("The 'bounds' argument must be 'NULL' or a list."))
    bounds <- list(D = 1000, g0 = c(0, 0.9))
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask, bounds = bounds),
                throws_error("Each component of 'bounds' must be a vector of length 2."))
    ## Testing error checking for 'fix' type.
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask, fix = c(D = 1000, g0 = 1)),
                throws_error("The 'fix' argument must be 'NULL' or a list."))
})

test_that("ss-related parameters set up correctly", {
    test.capt <- example$capt[c("bincapt", "ss")]
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask,
                         sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                         cutoff = 60, ss.link = "identity.link"),
                throws_error("ss.link must be either \"identity\" or \"log\""))
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask,
                         sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                         cutoff = 60, detfn = "hr"),
                gives_warning("Argument 'detfn' is being ignored as signal strength information is provided in 'capt'. A signal strength detection function has been fitted instead."))
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask,
                         sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10)),
                throws_error("Argument 'cutoff' is missing."))

})

test_that("exe.type argument is correct", {
    test.capt <- example$capt["bincapt"]
    ## Testing error checking for 'sv' type.
    expect_that(admbsecr(capt = test.capt, traps = example$traps,
                         mask = example$mask, exe.type = "diff"),
                throws_error("Argument 'exe.type' must be \"old\" or \"new\"."))
})

test_that("Extra components of 'sv' and 'fix' are removed", {
    test.capt <- example$capt["bincapt"]
    ## Testing warning checking 'sv' components.
    expect_that(test.fit <- admbsecr(capt = test.capt, traps = example$traps,
                                     mask = example$mask, sv = list(z = 5, g0 = 1),
                                     fix = list(g0 = 1)),
                gives_warning("Some parameters listed in 'sv' are not being used. These are being removed."))
    expect_that(test.fit, equals(example$fits$simple.hn))
    ## Testing warning checking 'fix' components.
    expect_that(test.fit <- admbsecr(capt = test.capt, traps = example$traps,
                                     mask = example$mask, fix = list(g0 = 1, foo = 0)),
                gives_warning("Some parameters listed in 'fix' are not being used. These are being removed."))
    expect_that(test.fit, equals(example$fits$simple.hn))
})
