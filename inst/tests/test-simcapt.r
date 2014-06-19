context("Testing sim.capt()")

test_that("simulation produces correct output", {
    set.seed(8173)
    test.capt <- sim.capt(traps = example$traps, mask = example$mask,
                             infotypes = c("bearing", "dist", "toa"),
                             detfn = "ss",
                             pars = list(D = 2500, b0.ss = 90, b1.ss = 4,
                                 sigma.ss = 10, kappa = 70, alpha = 4,
                                 sigma.toa = 0.002), cutoff = 60)
    ## Is a list.
    expect_that(test.capt, is_a("list"))
    ## Correct components.
    expect_that(names(test.capt), equals(c("bincapt", "ss", "bearing",
                                           "dist", "toa")))
    ## Checking dimensions.
    dims <- sapply(test.capt, dim)
    expect_that(all(dims[1, ] == nrow(example$capt$bincapt)), is_true())
    expect_that(all(dims[2, ] == nrow(example$traps)), is_true())
    ## Checking equality with example$capt.
    expect_that(test.capt, equals(example$capt))
})

test_that("simulation errors", {
    ## Specifying incorrect SS link.
    expect_that(sim.capt(traps = example$traps, mask = example$mask, detfn = "ss",
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10), ss.link = "identity.link",
                         cutoff= 0),
                throws_error("The argument 'ss.link' must be either \"identity\" or \"log\""))
    ## Trying to simulate SS using 'infotypes'.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         infotypes = c("ss", "toa"),
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10, sigma.toa = 0.002), cutoff= 0),
                throws_error("Signal strength information is simulated by setting argument 'detfn' to \"ss\"."))
    ## Cutoff provided when not required.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         cutoff = 0),
                gives_warning("The argument 'cutoff' is being ignored, as 'detfn' is not \"ss\"."))
    ## Signal strength link provided when not required.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         ss.link = "identity"),
                gives_warning("The argument 'ss.link' is being ignored, as 'detfn' is not \"ss\"."))
    ## Extra parameter.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5,
                             sigma.toa = 0.002)),
                throws_error("The following must be named components of the list 'pars': "))
    ## Missing parameter.
    expect_that(sim.capt(traps = traps, mask = mask,
                         pars = list(D = 2000, g0 = 0.75)),
                throws_error("The following must be named components of the list 'pars': "))
})
