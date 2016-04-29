context("Testing sim.capt()")

test_that("simulation produces correct output", {
    set.seed(8173)
    test.capt <- sim.capt(traps = example$traps, mask = example$mask,
                          infotypes = c("bearing", "dist", "toa"),
                          detfn = "ss",
                          pars = list(D = 2500, b0.ss = 90, b1.ss = 4,
                              sigma.ss = 10, kappa = 70,
                              alpha = 4, sigma.toa = 0.002),
                          ss.opts = list(cutoff = 60, directional = FALSE), test.detfn = TRUE)
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
    test.ss <- c(61.5810096413274, 85.7765377386673, 71.9380707339211, 0,
                 77.7630468936444, 82.4358785366641)
    expect_that(test.capt$ss[3, ], equals(test.ss))
    test.bearings <- c(0.593511290743994, 1.23589000772559, 2.18054479341382, 0,
                       5.59633425691473, 3.53634915849748)
    expect_that(test.capt$bearing[3, ], equals(test.bearings))
    test.dist <- c(5.86854382347385, 5.61710285563385, 4.29201576315865, 0,
                   3.53331938014272, 2.48135390226467)
    expect_that(test.capt$dist[3, ], equals(test.dist))
    test.toa <- c(0.022316230815507, 0.0111928740201518, 0.0158444826486618, 0,
                  0.0068018316556033, 0.00766521624617204)
    expect_that(test.capt$toa[3, ], equals(test.toa))
    ## Probably also best to compare to pre-specified data. Previous
    ## simulations used an older simulation function so this one is
    ## not validated. I think it's fine though.
    test.capt <- sim.capt(traps = example$traps, mask = example$mask,
                          detfn = "hn",
                          pars = list(D = 2500, g0 = 0.9, sigma = 5))
    capt <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0),
                   nrow = 3)
    expect_that(test.capt$bincapt[1:3, ], equals(capt))
    ## Testing simulation with heterogeneity in signal strengths.
    ## Same capture history analysed in test-fits.r.
    set.seed(6434)
    pars <- list(D = 500, b0.ss = 90, b1.ss = 4, sigma.b0.ss = 5, sigma.ss = 10)
    test.capt <- sim.capt(traps = example$traps, mask = example$mask, detfn = "ss",
                     pars = pars, ss.opts = list(cutoff = 60))
    test.ss <- c(75.8665334318406, 62.4175659997338, 0, 61.7592272206934, 0, 0)
    expect_that(test.capt$ss[3, ], equals(test.ss))
    ## Testing simulation for first-call data.
    ## Same capture history analysied in test-fits.r.
    set.seed(8298)
    traps <- cbind(c(0, 0, 21, 21, 42, 42), c(0, 42, 0, 42, 0, 42))
    mask <- create.mask(traps, buffer = 1250, spacing = 25)
    pars <- list(D = 5, b0.ss = 60, b1.ss = 0.1, sigma.ss = 5)
    lower.cutoff <- 52.5
    cutoff <- 55
    test.capt <- sim.capt(traps = traps, mask = mask, detfn = "ss",
                     infotypes = NULL, pars = pars,
                     ss.opts = list(cutoff = lower.cutoff,
                         ss.link = "identity"),
                     cue.rates = Inf, first.only = TRUE)
    test.ss <- c(0, 58.6358978925394, 53.0558496057121, 53.9825023770212, 0, 
                 65.2750802998676)
    expect_that(test.capt$ss[1103, ], equals(test.ss))
})

## Note simulation tests for mulitple-cue models is in test-fit.r

test_that("simulation errors", {
    ## Specifying incorrect SS link.
    expect_that(sim.capt(traps = example$traps, mask = example$mask, detfn = "ss",
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10),
                         ss.opts = list(cutoff = 0, ss.link = "identity.link")),
                throws_error("The argument 'ss.link' must be either \"identity\" or \"log\""))
    ## Trying to simulate SS using 'infotypes'.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         infotypes = c("ss", "toa"),
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10, sigma.toa = 0.002),
                         ss.opts = list(cutoff = 0)),
                throws_error("Signal strength information is simulated by setting argument 'detfn' to \"ss\"."))
    ## Cutoff provided when not required.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         ss.opts = list(cutoff = 0)),
                gives_warning("The argument 'ss.opts' is being ignored, as 'detfn' is not \"ss\"."))
    ## Signal strength link provided when not required.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         ss.opts = list(ss.link = "identity")),
                gives_warning("The argument 'ss.opts' is being ignored, as 'detfn' is not \"ss\"."))
    ## Extra parameter.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5,
                             sigma.toa = 0.002)),
                throws_error("The following must be named components of the list 'pars': "))
    ## Missing parameter.
    expect_that(sim.capt(traps = example$traps, mask = example$mask,
                         pars = list(D = 2000, g0 = 0.75)),
                throws_error("The following must be named components of the list 'pars': "))
})
