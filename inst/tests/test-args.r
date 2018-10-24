context("Testing fit.ascr() arguments")

test_that("fixing parameters", {
    simple.capt <- example$capt["bincapt"]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 0.9))
    ## Checking that g0 is not estimated.
    expect_that(get.par(fit, "g0"), is_equivalent_to(0.9))
    ## Checkint that phase is set to -1.
    expect_that(fit$phases$g0, is_equivalent_to(-1))
    ## Checking that all parameters have a phase.
    all.pars <- c(fit$D.betapars, fit$detpars, fit$suppars)
    phase.pars <- names(fit$phases)
    expect_that(sort(phase.pars), equals(sort(all.pars)))
    ## Checking that no other parameters are set to -1.
    active.phases <- c(fit$phases[phase.pars != "g0"],
                       recursive = TRUE)
    expect_that(all(active.phases > -1), is_true())
})

test_that("start values", {
    simple.capt <- example$capt["bincapt"]
    ## Fit original model.
    fit.start <- fit.ascr(capt = simple.capt, traps = example$traps,
                          mask = example$mask)
    ## Provide a single start value.
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, sv = list(D = 2145))
    ## Check that estimates are the same.
    relative.error <- max(abs((coef(fit.start) - coef(fit))/
                              coef(fit)))
    expect_that(relative.error < 1e-4, is_true())
    ## Check start value is passed correctly.
    expect_that(fit$args$sv$`D.(Intercept)`, is_equivalent_to(log(2145)))
})

test_that("parameter bounds", {
    simple.capt <- example$capt["bincapt"]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, bounds = list(D = c(0, 5000)))
    ## Check that bounds object is a list.
    expect_that(is.list(fit$args$bounds), is_true())
    ## Check that bounds for D set appropriately.
    expect_that(fit$args$bounds$`D.(Intercept)`, equals(c(log(1e-20), log(5000))))
    ## Check that bounds for g0 still set to defaults.
    expect_that(fit$args$bounds$g0, equals(c(0, 1)))
})

test_that("local integration", {
    simple.capt <- example$capt[c("bincapt", "toa")]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, local = TRUE)
    expect_that(fit$args$local, is_true())
})
