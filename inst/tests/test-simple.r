context("Testing simple")

test_that("simple fitting", {
    ## Fitting model.
    simple.capt <- example.capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask)
    ## Checking parameter values.
    pars.test <- c(2358.737165, 1.000000000, 5.262653889)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(345.89, 0.38032)
    relative.error <- max(abs((stdEr(fit)[-2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking estimate equivalence with secr.fit.
    library(secr)
    mask.secr <- convert.mask(example.mask)
    capt.secr <- convert.capt(example.capt, example.traps)
    options(warn = -1)
    fit.secr <- secr.fit(capthist = capt.secr, mask = mask.secr, trace = FALSE)
    options(warn = 0)
    coefs.secr <- numeric(n.pars)
    invlog <- function(x) exp(x)
    for (i in 1:n.pars){
        coefs.secr[i] <- eval(call(paste("inv", fit.secr$link[[i]], sep = ""),
                                   fit.secr$fit$par[i]))
    }
    relative.error <- max(abs((coef(fit) - coefs.secr)/coefs.secr))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("fixing parameters", {
    simple.capt <- example.capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 0.9))
    ## Checking that g0 is not estimated.
    expect_that(getpar(fit, "g0"), is_equivalent_to(0.9))
    ## Checkint that phase is set to -1.
    expect_that(fit$phases$g0, is_equivalent_to(-1))
    ## Checking that all parameters have a phase.
    all.pars <- c("D", fit$detpars, fit$suppars)
    phase.pars <- names(fit$phases)
    expect_that(sort(phase.pars), equals(sort(all.pars)))
    ## Checking that no other parameters are set to -1.
    active.phases <- c(fit$phases[phase.pars != "g0"],
                       recursive = TRUE)
    expect_that(all(active.phases > -1), is_true())

})

test_that("start values", {
    simple.capt <- example.capt["bincapt"]
    ## Fit original model.
    fit.start <- admbsecr(capt = simple.capt, traps = example.traps,
                          mask = example.mask)
    ## Provide a single start value.
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask, sv = list(D = 2145))
    ## Check that estimates are the same.
    relative.error <- max(abs((coef(fit.start) - coef(fit))/
                              coef(fit)))
    expect_that(relative.error < 1e-4, is_true())
    ## Check start value is passed correctly.
    expect_that(fit$sv$D, is_equivalent_to(2145))
})

test_that("parameter bounds", {
    simple.capt <- example.capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask, bounds = list(D = c(0, 5000)))
    ## Check that bounds object is a list.
    expect_that(is.list(fit$bounds), is_true())
    ## Check that bounds for D set appropriately.
    expect_that(fit$bounds$D, equals(c(0, 5000)))
    ## Check that bounds for g0 still set to defaults.
    expect_that(fit$bounds$g0, equals(c(0, 1)))
})
