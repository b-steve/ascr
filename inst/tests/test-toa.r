context("Testing toa")

test_that("simple simulation and fitting works correctly", {
    ## Fitting model.
    toa.capt <- example.capt[c("bincapt", "toa")]
    fit <- admbsecr(capt = toa.capt, traps = example.traps,
                    fix = list(g0 = 1), mask = example.mask)                
    ## Checking parameter values.
    pars.test <- c(2238.943787, 5.434543258, 0.00184227498)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(285.68, 0.30612, 0.00018903)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
})
