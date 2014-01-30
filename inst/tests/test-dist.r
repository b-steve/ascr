context("Testing dist")

test_that("dist fitting", {
    ## Fitting model.
    dist.capt <- example.capt[c("bincapt", "dist")]
    fit <- admbsecr(capt = dist.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2477.446007, 5.104500978, 4.053234339)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(260.88, 0.18057, 0.45189)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("alpha"))
})
