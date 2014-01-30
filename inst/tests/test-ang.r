context("Testing ang")

test_that("bearing fitting", {
    ## Fitting model.
    ang.capt <- example.capt[c("bincapt", "ang")]
    fit <- admbsecr(capt = ang.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2394.109251, 5.214332270, 67.65495606)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(248.67, 0.17466, 9.5809)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("kappa"))
})
