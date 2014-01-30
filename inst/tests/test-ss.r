context("Testing ss")

test_that("ss fitting", {
    ## Fitting model.
    ss.capt <- example.capt[c("bincapt", "ss")]
    fit <- admbsecr(capt = ss.capt, traps = example.traps,
                    mask = example.mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    cutoff = 60)
    ## Checking parameter values.
    pars.test <- c(2631.276509, 90.94073454, 4.396885247, 10.552865862)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(306.31, 1.6572, 0.27738, 0.64280)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
})
