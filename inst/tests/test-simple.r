context("Testing simple")

test_that("simple simulation and fitting works correctly", {
    ## Fitting model.
    fit <- admbsecr(capt = simple.capt, traps = simple.traps,
                    mask = simple.mask)                
    ## Checking parameter values.
    pars.test <- c(3094.6992410, 0.9073728, 4.2213375)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(525.48, 0.10928, 0.46901)
    relative.error <- max(abs((stdEr(fit)[1:3] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
})
