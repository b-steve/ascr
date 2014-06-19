context("Testing model fits")

test_that("simple fitting -- half normal", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example$traps,
                    mask = example$mask)
    ## Checking parameter values.
    pars.test <- c(2429.62882766283, 0.999999985691245, 5.36419062611109)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(370.82, 0.41458)
    relative.error <- max(abs((stdEr(fit)[-2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking get.par().
    expect_that(get.par(fit, "D"), equals(fit$coefficients["D"]))
    expect_that(get.par(fit, c("D", "sigma")),
                equals(fit$coefficients[c("D", "sigma")]))
    expect_that(get.par(fit, "esa"), equals(fit$coefficients["esa"]))
    expect_that(get.par(fit, c("D", "esa")),
                equals(fit$coefficients[c("D", "esa")]))
    expect_that(get.par(fit, "all"), equals(coef(fit, c("fitted", "derived"))))
    expect_that(get.par(fit, "fitted"), equals(coef(fit, "fitted")))
    ## Testing some generic functions.
    expect_that(is.list(summary(fit)), is_true())
    expect_that(confint(fit), is_a("matrix"))
    ## Checking hess argument.
    fit.nohess <- admbsecr(capt = simple.capt, traps = example$traps,
                           mask = example$mask, hess = FALSE)
    expect_that(coef(fit.nohess), equals(coef(fit)))
    ## Checking estimate equivalence with secr.fit.
    library(secr)
    mask.secr <- convert.mask(example$mask)
    capt.secr <- convert.capt(example$capt["bincapt"], example$traps)
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

test_that("simple fitting -- hazard rate", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, detfn = "hr")
    ## Checking parameter values.
    pars.test <- c(2647.32709086757, 0.862838540084305, 7.29740531390158, 
                   6.98650302516657)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(410.24, 0.076434, 0.68136, 2.247)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma", "z")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking estimate equivalence with secr.fit.
    library(secr)
    mask.secr <- convert.mask(example$mask)
    capt.secr <- convert.capt(example$capt["bincapt"], example$traps)
    options(warn = -1)
    fit.secr <- secr.fit(capthist = capt.secr, mask = mask.secr, detectfn = 1,
                         trace = FALSE)
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

test_that("bearing fitting", {
    ## Fitting model.
    bearing.capt <- example$capt[c("bincapt", "bearing")]
    fit <- admbsecr(capt = bearing.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2326.9511045605, 5.50755363089523, 58.2671714153125)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(245.34, 0.20411, 7.7149)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("kappa"))
})

test_that("dist fitting", {
    ## Fitting model.
    dist.capt <- example$capt[c("bincapt", "dist")]
    fit <- admbsecr(capt = dist.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2477.20317036246, 5.30067437855294, 6.23928789994028)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(249.42, 0.17043, 0.80207)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("alpha"))
    ## Testing get.par() with fixed values.
    expect_that(get.par(fit, "D"), equals(fit$coefficients["D"]))
    expect_that(get.par(fit, "g0"),
                is_equivalent_to(fit$args$sv[["g0"]]))
    expect_that(get.par(fit, c("g0", "sigma")),
                        is_equivalent_to(c(1, fit$coefficients["sigma"])))
    expect_that(get.par(fit, c("sigma", "g0")),
                is_equivalent_to(c(fit$coefficients["sigma"], 1)))
})

test_that("ss fitting", {
    ## Fitting model.
    ss.capt <- example$capt[c("bincapt", "ss")]
    fit <- admbsecr(capt = ss.capt, traps = example$traps,
                    mask = example$mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    cutoff = 60)
    ## Checking parameter values.
    pars.test <- c(2440.99968246751, 88.2993844240013, 3.7663100027822, 
                   10.8142267688236)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(311.63, 1.8566, 0.27901, 0.67581)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking get.par() with SS stuff.
    expect_that(get.par(fit, "all"), equals(coef(fit, c("fitted", "derived"))))
    expect_that(get.par(fit, "b0.ss"), equals(fit$coefficients["b0.ss"]))
})

test_that("toa fitting", {
    ## Fitting model.
    toa.capt <- example$capt[c("bincapt", "toa")]
    fit <- admbsecr(capt = toa.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2348.61001583513, 5.47656668927085, 0.00209625029558564)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(305.38, 0.32532, 0.00020158)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
    ## Checking get.par() with supplementary info parameters.
    expect_that(get.par(fit, "sigma.toa"), equals(fit$coefficients["sigma.toa"]))
    expect_that(get.par(fit, c("esa", "sigma.toa")),
                equals(fit$coefficients[c("esa", "sigma.toa")]))
})

test_that("joint ss/toa fitting", {
    joint.capt <- example$capt[c("bincapt", "ss", "toa")]
    fit <- admbsecr(capt = joint.capt, traps = example$traps,
                    mask = example$mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    cutoff = 60)
    ## Checking parameter values.
    pars.test <- c(2439.05009421856, 89.6628138603152, 3.85848602239484, 
                   10.2843370178208, 0.00212227986424938)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(260.1, 1.635, 0.23517, 0.55907, 0.00019848)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
})

test_that("joint bearing/dist fitting", {
    joint.capt <- example$capt[c("bincapt", "bearing", "dist")]
    fit <- admbsecr(capt = joint.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2378.76835878674, 5.43409577354988, 57.0835463656355, 
                   5.17899019377086)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(235.59, 0.16267, 7.4596, 0.49989)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals(c("kappa", "alpha")))
})

test_that("multiple calls fitting", {
    ## Checking fit works.
    simple.capt <- example$capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1),
                    call.freqs = c(9, 10, 11))
    pars.test <- c(2429.62646548165, 5.36419300496923, 10, 0.0555640833502627, 
                   242.962646548165)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit, c("fitted", "derived")) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    expect_that(confint(fit),
                throws_error("Standard errors not calculated; use boot.admbsecr()"))
    expect_that(summary(fit), is_a("list"))
    ## Checking hess argument.
    fit.hess <- admbsecr(capt = simple.capt, traps = example$traps,
                         mask = example$mask, fix = list(g0 = 1),
                         call.freqs = c(9, 10, 11), hess = TRUE)
    expect_that(coef(fit.hess), equals(coef(fit)))
    expect_that(is.na(stdEr(fit.hess, "all")["Da"]), is_true())
    ses.test <- c(370.82, 0.41458)
    relative.error <- max(abs((stdEr(fit.hess)[1:2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})
