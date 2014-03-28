context("Testing model fits")

test_that("simple fitting -- half normal", {
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

test_that("simple fitting -- hazard rate", {
    ## Fitting model.
    simple.capt <- example.capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask, detfn = "hr")
    ## Checking parameter values.
    pars.test <- c(2601.486636, 0.9283243400, 6.782970400, 6.409336100)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(348.14, 0.074121, 0.71436, 1.4893)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma", "z")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking estimate equivalence with secr.fit.
    library(secr)
    mask.secr <- convert.mask(example.mask)
    capt.secr <- convert.capt(example.capt, example.traps)
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
    bearing.capt <- example.capt[c("bincapt", "bearing")]
    fit <- admbsecr(capt = bearing.capt, traps = example.traps,
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
    ## Checking get.par() with SS stuff.
    expect_that(get.par(fit, "all"), equals(coef(fit, c("fitted", "derived"))))
    expect_that(get.par(fit, "b0.ss"), equals(fit$coefficients["b0.ss"]))
})

test_that("toa fitting", {
    ## Fitting model.
    toa.capt <- example.capt[c("bincapt", "toa")]
    fit <- admbsecr(capt = toa.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2238.943787, 5.434543258, 0.00184227498)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(285.68, 0.30612, 0.00019294)
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
    joint.capt <- example.capt[c("bincapt", "ss", "toa")]
    fit <- admbsecr(capt = joint.capt, traps = example.traps,
                    mask = example.mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    cutoff = 60)
    ## Checking parameter values.
    pars.test <- c(2518.778360, 91.11204602, 4.312022549,
                   10.85171521, 0.001954810264)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(266.32, 1.5753, 0.24763, 0.56685, 0.00020987)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
})

test_that("joint bearing/dist fitting", {
    joint.capt <- example.capt[c("bincapt", "bearing", "dist")]
    fit <- admbsecr(capt = joint.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2476.527851, 5.105681597, 70.35386170, 3.941970659)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(249.77, 0.15269, 9.4000, 0.34958)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("g0", "sigma")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals(c("kappa", "alpha")))
})
      
test_that("multiple calls fitting", {
    simple.capt <- example.capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example.traps,
                    mask = example.mask, fix = list(g0 = 1),
                    call.freqs = c(9, 10, 11))
    pars.test <- c(2358.737165, 5.262653889, 10, 0.05384236, 235.87375890)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit, c("fitted", "derived")) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    expect_that(confint(fit),
                throws_error("Standard errors not calculated; use boot.admbsecr()"))
    expect_that(summary(fit), is_a("list"))
})
