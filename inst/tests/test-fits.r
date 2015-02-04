context("Testing model fits")

test_that("simple fitting -- half normal", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example$traps,
                    mask = example$mask)
    args.simple <- list(capt = simple.capt, traps = example$traps,
                        mask = example$mask)
    ## Checking parameter values.
    pars.test <- c(2267.73950296093, 0.999999965720685, 5.39011185815489)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(351.86, 0.42008)
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
    suppressPackageStartupMessages(library(secr))
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
    ## Checking fitting with local integration.
    fit.loc <- admbsecr(capt = simple.capt, traps = example$traps,
                        mask = example$mask, local = TRUE)
    args.loc <- list(capt = simple.capt, traps = example$traps,
                     mask = example$mask, local = TRUE)
    pars.test <- c(2267.76041594695, 0.999999965778877, 5.390081422004)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit.loc) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ses.test <- c(351.84, 0.42003)
    relative.error <- max(abs((stdEr(fit.loc)[-2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking parallel function.
    fits <- par.admbsecr(n.cores = 2, args.simple, args.loc)
    expect_that(fits[[1]], is_identical_to(fit))
    expect_that(fits[[2]], is_identical_to(fit.loc))
    fits.list <- par.admbsecr(n.cores = 2,
                              arg.list = list(args.simple, args.loc))
    expect_that(fits.list, is_identical_to(fits))
})

test_that("simple fitting -- hazard rate", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- admbsecr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, detfn = "hr")
    ## Checking parameter values.
    pars.test <- c(2665.01820347807, 0.879169466112046, 7.04615764143495, 
                   7.47669207224367)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(406.39, 0.063895, 0.64907, 2.7843)
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
    set.seed(1512)
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
    pars.test <- c(2171.50549951556, 5.53458073287418, 53.8758189647376)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(236.34, 0.21195, 7.1895)
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
    pars.test <- c(2308.372990408, 5.33172290979062, 5.61787904103153)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(240.55, 0.17895, 0.7112)
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
                    ss.opts = list(cutoff = 60))
    ## Checking parameter values.
    pars.test <- c(2657.27125943043, 89.0711498188406, 4.12645016662117, 
                   9.5214521578511)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(302.72, 1.519, 0.27205, 0.61166)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(length(fit$suppars), equals(0))
    ## Checking get.par() with SS stuff.
    expect_that(get.par(fit, "all")[c(-4, -5)], equals(coef(fit, c("fitted", "derived"))))
    expect_that(get.par(fit, "b0.ss"), equals(fit$coefficients["b0.ss"]))
    ## Testing parameter phases.
    fit.phase <- admbsecr(capt = ss.capt, traps = example$traps,
                    mask = example$mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    ss.opts = list(cutoff = 60),
                    phases = list(sigma.ss = 2, b0.ss = 3))
    relative.error <- max(abs((coef(fit.phase) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing log-link function.
    fit.log <- admbsecr(capt = ss.capt, traps = example$traps,
                        mask = example$mask,
                        ss.opts = list(cutoff = 60, ss.link = "log"),
                        phase = list(b0.ss = 2), hess = FALSE)
    pars.test <- c(2593.90798693611, 4.543225020375, 0.0652746373675036, 
                   9.07891276875498)
    relative.error <- max(abs((coef(fit.log) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing spherical spreading.
    fit.spherical <- admbsecr(capt = ss.capt, traps = example$traps,
                              mask = example$mask,
                              sv = list(D = 2000, b0.ss = 100, b1.ss = 4, sigma.ss = 15),
                              ss.opts = list(cutoff = 60, ss.link = "spherical"),
                              phase = list(b0.ss = 2), hess = FALSE)
    pars.test <- c(2598.87241252612, 91.3484275416721, 2.45078362490863, 
                   8.96224093008298)
    relative.error <- max(abs((coef(fit.spherical) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("toa fitting", {
    ## Fitting model.
    toa.capt <- example$capt[c("bincapt", "toa")]
    fit <- admbsecr(capt = toa.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2007.91805398409, 5.80251359291497, 0.00177369359452944)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(264.45, 0.33865, 0.00017356)
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
                    ss.opts = list(cutoff = 60))
    ## Checking parameter values.
    pars.test <- c(2433.56128756032, 90.0987641388252, 3.98404182942142, 
                   9.34865299092028, 0.00193398569605552)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(252.43, 1.4504, 0.22561, 0.52189, 0.00018498)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss")))
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
})

test_that("joint bearing/dist fitting", {
    joint.capt <- example$capt[c("bincapt", "bearing", "dist")]
    fit <- admbsecr(capt = joint.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1))
    ## Checking parameter values.
    pars.test <- c(2264.81927702306, 5.39436567560547, 51.1182707112547, 
                   5.02663804138096)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(231.18, 0.16661, 6.8418, 0.49763)
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
    pars.test <- c(2267.7394754986, 5.39011188311111, 10, 0.0560029, 
                   226.77394754986)
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
    ses.test <- c(351.86, 0.42008)
    relative.error <- max(abs((stdEr(fit.hess)[1:2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("directional call fitting", {
    joint.capt <- example$capt[c("bincapt", "ss", "toa")]
    joint.capt <- lapply(joint.capt, function(x) x[41:60, ])
    fit <- admbsecr(capt = joint.capt, traps = example$traps,
                    sv = list(b0.ss = 90, b1.ss = 4, b2.ss = 0.1, sigma.ss = 10),
                    mask = example$mask, ss.opts = list(cutoff = 60))
    ## Checking parameter values.
    pars.test <- c(339.357586915871, 89.7928396498161, 3.56048335373993, 
                   7.60947487477794e-06, 8.22625149900229, 0.00186058818783199)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(87.157, 3.1701, 0.62158, 0.40353, 1.0747, 0.0004036)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss")))
    ## Checking R's ESA calculation.
    relative.error <- (admbsecr:::p.dot(fit, esa = TRUE) - coef(fit, "esa"))/coef(fit, "esa")
    expect_that(relative.error < 1e-4, is_true())
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
    ## Checking fitting with local integration.
    fit <- admbsecr(capt = joint.capt, traps = example$traps,
                    sv = list(b0.ss = 90, b1.ss = 4, b2.ss = 0.1, sigma.ss = 10),
                    mask = example$mask, ss.opts = list(cutoff = 60), local = TRUE)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("fitting heterogeneity in source strengths", {
    set.seed(6434)
    pars <- list(D = 500, b0.ss = 90, b1.ss = 4, sigma.b0.ss = 5, sigma.ss = 10)
    capt <- sim.capt(traps = example$traps, mask = example$mask, detfn = "ss",
                     pars = pars, ss.opts = list(cutoff = 60))
    ## First line used to be:
    ## 0.00000  0.00000  0.00000 65.49092   0.00000 79.19198
    ## For reasons I cannot figure out, it has changed to:
    ## 0.00000  0.00000  0.00000  0.00000  0.00000 78.12121
    ## I hope this is not an issue with sim.capt(); though I ran the
    ## above code with older versions and they gave the same capture
    ## history. I don't know what's going on.
    fit <- admbsecr(capt = capt, traps = example$traps,
                    mask = example$mask, sv = pars,
                    phases = list(b0.ss = 2, sigma.b0.ss = 3,
                        sigma.ss = 4),
                    ss.opts = list(cutoff = 60, het.source = TRUE,
                        n.het.source.quadpoints = 5), hess = FALSE,
                    local = TRUE, optim.opts = list(cbs = 1e8, gbs = 1e8))
    pars.test <- c(731.363920521785, 90.5797487316036, 4.71397727982896, 
                   5.633818295533, 8.09256566662896)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
})
