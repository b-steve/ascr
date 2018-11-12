context("Testing model fits")

test_that("simple fitting -- half normal", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
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
    expect_that(get.par(fit, "esa"), equals(fit$coefficients["esa.1"]))
    expect_that(get.par(fit, c("D", "esa")),
                equals(fit$coefficients[c("D", "esa.1")]))
    expect_that(get.par(fit, "all"), equals(coef(fit, c("fitted", "derived"))))
    expect_that(get.par(fit, "fitted"), equals(coef(fit, "fitted")))
    ## Testing some generic functions.
    expect_that(is.list(summary(fit)), is_true())
    expect_that(confint(fit), is_a("matrix"))
    ## Checking hess argument.
    fit.nohess <- fit.ascr(capt = simple.capt, traps = example$traps,
                           mask = example$mask, hess = FALSE)
    expect_that(coef(fit.nohess), equals(coef(fit)))
    ## Checking estimate equivalence with secr.fit.
    suppressPackageStartupMessages(library(secr))
    mask.secr <- convert.mask(example$mask)
    capt.secr <- convert.capt.to.secr(example$capt["bincapt"], example$traps)
    fit.secr <- suppressWarnings(secr.fit(capthist = capt.secr, mask = mask.secr, trace = FALSE))
    coefs.secr <- numeric(n.pars)
    invlog <- function(x) exp(x)
    for (i in 1:n.pars){
        coefs.secr[i] <- eval(call(paste("inv", fit.secr$link[[i]], sep = ""),
                                   fit.secr$fit$par[i]))
    }
    relative.error <- max(abs((coef(fit) - coefs.secr)/coefs.secr))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking fitting with local integration.
    fit.loc <- fit.ascr(capt = simple.capt, traps = example$traps,
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
    fits <- par.fit.ascr(n.cores = 2, args.simple, args.loc)
    expect_that(summary(fits[[1]]), is_identical_to(summary(fit)))
    expect_that(summary(fits[[2]]), is_identical_to(summary(fit.loc)))
})

test_that("simple fitting -- hazard rate", {
    ## Fitting model.
    simple.capt <- example$capt["bincapt"]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, detfn = "hr", sv = list(z = 5))
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
    capt.secr <- convert.capt.to.secr(example$capt["bincapt"], example$traps)
    set.seed(1512)
    fit.secr <- suppressWarnings(secr.fit(capthist = capt.secr, mask = mask.secr, detectfn = 1,
                                          trace = FALSE))
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
    fit <- fit.ascr(capt = bearing.capt, traps = example$traps,
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
    fit <- fit.ascr(capt = dist.capt, traps = example$traps,
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
    fit <- fit.ascr(capt = ss.capt, traps = example$traps,
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
    fit.phase <- fit.ascr(capt = ss.capt, traps = example$traps,
                    mask = example$mask,
                    sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10),
                    ss.opts = list(cutoff = 60),
                    phases = list(sigma.ss = 2, b0.ss = 3))
    relative.error <- max(abs((coef(fit.phase) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing log-link function.
    fit.log <- fit.ascr(capt = ss.capt, traps = example$traps,
                        mask = example$mask,
                        ss.opts = list(cutoff = 60, ss.link = "log"),
                        phase = list(b0.ss = 2), hess = FALSE)
    pars.test <- c(2593.90798693611, 4.543225020375, 0.0652746373675036, 
                   9.07891276875498)
    relative.error <- max(abs((coef(fit.log) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing spherical spreading.
    fit.spherical <- fit.ascr(capt = ss.capt, traps = example$traps,
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
    fit <- fit.ascr(capt = toa.capt, traps = example$traps,
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
                equals(fit$coefficients[c("sigma.toa", "esa.1")]))
})

test_that("joint ss/toa fitting", {
    joint.capt <- example$capt[c("bincapt", "ss", "toa")]
    fit <- fit.ascr(capt = joint.capt, traps = example$traps,
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
    fit <- fit.ascr(capt = joint.capt, traps = example$traps,
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
    ## Checking things work with survey.length != 1.
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1),
                    cue.rates = c(9, 10, 11), survey.length = 2)
    pars.test <- c(113.38697377493, 5.39011188311111)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit, c("Da", "sigma")) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1),
                    cue.rates = c(9, 10, 11), survey.length = 1)
    pars.test <- c(2267.7394754986, 5.39011188311111, 10, 0.0560029, 
                   226.77394754986, 2267.7394754986)
    n.pars <- length(pars.test)
    relative.error <- max(abs((coef(fit, c("fitted", "derived")) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    expect_that(confint(fit),
                throws_error("Standard errors not calculated; use boot.ascr()"))
    expect_that(summary(fit), is_a("list"))
    ## Checking hess argument.
    fit.hess <- fit.ascr(capt = simple.capt, traps = example$traps,
                         mask = example$mask, fix = list(g0 = 1),
                         cue.rates = c(9, 10, 11), survey.length = 1, hess = TRUE)
    expect_that(coef(fit.hess), equals(coef(fit)))
    expect_that(is.na(stdEr(fit.hess, "all")["Da"]), is_true())
    ses.test <- c(351.86, 0.42008)
    relative.error <- max(abs((stdEr(fit.hess)[1:2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking simulation generates multiple cues per location.
    set.seed(2125)
    test.capt <- sim.capt(fit, keep.locs = TRUE)
    expect_that(nrow(test.capt$popn.locs), equals(917))
    expect_that(length(unique(test.capt$popn.locs[, 1])), equals(92))
    expect_that(nrow(test.capt$capt.locs), equals(134))
    expect_that(length(unique(test.capt$capt.locs[, 1])), equals(24))
    set.seed(3429)
    fit.boot <- boot.ascr(fit = fit, N = 10, prog = FALSE)
    ses.test <- c(886.9398, 1.1567)
    relative.error <- max(abs((stdEr(fit.boot)[1:2] - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("directional call fitting", {
    joint.capt <- example$capt[c("bincapt", "ss", "toa")]
    joint.capt <- lapply(joint.capt, function(x) x[41:60, ])
    fit <- fit.ascr(capt = joint.capt, traps = example$traps,
                    sv = list(b0.ss = 90, b1.ss = 4, b2.ss = 0.1, sigma.ss = 10),
                    mask = example$mask, ss.opts = list(cutoff = 60))
    ## Checking parameter values.
    pars.test <- c(339.357584343541, 89.7928412454349, 3.56047013041414, 
                   1.363056e-05, 8.22625124933556, 0.00186058803701272)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking standard errors.
    ses.test <- c(87.157, 3.1701, 0.6174500, 0.79429, 1.0747, 0.0004036)
    relative.error <- max(abs((stdEr(fit) - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Checking detection parameters.
    expect_that(fit$detpars, equals(c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss")))
    ## Checking R's ESA calculation.
    relative.error <- (ascr:::p.dot(fit, esa = TRUE) - coef(fit, "esa"))/coef(fit, "esa")
    expect_that(relative.error < 1e-4, is_true())
    ## Checking supplementary parameters.
    expect_that(fit$suppars, equals("sigma.toa"))
})

test_that("fitting heterogeneity in source strengths", {
    set.seed(6434)
    pars <- list(D = 500, b0.ss = 90, b1.ss = 4, sigma.b0.ss = 5, sigma.ss = 10)
    capt <- sim.capt(traps = example$traps, mask = example$mask, detfn = "ss",
                     pars = pars, ss.opts = list(cutoff = 60))
    ## Using only the first ten observations.
    capt <- lapply(capt, function(x) x[1:10, ])
    ## First line used to be:
    ## 0.00000  0.00000  0.00000 65.49092   0.00000 79.19198
    ## For reasons I cannot figure out, it has changed to:
    ## 0.00000  0.00000  0.00000  0.00000  0.00000 78.12121
    ## I hope this is not an issue with sim.capt(); though I ran the
    ## above code with older versions and they gave the same capture
    ## history. I don't know what's going on.
    fit <- fit.ascr(capt = capt, traps = example$traps,
                    mask = example$mask, sv = pars,
                    phases = list(b0.ss = 2, sigma.b0.ss = 3,
                        sigma.ss = 4),
                    ss.opts = list(cutoff = 60, het.source = TRUE,
                        n.het.source.quadpoints = 5), hess = FALSE,
                    local = TRUE)
    pars.test <- c(221.900466410203, 88.8289943264761, 4.07833236335682, 
                   7.11871056989856, 6.09904893903392)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("first-call signal strength models", {
    set.seed(8298)
    traps <- cbind(c(0, 0, 21, 21, 42, 42), c(0, 42, 0, 42, 0, 42))
    mask <- create.mask(traps, buffer = 1250, spacing = 25)
    pars <- list(D = 5, b0.ss = 60, b1.ss = 0.1, sigma.ss = 5)
    lower.cutoff <- 52.5
    cutoff <- 55
    capt <- sim.capt(traps = traps, mask = mask, detfn = "ss",
                     infotypes = NULL, pars = pars,
                     ss.opts = list(cutoff = lower.cutoff,
                         ss.link = "identity"),
                     cue.rates = Inf, first.only = TRUE)
    fit <- fit.ascr(capt = capt, traps = traps, mask = mask,
                    ss.opts = list(cutoff = cutoff,
                                   lower.cutoff = lower.cutoff),
                    hess = TRUE)
    pars.test <- c(2.6061075, 62.0208816041957, 0.11497237446662, 
                   6.24356345437168)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-3, is_true())
    ## Testing for an error if identity ss.link is not used.
    expect_that(fit <-  fit.ascr(capt = capt, traps = traps, mask = mask,
                     ss.opts = list(cutoff = cutoff,
                         lower.cutoff = lower.cutoff, ss.link = "log"), hess = FALSE),
                throws_error("First-call models are only implemented for ss.link = \"identity\"."))
})

test_that("Multi-session models", {
    fit <- fit.ascr(multi.example$capt, multi.example$traps, multi.example$mask,
                    sv = list(kappa = 100))
    pars.test <- c(2525.4060484, 0.9849349, 2.8395646, 120.0353505)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-3, is_true())
    set.seed(2987)
    boot.fit <- boot.ascr(fit, N = 10, prog = FALSE)
    ses <- stdEr(boot.fit, "all")
    ses.test <- c(325.117358182622, 0.0708826944902274,
                  0.182620570101098, 14.6858155989752,
                  0.00190046466718596, 0.00137547230211622,
                  0.121005550405832, 7.78842978292226,
                  0.0650553077849548, 0.233809361964724)
    relative.error <- max(abs((ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("Inhomogeneous density estimation", {
    simple.capt <- example$capt[1]
    df <- data.frame(x = example$mask[, 1]/1000, y = example$mask[, 2]/1000)
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1),
                    ihd.opts = list(model = ~ x + y, covariates = df))
    pars.test <- c(7.72626004368, 0.0562992293496, 0.088380677303, 5.3875700735678)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-3, is_true())
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, ihd.opts = list(model = ~ s(x, k = 3)))
    pars.test <- c(5.00442908703, 12.0411342793, 0.0626406907173, 0.999999999974208,
                   5.08547672171856)
    relative.error <- max(abs((coef(fit) - pars.test)/pars.test))
    expect_that(relative.error < 1e-3, is_true())
})
