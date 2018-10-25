context("Testing bootstrapping")

test_that("simple model bootstrapping", {
    set.seed(8871)
    ## Carrying out bootstrap.
    boot.fit <- boot.ascr(example$fits$simple.hn, N = 5, prog = FALSE)
    ## Mean bootstrap parmeter values.
    means.test <- c(7.84545539423, 1.645086151006, 7.84545539423, 5.19314936078642, 
                    0.05289568, 2579.83223503568)
    boot.means <- apply(boot.fit$boot$boots, 2, mean)
    relative.error <- max(abs((boot.means - means.test)/means.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Coefficients.
    expect_that(coef(boot.fit), is_identical_to(coef(example$fits$simple.hn)))
    ## Standard errors.
    ses.test <- c(399.379618344681, 0.392493876071427, 0.00651806741181771, 
                  0.160164821510471, 0.0748244587754581)
    boot.ses <- stdEr(boot.fit, "all")
    relative.error <- max(abs((boot.ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Monte Carlo error calculation.
    se.mces.test <- c(0.0421420801420685, 0.0190421463886411, 0.0421420801420685, 
                      0.102045193521712, 0.00173380322371648, 99.0170592323127)
    boot.se.mces <- ascr:::get.mce(boot.fit, "se")
    relative.error <- max(abs((boot.se.mces - se.mces.test)/se.mces.test))
    expect_that(relative.error < 1e-4, is_true())
    bias.mces.test <- c(0.0636401620654382, 0.0298417966775328, 0.0636401620654382, 
                        0.156547059873566, 0.002599935557669, 158.490221019319)
    boot.bias.mces <- ascr:::get.mce(boot.fit, "bias")
    relative.error <- max(abs((boot.bias.mces - bias.mces.test)/bias.mces.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing parallel bootstrap.
    library(parallel)
    set.seed(8871)
    boot.fit.par <- boot.ascr(example$fits$simple.hn, N = 5, prog = FALSE,
                                  n.cores = detectCores())
    expect_that(boot.fit.par, is_identical_to(boot.fit))
})

## Note bootstrapping of multiple-cue models is in test-fits.r
## Also bootstrapping of multiple-session models is there too.

test_that("cue-rate model bootstrapping", {
    set.seed(6324)          
    fit <- fit.ascr(capt = lapply(lightfooti$capt, function(x) x[1:20, ]),
                    traps = lightfooti$traps, mask = lightfooti$mask,
                    cue.rates = lightfooti$freqs, survey.length = 1/5,
                    ss.opts = list(cutoff = lightfooti$cutoff))
    boot.fit <- boot.ascr(fit = fit, N = 10, prog = FALSE)
    ses.test <- c(305.466627386732, 5.03347019909088,
                  0.861438881545055, 1.24394224402644,
                  0.00028875810159441, 0.0982092751647983,
                  0.0384548858779046, 45.1831254619258,
                  305.466627386732, 0.399775773640889,
                  0.0291562697985857, 0.223637444167788,
                  0.108387339783415, 0.350108662810472)
    ses <- stdEr(boot.fit, "all")
    relative.error <- max(abs((ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("bootstrapping helpers", {
    ## MCE extraction.
    expect_that(ascr:::get.mce(example$fits$boot.simple.hn, "se"),
                equals(example$fits$boot.simple.hn$boot$se.mce))
    ## Variance-covariance matrix extraction.
    expect_that(sort(vcov(example$fits$boot.simple.hn, "all")[c(4, 5, 4, 2, 3, 1), c(4, 5, 4, 2, 3, 1)]),
                equals(sort(example$fits$boot.simple.hn$boot$vcov)))
    ## Extraction of bias.
    expect_that(get.bias(example$fits$boot.simple.hn, c("D", "esa")),
                equals(example$fits$boot.simple.hn$boot$bias[c("D", "esa.1")]))
    ## Coefficient extraction with bias correction.
    expect_that(coef(example$fits$boot.simple.hn, correct.bias = TRUE),
                equals(coef(example$fits$boot.simple.hn) - get.bias(example$fits$boot.simple.hn)))
    ## Confidence interval methods.
    ## Default.
    cis.test <- matrix(c(1420.53146733787, 4.3550748265765, 3114.94748365933, 6.42514893964572),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Default linked, also making sure qqplot works.
    cis.test <- matrix(c(1559.57881942082, 4.45279026517435, 3297.45586737609, 6.52474165237109),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, linked = TRUE, qqplot = TRUE, ask = FALSE)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Default with bias correction.
    cis.test <- matrix(c(1402.08499580554, 4.33572108417459, 3096.50101212701, 6.40579519724382),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "default.bc")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Basic.
    cis.test <- matrix(c(1421.96795794811, 4.2850913343052, 3065.72978328319, 6.31246947068946),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "basic")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Percentile.
    cis.test <- matrix(c(1469.74916771401, 4.46775429553276, 3113.51099304909, 6.49513243191702),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "percentile")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
})
