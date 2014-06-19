context("Testing bootstrapping")

test_that("simple model bootstrapping", {
    set.seed(8871)
    ## Carrying out bootstrap.
    boot.fit <- boot.admbsecr(example$fits$simple.hn, N = 5, prog = FALSE)
    ## Mean bootstrap parmeter values.
    means.test <- c(7.834395656476, 1.676563641196, 2546.07521272036, 
                    5.3532285309044, 0.0554772512145224)
    boot.means <- apply(boot.fit$boot$boots, 2, mean)
    relative.error <- max(abs((boot.means - means.test)/means.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Coefficients.
    expect_that(coef(boot.fit), is_identical_to(coef(example$fits$simple.hn)))
    ## Standard errors.
    ses.test <- c(349.407607902095, 0.290013693505804, 0.00500991908310524, 
                  0.142526039268406, 0.0528612036412651)
    boot.ses <- stdEr(boot.fit, "all")
    relative.error <- max(abs((boot.ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Monte Carlo error calculation.
    se.mces.test <- c(0.0463009660702605, 0.0188408387410946, 108.747094898535, 
                      0.105149289710528, 0.00184488909656817)
    boot.se.mces <- get.mce(boot.fit, "se")
    relative.error <- max(abs((boot.se.mces - se.mces.test)/se.mces.test))
    expect_that(relative.error < 1e-4, is_true())
    bias.mces.test <- c(0.0560663040454648, 0.0207958537486876, 137.685386201599, 
                        0.114050114609312, 0.00196952221997958)
    boot.bias.mces <- get.mce(boot.fit, "bias")
    relative.error <- max(abs((boot.bias.mces - bias.mces.test)/bias.mces.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Testing parallel bootstrap.
    library(parallel)
    set.seed(8871)
    boot.fit.par <- boot.admbsecr(example$fits$simple.hn, N = 5, prog = FALSE,
                                  n.cores = detectCores())
    expect_that(boot.fit.par, is_identical_to(boot.fit))
})

test_that("bootstrapping helpers", {
    ## MCE extraction.
    expect_that(get.mce(example$fits$boot.simple.hn, "se"),
                equals(example$fits$boot.simple.hn$boot$se.mce))
    ## Variance-covariance matrix extraction.
    expect_that(sort(vcov(example$fits$boot.simple.hn, "all")),
                equals(sort(example$fits$boot.simple.hn$boot$vcov)))
    ## Extraction of bias.
    expect_that(get.bias(example$fits$boot.simple.hn, c("D", "esa")),
                equals(example$fits$boot.simple.hn$boot$bias[c("D", "esa")]))
    ## Coefficient extraction with bias correction.
    expect_that(coef(example$fits$boot.simple.hn, correct.bias = TRUE),
                equals(coef(example$fits$boot.simple.hn) - get.bias(example$fits$boot.simple.hn)))
    ## Confidence interval methods.
    ## Default.
    cis.test <- matrix(c(1613.95919814479, 4.37162348430638, 3245.29373281852, 6.35676252563208),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Default linked, also making sure qqplot works.
    cis.test <- matrix(c(1724.78300492811, 4.46716749063425, 3422.50865465531, 6.4413449137264),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, linked = TRUE, qqplot = TRUE, ask = FALSE)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Default with bias correction.
    cis.test <- matrix(c(1594.05673087843, 4.36168282533297, 3225.39126555216, 6.34682186665868),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "default.bc")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Basic.
    cis.test <- matrix(c(1590.06206614825, 4.24843098809612, 3166.72911351503, 6.22751890142145),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "basic")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Percentile.
    cis.test <- matrix(c(1692.52381744828, 4.50086710851701, 3269.19086481506, 6.47995502184234),
                       nrow = 2)
    cis <- confint(example$fits$boot.simple.hn, method = "percentile")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
})
