context("Testing bootstrapping")

test_that("simple model bootstrapping", {
    set.seed(8871)
    boot.fit <- boot.admbsecr(simple.hn.fit, N = 5, prog = FALSE)
    means.test <- c(7.6667469, 1.7383131, 2165.68, 5.705, 0.0616552)
    boot.means <- apply(boot.fit$boot$boots, 2, mean)
    relative.error <- max(abs((boot.means - means.test)/means.test))
    expect_that(relative.error < 1e-4, is_true())
    expect_that(coef(boot.fit), is_identical_to(coef(simple.hn.fit)))
    ses.test <- c(387.9411, 0.4969067, 0.008813962, 0.1886847, 0.08700770)
    boot.ses <- stdEr(boot.fit, "all")
    relative.error <- max(abs((boot.ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
    se.mces.test <- c(0.05604748, 0.02459689, 108.3314, 0.1411242, 0.002527854)
    boot.se.mces <- get.mce(boot.fit, "se")
    relative.error <- max(abs((boot.se.mces - se.mces.test)/se.mces.test))
    expect_that(relative.error < 1e-4, is_true())
    bias.mces.test <- c(0.07585980, 0.03489304, 155.8189, 0.1994087, 0.003538833)
    boot.bias.mces <- get.mce(boot.fit, "bias")
    relative.error <- max(abs((boot.bias.mces - bias.mces.test)/bias.mces.test))
    expect_that(relative.error < 1e-4, is_true())
})

test_that("bootstrapping helpers", {
    ## MCE extraction.
    expect_that(get.mce(boot.simple.hn.fit, "se"),
                equals(boot.simple.hn.fit$boot$se.mce))
    ## Variance-covariance matrix extraction.
    expect_that(sort(vcov(boot.simple.hn.fit, "all")),
                equals(sort(boot.simple.hn.fit$boot$vcov)))
    ## Confidence interval methods.
    ## Default.
    cis.test <- matrix(c(1597.158022, 4.309298, 3120.317156, 6.216009), nrow = 2)
    cis <- confint(boot.simple.hn.fit)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Default linked, also making sure qqplot works.
    cis.test <- matrix(c(1702.471189, 4.400571, 3267.98071, 6.29362), nrow = 2)
    cis <- confint(boot.simple.hn.fit, linked = TRUE, qqplot = TRUE, ask = FALSE)
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Basic (also make sure qqplot works).
    cis.test <- matrix(c(1581.378741, 4.260765, 3043.377993, 6.060344), nrow = 2)
    cis <- confint(boot.simple.hn.fit, method = "basic")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
    ## Percentile.
    cis.test <- matrix(c(1674.097184, 4.464963, 3136.096437, 6.264543), nrow = 2)
    cis <- confint(boot.simple.hn.fit, method = "percentile")
    relative.error <- max(abs((cis - cis.test)/cis.test))
    expect_that(relative.error < 1e-4, is_true())
})
