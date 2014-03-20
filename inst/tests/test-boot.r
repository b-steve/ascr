context("Testing bootstrapping")

test_that("simple model bootstrapping", {
    set.seed(8871)
    t.ser <- system.time({
        boot.fit <- boot.admbsecr(simple.hn.fit, N = 2, n.cores = 1)
    })
    set.seed(8871)
    t.par <- system.time({
        boot.fit <- boot.admbsecr(simple.hn.fit, N = 20, n.cores = 3)
    })
    means.test <- c(2165.6640677, 5.7049935, 0.0616552)
    boot.means <- apply(boot.fit$boot, 2, mean)
    relative.error <- max(abs((boot.means - means.test)/means.test))
    expect_that(relative.error < 1e-4, is_true())
    expect_that(coef(boot.fit), is_identical_to(coef(simple.hn.fit)))
    ses.test <- c(387.9329, 0.4969055, 0.008813962)
    boot.ses <- stdEr(boot.fit, "all")
    relative.error <- max(abs((boot.ses - ses.test)/ses.test))
    expect_that(relative.error < 1e-4, is_true())
})
