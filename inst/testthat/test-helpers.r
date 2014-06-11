context("Testing helper functions")

test_that("testing get.par", {
    ## All parameters estimated.
    expect_that(get.par(example$fits$simple.hr, "all"),
                is_identical_to(coef(example$fits$simple.hr, c("fitted", "derived"))))
    expect_that(get.par(example$fits$simple.hr, c("D", "esa")),
                is_identical_to(coef(example$fits$simple.hr, c("fitted", "derived"))[c(1, 5)]))
    ## Fixed parameter.
    expect_that(get.par(example$fits$simple.hn, "g0"), is_equivalent_to(1))
    expect_that(get.par(example$fits$simple.hn, c("sigma", "g0")),
                is_equivalent_to(c(coef(example$fits$simple.hn)[2], 1)))
    ## Supplementary parameters.
    expect_that(get.par(example$fits$bearing.hn, c("kappa", "D", "esa", "g0")),
                is_equivalent_to(c(coef(example$fits$bearing.hn, "all")[c(3, 1, 4)], 1)))
})

test_that("testing calculation of probability detection surface", {
    esa.test <- sum(admbsecr:::p.dot(example$fits$simple.hn))*attr(example$mask, "area")
    esa <- get.par(example$fits$simple.hn, "esa")
    relative.error <- (esa.test - esa)/esa
    expect_that(relative.error < 1e-4, is_true())
})
