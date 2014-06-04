context("Testing helper functions")

test_that("testing get.par", {
    ## All parameters estimated.
    expect_that(get.par(simple.hr.fit, "all"),
                is_identical_to(coef(simple.hr.fit, c("fitted", "derived"))))
    expect_that(get.par(simple.hr.fit, c("D", "esa")),
                is_identical_to(coef(simple.hr.fit, c("fitted", "derived"))[c(1, 5)]))
    ## Fixed parameter.
    expect_that(get.par(simple.hn.fit, "g0"), is_equivalent_to(1))
    expect_that(get.par(simple.hn.fit, c("sigma", "g0")),
                is_equivalent_to(c(coef(simple.hn.fit)[2], 1)))
    ## Supplementary parameters.
    expect_that(get.par(bearing.hn.fit, c("kappa", "D", "esa", "g0")),
                is_equivalent_to(c(coef(bearing.hn.fit, "all")[c(3, 1, 4)], 1)))
})

test_that("testing calculation of probability detection surface", {
    esa.test <- sum(admbsecr:::p.dot(simple.hn.fit))*attr(example.mask, "area")
    esa <- get.par(simple.hn.fit, "esa")
    relative.error <- (esa.test - esa)/esa
    expect_that(relative.error < 1e-4, is_true())
})
