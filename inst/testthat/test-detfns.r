context("Testing detection probability calculations")

test_that("half normal", {
    detfn <- admbsecr:::get.detfn("hn")
    test.hn <- function(d, g0, sigma){
        g0*exp(-d^2/(2*sigma^2))
    }
    ds <- 0:50
    pars <- list(g0 = 0.6, sigma = 5)
    probs <- detfn(ds, pars)
    test.probs <- test.hn(ds, 0.6, 5)
    expect_that(probs, equals(test.probs))
    pars$z <- 3
    expect_that(detfn(ds, pars),
                throws_error("Argument 'pars' must have named components"))
})

test_that("hazard rate", {
    detfn <- admbsecr:::get.detfn("hr")
    test.hr <- function(d, g0, sigma, z){
        g0*(1 - exp(-(d/sigma)^-z))
    }
    ds <- 0:50
    pars <- list(g0 = 0.6, sigma = 7, z = 3)
    probs <- detfn(ds, pars)
    test.probs <- test.hr(ds, 0.6, 7, 3)
    expect_that(probs, equals(test.probs))
    pars <- list(g0 = 0.6, sigma = 7)
    expect_that(detfn(ds, pars),
                throws_error("Argument 'pars' must have named components"))
})

test_that("threshold", {
    detfn <- admbsecr:::get.detfn("th")
    test.th <- function(d, scale, shape){
        0.5 - 0.5*admbsecr:::erf(d/scale - shape)
    }
    ds <- 0:50
    pars <- list(shape = 6, scale = 3)
    probs <- detfn(ds, pars)
    test.probs <- test.th(ds, 3, 6)
    expect_that(probs, equals(test.probs))
    pars <- list(g0 = 0.6, shape = 7)
    expect_that(detfn(ds, pars),
                throws_error("Argument 'pars' must have named components"))
})

test_that("log-link threshold", {
    detfn <- admbsecr:::get.detfn("lth")
    test.lth <- function(d, scale, shape.1, shape.2){
       0.5 - 0.5*admbsecr:::erf(shape.1 - exp(shape.2 - scale*d))
    }
    ds <- 0:50
    pars <- list(shape.1 = 2, scale = 0.1, shape.2 = 3)
    probs <- detfn(ds, pars)
    test.probs <- test.lth(ds, 0.1, 2, 3)
    expect_that(probs, equals(test.probs))
    pars <- list(g0 = 0.6, shape = 7)
    expect_that(detfn(ds, pars),
                throws_error("Argument 'pars' must have named components"))
})


