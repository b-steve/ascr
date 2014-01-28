context("Testing simulation errors")

test_that("Test output", {
    set.seed(8172)
    traps <- cbind(c(0, 1, 0, 1),
                   c(0, 0, 1, 1))
    colnames(traps) <- c("x", "y")
    mask <- create.mask(traps, buffer = 30)
    capt <- sim.capt(traps = traps, mask = mask,
                     pars = list(D = 2000, g0 = 0.75, sigma = 5))
    capt.test <- list(bincapt = structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0,
                          1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0,
                          1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
                          1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0,
                          1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
                          0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0,
                          1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
                          1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                          1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1,
                          1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0,
                          1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1,
                          1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0,
                          0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1,
                          1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0,                          
                          0, 0, 1, 0, 1, 0, 1, 0, 1, 0), .Dim = c(61L,
                                                             4L)))
    expect_that(capt, is_a("list"))
    expect_that(length(capt), equals(1))
    expect_that(names(capt), equals("bincapt"))
    expect_that(dim(capt$bincapt), equals(c(61, 4)))
    expect_that(capt, equals(capt.test))
})

test_that("Test errors", {
    set.seed(8172)
    traps <- cbind(c(0, 1, 0, 1),
                   c(0, 0, 1, 1))
    colnames(traps) <- c("x", "y")
    mask <- create.mask(traps, buffer = 30)
    ## Specifying incorrect SS link.
    expect_that(sim.capt(traps = traps, mask = mask, detfn = "ss",
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10), ss.link = "identity.link",
                         cutoff= 0),
                throws_error("The argument 'ss.link' must be either \"identity\" or \"log\""))
    ## Trying to simulate SS using 'infotypes'.
    expect_that(sim.capt(traps = traps, mask = mask, infotypes = c("ss", "toa"),
                         pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                             sigma.ss = 10, sigma.toa = 0.002), cutoff= 0),
                throws_error("Signal strength information is simulated by setting argument 'detfn' to \"ss\"."))
    ## Cutoff provided when not required.
    expect_that(sim.capt(traps = traps, mask = mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         cutoff = 0),
                gives_warning("The argument 'cutoff' is being ignored, as 'detfn' is not \"ss\"."))
    ## Signal strength link provided when not required.
    expect_that(sim.capt(traps = traps, mask = mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5),
                         ss.link = "identity"),
                gives_warning("The argument 'ss.link' is being ignored, as 'detfn' is not \"ss\"."))
    ## Extra parameter.
    expect_that(sim.capt(traps = traps, mask = mask,
                         pars = list(D = 2000, g0 = 0.75, sigma = 5,
                             sigma.toa = 0.002)),
                throws_error("The following must be named components of the list 'pars': "))
    ## Missing parameter.
    expect_that(sim.capt(traps = traps, mask = mask,
                         pars = list(D = 2000, g0 = 0.75)),
                throws_error("The following must be named components of the list 'pars': "))
})

load_all("~/admbsecr")
testfun <- function(x = 5){
    print(missing(x))
    x
}
