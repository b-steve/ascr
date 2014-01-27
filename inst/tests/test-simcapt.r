context("Testing simulation code")

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
    expect_that({sim.capt(traps = traps, mask = mask, detfn = "ss",
                          pars = list(D = 2000, ss.b0 = 20, ss.b1 = 5,
                              sigma.ss = 10), ss.link = "identity.link",
                          cutoff= 0)},
                throws_error("The argument 'ss.link' must be either \"identity\" or \"log\""))
})
