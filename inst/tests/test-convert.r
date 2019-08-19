context("Testing conversion functions")

## Note convert.mask() and convert.capt.to.secr() are tested in test-fits.r.

test_that("Converting/creating capture histories", {
    dets <- lightfooti$dets
    traps <- lightfooti$traps
    ## Testing old allocation method.
    capt.test <- convert.pamguard(dets = dets, mics = traps, new.allocation = FALSE)
    ## New convert.pamguard() retains original time data.
    capt.test$toa[capt.test$toa > 0] <- capt.test$toa[capt.test$toa > 0] - (capt.test$toa[1, 2] - 1)
    expect_that(capt.test, is_identical_to(lightfooti$capt))
    ## Testing new allocation method.
    capt.test.new <- convert.pamguard(dets = dets, mics = traps)
    expect_that(capt.test.new, is_identical_to(lightfooti$capt.new))
    mask.test <- create.mask(traps = traps, buffer = 40)
    expect_that(mask.test, is_identical_to(lightfooti$mask))
    ## Creating a multi-session captures data frame.
    captures <- data.frame(session = c(1, 1, 1, 1, 3, 3, 3, 3, 3),
                           id = c(1, 2, 1, 2, 1, 2, 3, 4, 4),
                           occasion = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
                           trap = c(2, 1, 3, 2, 4, 3, 2, 6, 5),
                           bearing = c(5.59814552863741, 1.55785557290341, 0.226989015597898,
                                       6.0614780856239, 4.45397150613391, 4.11997712380155,
                                       5.63614140398099, 4.54698079620703, 3.70723587716335))
    traps <- list(data.frame(x = 1:3, y = rep(0, 3)),
                  data.frame(x = 1:3, y = rep(0, 3)),
                  data.frame(x = rep(1:4, 2), y = rep(c(0, 1), each = 4)),
                  data.frame(x = 1:3, y = rep(0, 3)))
    capt <- create.capt(captures, traps = traps)
    mask <- create.mask(traps, 10)
    ## Fitting a model for testing purposes.
    fit <- fit.ascr(capt = capt, traps = traps, mask = mask, fix = list(g0 = 1, kappa = 0.2))
    expect_equal(coef(fit), expected = c(D = 2393.5239726, sigma = 0.5470687), tol = 1e-5)
})
