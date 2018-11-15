library(ascr)
set.seed(8173)
example.traps <- cbind(rep(c(0, 5), each = 3), rep(c(0, 5, 10), 2))
example.mask <- create.mask(traps = example.traps, buffer = 30, spacing = 1)
example.capt <- sim.capt(traps = example.traps, mask = example.mask,
                         infotypes = c("bearing", "dist", "toa"),
                         detfn = "ss",
                         pars = list(D = 2500, b0.ss = 90, b1.ss = 4,
                             sigma.ss = 10, kappa = 70, alpha = 4,
                             sigma.toa = 0.002), ss.opts = list(cutoff = 60))
simple.capt <- example.capt["bincapt"]
bearing.capt <- example.capt[c("bincapt", "bearing")]
simple.hn.fit <- fit.ascr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, fix = list(g0 = 1))
simple.hr.fit <- fit.ascr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, detfn = "hr", sv = list(z = 5))
bearing.hn.fit <- fit.ascr(capt = bearing.capt, traps = example.traps,
                           mask = example.mask, fix = list(g0 = 1))
boot.simple.hn.fit <- boot.ascr(simple.hn.fit, N = 500, n.cores = 4)
example <- list(capt = example.capt, traps = example.traps, mask = example.mask,
                cutoff = 60, fits = list(simple.hn = simple.hn.fit,
                                 simple.hr = simple.hr.fit,
                                 bearing.hn = bearing.hn.fit,
                                 boot.simple.hn = boot.simple.hn.fit))
traps1 <- example$traps
mask1 <- create.mask(traps1, buffer = 15)
capt1 <- sim.capt(traps = traps1, mask = mask1, infotypes = "bearing",
                  pars = list(D = 2500, g0 = 0.9, sigma = 3, kappa = 50))
traps2 <- (example$traps + 20)[1:3, ]
mask2 <- create.mask(traps2, buffer = 15)
capt2 <- sim.capt(traps = traps2, mask = mask2, infotypes = "bearing",
                  pars = list(D = 2500, g0 = 0.9, sigma = 3, kappa = 50))
traps <- list(traps1, traps2)
mask <- list(mask1, mask2)
capt <- list(capt1, capt2)
multi.example <- list(capt = capt, traps = traps, mask = mask)

save(example, multi.example, file = "../../data/example.RData")
