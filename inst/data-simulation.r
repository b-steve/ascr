library(admbsecr)
set.seed(8173)
example.traps <- cbind(rep(c(0, 5), each = 3), rep(c(0, 5, 10), 2))
example.mask <- create.mask(traps = example.traps, buffer = 30, spacing = 1)
example.capt <- sim.capt(traps = example.traps, mask = example.mask,
                         infotypes = c("bearing", "dist", "toa"),
                         detfn = "ss",
                         pars = list(D = 2500, b0.ss = 90, b1.ss = 4,
                             sigma.ss = 10, kappa = 70, alpha = 4,
                             sigma.toa = 0.002), cutoff = 60)
simple.capt <- example.capt["bincapt"]
bearing.capt <- example.capt[c("bincapt", "bearing")]
simple.hn.fit <- admbsecr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, fix = list(g0 = 1))
simple.hr.fit <- admbsecr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, detfn = "hr")
bearing.hn.fit <- admbsecr(capt = bearing.capt, traps = example.traps,
                           mask = example.mask, fix = list(g0 = 1))
boot.simple.hn.fit <- boot.admbsecr(simple.hn.fit, N = 500, n.cores = 4)
example <- list(capt = example.capt, traps = example.traps, mask = example.mask,
                cutoff = 60, fits = list(simple.hn = simple.hn.fit,
                                 simple.hr = simple.hr.fit,
                                 bearing.hn = bearing.hn.fit,
                                 boot.simple.hn = boot.simple.hn.fit))
save(example, file = "/scratch/admbsecr/data/example.RData")
