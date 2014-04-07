library(admbsecr)
set.seed(8173)
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
save(example.traps, example.mask, example.capt,
     file = "/scratch/admbsecr/data/example_data.RData")
save(simple.hn.fit, simple.hr.fit, bearing.hn.fit, boot.simple.hn.fit,
     file = "/scratch/admbsecr/data/example_fits.RData")
