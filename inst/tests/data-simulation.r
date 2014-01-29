set.seed(8173)
example.capt <- sim.capt(traps = example.traps, mask = example.mask,
                         infotypes = c("ang", "dist", "toa"),
                         detfn = "ss",
                         pars = list(D = 2500, b0.ss = 90, b1.ss = 4,
                             sigma.ss = 10, kappa = 70, alpha = 4,
                             sigma.toa = 0.002), cutoff = 60)    
##save(example.traps, example.mask, example.capt,
##     file = "~/admbsecr/data/simdata.RData")
