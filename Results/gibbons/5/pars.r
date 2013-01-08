## Redoing simulation set 4, but estimating g0.
## Hopefully this will show this is why there is no bias in simulation set 3.
set.seed <- 7862
buffer <- 6000
mask.spacing <- 50
trap.spacing <- 500

## True parameter values.
D <- 0.0416
g0 <- 1
sigma <- 1250
kappa <- 70
truepars <- c(D = D, g0 = g0, sigma = sigma, kappa = kappa)
detectpars <- list(g0 = g0, sigma = sigma)

