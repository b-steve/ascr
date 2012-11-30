## Using identity link function for signal strengths.

seed <- 5253
D <- 4450
g0 <- 0.99999
sigma <- 5.60
sigmatoa <- 0.002
ssb0 <- 170
ssb1 <- -2.50
sigmass <- 7.00
cutoff <- 150
truepars <- c(D = D, g0 = g0, sigma = sigma, sigmatoa = sigmatoa,
              ssb0 = ssb0, ssb1 = ssb1, sigmass = sigmass)
detectpars <- list(beta0 = ssb0, beta1 = ssb1, sdS = sigmass, cutval = cutoff)
