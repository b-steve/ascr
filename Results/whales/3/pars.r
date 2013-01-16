## Same as (2) with alpha set at 9.8. (32% CV)
seed <- 8465
D <- 1.7
g01 <- 0.99999
sigma1 <- 210
g02 <- 0.30
sigma2 <- 320
alpha <- 9.8
truepars <- c(D = D, g01 = g01, sigma1 = sigma1, g02 = g02,
              sigma2 = sigma2, alpha = alpha)
detectpars <- list(g01 = g01, sigma1 = sigma1, g02 = g02,
                   sigma2 = sigma2, alpha = alpha)
