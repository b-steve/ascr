distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

secrlikelihood <- function (beta, capthist, mask, dist = NULL, trace = FALSE) {
  ## Compute negative log likelihood for halfnormal proximity model
  ## Murray Efford 2011-02-27
  
  ## Inputs
  ##     beta  -  parameter vector on link scale (D, g0, sigma)
  ##     capthist - capthist object (n x S x K binary array)
  ##     mask - mask object (M x 2 matrix of x-y coords, with attribute 'area')
  ##     dist - K x M matrix of distances between each detector and each mask point
  ##     trace - logical TRUE for output of each evaluation
  
  ## where K = number of detectors, M = number of mask points
  
  if (!all(capthist %in% c(0,1)))
    stop ('secrlikelihood requires binary data')
  
  ## 'real' parameter values
  D <- exp(beta[1])
  ## if D is modelled, expand here to a vector of length
  ## equal to number of rows in mask
  g0 <- invlogit(beta[2])
  sigma <- exp(beta[3])
  
  n <- nrow(capthist)        ## number observed
  S <- ncol(capthist)        ## number of occasions
  A <- attr( mask, 'area')   ## area of one cell
  
  ## distances if not passed as argument
  if (is.null(dist)) {
    traps <- traps(capthist)
    dist <- distances(traps, mask)
  }
  
  ## precompute probabilities for each detector and each mask point
  ## gk is K x M matrix
  gk <- g0 * exp(-dist^2 / 2 / sigma^2)
  log.gk <- log(gk)
  log.gk1 <- log(1-gk)
  
  ## probability of being caught at least once if at mask site m
  p.m <- 1 - apply(1-gk, 2, prod) ^ S  ## vector of length M
  sumDp <- sum(p.m * D)                ## scalar
  
  L1 <- sum ( apply(capthist, 1, Dprwi) )
  L2 <- - n * log(sumDp)
  L3 <- dpois (n, sumDp * A, log = TRUE)
  LL <- L1+L2+L3
  if (trace) {
    cat (D, g0, sigma, LL, '\n')
    flush.console()
  }
  if (!is.finite(LL))
    1e10
  else
    -LL   ## return negative log likelihood
}

Dprwi <- function (wi) {
  ## wi is S x K, log.gk is K x M, prwi.s is S x M
  prwi.s <- wi %*% log.gk + (1-wi) %*% log.gk1
  ## sum log(p) over occasions, result is M-vector
  prwi.s <- apply(prwi.s,2,sum)
  log (sum(D * exp(prwi.s)))
}
