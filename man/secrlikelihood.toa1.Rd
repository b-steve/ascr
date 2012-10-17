\name{secrlikelihood.toa1}
\alias{secrlikelihood.toa1}
\title{Likelihood function for SECR with TOA information.}
\usage{
  secrlikelihood.toa1(beta, capthist, mask, dists = NULL,
    ssqtoa = NULL, trace = FALSE)
}
\arguments{
  \item{beta}{beta.}

  \item{capthist}{capthist.}

  \item{mask}{mask.}

  \item{dists}{dists.}

  \item{ssqtoa}{ssqtoa.}

  \item{trace}{trace.}
}
\value{
  The negative log-likelihood.
}
\description{
  Calculates the likelihood for a SECR model with TOA
  information given a set of parameter values.
}
\author{
  David Borchers, Murray Efford, Ben Stevenson
}

