\name{secrlikelihood.angs}
\alias{secrlikelihood.angs}
\title{Likelihood function for SECR with angle information.}
\usage{
  secrlikelihood.angs(beta, capthist, mask, dists = NULL,
    angs = NULL, trace = FALSE)
}
\arguments{
  \item{beta}{beta.}

  \item{capthist}{capthist.}

  \item{mask}{mask.}

  \item{dists}{dists.}

  \item{angs}{angs.}

  \item{trace}{trace.}
}
\value{
  The negative log-likelihood.
}
\description{
  Calculates the likelihood for a SECR model with angle
  information given a set of parameter values.
}
\author{
  Darren Kidney, Murray Efford, Ben Stevenson
}

