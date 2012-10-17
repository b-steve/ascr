\name{secrlikelihood}
\alias{secrlikelihood}
\title{Likelihood function for simple SECR.}
\usage{
  secrlikelihood(beta, capthist, mask, dist = NULL,
    trace = FALSE)
}
\arguments{
  \item{beta}{beta.}

  \item{capthist}{capthist.}

  \item{mask}{mask.}

  \item{dist}{dist.}

  \item{trace}{trace.}
}
\value{
  The negative log-likelihood.
}
\description{
  Calculates the likelihood for a simple SECR model given a
  set of parameter values.
}
\author{
  Murray Efford, David Borchers, Ben Stevenson
}

