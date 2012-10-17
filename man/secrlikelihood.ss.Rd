\name{secrlikelihood.ss}
\alias{secrlikelihood.ss}
\title{Likelihood function for SECR with signal strength information.}
\usage{
  secrlikelihood.ss(beta, capthist, mask, dists = NULL,
    cutoff, trace = FALSE)
}
\arguments{
  \item{beta}{beta.}

  \item{capthist}{capthist.}

  \item{mask}{mask.}

  \item{dists}{dists.}

  \item{cutoff}{cutoff.}

  \item{trace}{trace.}
}
\value{
  The negative log-likelihood.
}
\description{
  Calculates the likelihood for a SECR model with signal
  strength information given a set of parameter values.
}
\author{
  Ben Stevenson, Murray Efford
}

