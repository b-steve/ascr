\name{get.bias}
\alias{get.bias}
\title{Extracting estimated biases}
\usage{
get.bias(fit, pars = "fitted")
}
\arguments{
  \item{fit}{A fitted model from \link{admbsecr}.}

  \item{pars}{A character vector containing either
  parameter names, or a subset of \code{"all"},
  \code{"derived"}, \code{"fitted"}, and \code{"linked"};
  \code{"fitted"} corresponds to the parameters of
  interest, \code{"derived"} corresponds to quantities that
  are functions of these parameters (e.g., the effective
  survey area or animal density from an acoustic survey),
  and \code{"linked"} corresponds to the parameters AD
  Model Builder has maximised the likelihood over.}
}
\description{
Extracts bias in parameter estimates, calculated using the
bootstrap procedure carried out by \link{boot.admbsecr}.
}
\seealso{
\link{boot.admbsecr} for the bootstrap procedure.

\link{get.mce} for Monte Carlo error the biases are subject
to.
}

