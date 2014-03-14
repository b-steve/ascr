\name{stdEr.admbsecr}
\alias{stdEr.admbsecr}
\title{Extract standard errors from an admbsecr model fit}
\usage{
\method{stdEr}{admbsecr}(object, pars = "fitted", ...)
}
\arguments{
  \item{object}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{pars}{A character string; either \code{"all"},
  \code{"fitted"}, or \code{"derived"}. If \code{"all"},
  all estimated parameters are returned. If
  \code{"fitted"}, only those parameters for which the
  likelihood was maximised over are returned. If
  \code{"derived"}, only those parameters (e.g., the
  effective survey area) which are functions of fitted
  parameters are returned.}

  \item{...}{Other parameters (for S3 generic
  compatibility).}
}
\description{
Extracts standard errors for estimated and derived
parameters from a model fitted using
\link[admbsecr]{admbsecr}.
}
\examples{
stdEr(simple.hn.fit)
stdEr(simple.hn.fit, pars = "all")
stdEr(simple.hn.fit, pars = "derived")
}

