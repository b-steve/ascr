\name{stdEr.admbsecr}
\alias{stdEr.admbsecr}
\title{Extract standard errors from an admbsecr model fit}
\usage{
\method{stdEr}{admbsecr}(object, pars = "fitted", ...)
}
\arguments{
  \item{object}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{pars}{A character containing a subset of
  \code{"all"}, \code{"derived"}, \code{"fitted"}, and
  \code{"linked"}; \code{"fitted"} corresponds to the
  parameters of interest, \code{"derived"} corresponds to
  quantities that are functions of these parameters (e.g.,
  the effective survey area or animal density from an
  acoustic survey), and \code{"linked"} corresponds to the
  parameters AD Model Builder has maximised the likelihood
  over.}

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

