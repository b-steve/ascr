\name{stdEr.admbsecr.boot}
\alias{stdEr.admbsecr.boot}
\title{Extract standard errors from a bootstrapped admbsecr model object}
\usage{
\method{stdEr}{admbsecr.boot}(object, pars = "fitted", ...)
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
Extracts standard errors for parameters in a model fitted
using \link[admbsecr]{admbsecr}, with a bootstrap procedure
carried out using \link[admbsecr]{boot.admbsecr}.
}

