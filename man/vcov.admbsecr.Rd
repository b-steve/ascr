\name{vcov.admbsecr}
\alias{vcov.admbsecr}
\title{Extract the variance-covariance matrix from an admbsecr model
object}
\usage{
\method{vcov}{admbsecr}(object, pars = "fitted", ...)
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
Extracts the variance-covariance matrix for parameters in a
model fitted using \link[admbsecr]{admbsecr}.
}
\examples{
vcov(simple.hn.fit)
vcov(simple.hn.fit, pars = "all")
vcov(simple.hn.fit, pars = "derived")
}

