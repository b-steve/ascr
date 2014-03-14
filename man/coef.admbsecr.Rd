\name{coef.admbsecr}
\alias{coef.admbsecr}
\title{Extract admbsecr model coefficients}
\usage{
\method{coef}{admbsecr}(object, pars = "fitted", ...)
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
Extracts estimated and derived parameters from a model
fitted using \link{admbsecr}.
}
\examples{
coef(simple.hn.fit)
coef(simple.hn.fit, pars = "all")
coef(simple.hn.fit, pars = "derived")
}

