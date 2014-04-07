\name{coef.admbsecr}
\alias{coef.admbsecr}
\alias{coef.admbsecr.boot}
\title{Extract admbsecr model coefficients}
\usage{
\method{coef}{admbsecr}(object, pars = "fitted", ...)

\method{coef}{admbsecr.boot}(object, pars = "fitted", correct.bias = FALSE,
  ...)
}
\arguments{
  \item{object}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{pars}{A character vector containing either
  parameter names, or a subset of \code{"all"},
  \code{"derived"}, \code{"fitted"}, and \code{"linked"};
  \code{"fitted"} corresponds to the parameters of
  interest, \code{"derived"} corresponds to quantities that
  are functions of these parameters (e.g., the effective
  survey area or animal density from an acoustic survey),
  and \code{"linked"} corresponds to the parameters AD
  Model Builder has maximised the likelihood over.}

  \item{...}{Other parameters (for S3 generic
  compatibility).}

  \item{correct.bias}{Logical, if \code{TRUE}, estimated
  biases are subtracted from estimated parameter values.}
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

