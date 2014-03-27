\name{confint.admbsecr}
\alias{confint.admbsecr}
\alias{confint.admbsecr.boot}
\title{Confidence intervals for admbsecr model parameters}
\usage{
\method{confint}{admbsecr}(object, parm = "fitted", level = 0.95, ...)

\method{confint}{admbsecr.boot}(object, parm = "fitted", level = 0.95,
  method = "default", ...)
}
\arguments{
  \item{parm}{A character vector specifying which
  parameters are to be given confidence intervals.}

  \item{method}{A character string specifying the method
  used to calculate the confidence intervals. See 'Details'
  below.}

  \item{object}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{...}{Other parameters (for S3 generic
  compatibility).}

  \item{level}{the confidence level required.}
}
\description{
Computes confidence intervals for one or more parameters
estimated in an admbsecr model object.
}
\details{
Options for the argument \code{method} are as follows:
\code{"default"} for intervals based on a normal
approximation using the calculated standard errors (for
objects of class \code{admbsecr.boot}, these standard
errors are calculated from the bootstrap procedure);
"linked" for intervals that are calculated on the
parameters' link scales, then transformed back onto their
"real" scales; \code{"basic"} for the so-called "basic"
bootstrap method; and \code{"percentile"} for intervals
calculated using the bootstrap percentile method (the
latter two are only available for objects of class
\code{admbsecr.boot}).
}

