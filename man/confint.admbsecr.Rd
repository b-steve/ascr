\name{confint.admbsecr}
\alias{confint.admbsecr}
\alias{confint.admbsecr.boot}
\title{Confidence intervals for admbsecr model parameters}
\usage{
\method{confint}{admbsecr}(object, parm = "fitted", level = 0.95,
  linked = FALSE, ...)

\method{confint}{admbsecr.boot}(object, parm = "fitted", level = 0.95,
  method = "default", linked = FALSE, qqplot = FALSE, ask = TRUE, ...)
}
\arguments{
  \item{parm}{A character vector containing either
  parameter names, or a subset of \code{"all"},
  \code{"derived"}, \code{"fitted"}, and \code{"linked"};
  \code{"fitted"} corresponds to the parameters of
  interest, \code{"derived"} corresponds to quantities that
  are functions of these parameters (e.g., the effective
  survey area or animal density from an acoustic survey),
  and \code{"linked"} corresponds to the parameters AD
  Model Builder has maximised the likelihood over.}

  \item{linked}{Logical, if \code{TRUE}, intervals for
  fitted parameters are calculated on their link scales,
  then transformed back onto their "real" scales.}

  \item{method}{A character string specifying the method
  used to calculate the confidence intervals. See 'Details'
  below.}

  \item{qqplot}{Logical, if \code{TRUE} and \code{method}
  is \code{"default"} then a normal QQ plot is plotted. The
  default method is based on a normal approximation; this
  plot tests its validity.}

  \item{ask}{Logical, if \code{TRUE}, hitting return will
  show the next plot.}

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
\code{"default.bc"} is a bias-corrected version of
\code{default}, whereby the estimated bias is subtracted
from each confidence limit; \code{"basic"} for the
so-called "basic" bootstrap method; and \code{"percentile"}
for intervals calculated using the bootstrap percentile
method (the latter three are only available for objects of
class \code{admbsecr.boot}; see Davison & Hinkley, 1997,
for details).

For method \code{"default"} with objects of class
\code{admbsecr.boot}, the appropriateness of the normal
approximation can be evaluated by setting \code{qqnorm} to
\code{TRUE}. If this indicates a poor fit, set
\code{linked} to \code{TRUE} and evaluate the QQ plot to
see if this yields an improvement (see Davison & Hinkley,
1997, pp. 194, for details).
}
\references{
Davison, A. C., and Hinkley, D. V. (1997) \emph{Bootstrap
methods and their application}. Cambridge: Cambridge
University Press.
}

