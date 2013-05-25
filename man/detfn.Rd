\name{detfn}
\alias{detfn}
\alias{detfn.hn}
\alias{detfn.hr}
\alias{detfn.ss}
\alias{detfn.sstoa}
\alias{detfn.th}
\title{Estimated detection probability from a fitted model.}
\usage{
  detfn(fit, ...)

  \method{detfn}{hn} (fit, d, ...)

  \method{detfn}{hr} (fit, d, ...)

  \method{detfn}{ss} (fit, d, ...)

  \method{detfn}{sstoa} (fit, d, ...)

  \method{detfn}{th} (fit, d, ...)
}
\arguments{
  \item{fit}{a fitted model returned by
  \code{\link[admbsecr]{admbsecr}}.}

  \item{...}{additional arguments.}

  \item{d}{vector of distances from which probabilities are
  calculated.}
}
\description{
  Calculates the probability of detection for given
  distances from a fit returned by \code{admbsecr()}.
}

