\name{calc.detfn}
\alias{calc.detfn}
\alias{calc.detfn.hn}
\alias{calc.detfn.hr}
\alias{calc.detfn.ss}
\alias{calc.detfn.sstoa}
\alias{calc.detfn.th}
\title{Estimated detection probability from a fitted model.}
\usage{
  calc.detfn(fit, ...)

  \method{calc.detfn}{hn} (fit, d, ...)

  \method{calc.detfn}{hr} (fit, d, ...)

  \method{calc.detfn}{ss} (fit, d, prob = TRUE, ...)

  \method{calc.detfn}{sstoa} (fit, d, prob = TRUE, ...)

  \method{calc.detfn}{th} (fit, d, ...)
}
\arguments{
  \item{fit}{a fitted model returned by
  \code{\link[admbsecr]{admbsecr}}.}

  \item{...}{additional arguments.}

  \item{d}{vector of distances from which probabilities are
  calculated.}

  \item{prob}{logical, if \code{TRUE}, capture probability
  is returned. If \code{FALSE}, expected signal strengths
  are returned.}
}
\description{
  Calculates the probability of detection for given
  distances from a fit returned by \code{admbsecr()}.
}

