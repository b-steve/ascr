\name{sim.capt}
\alias{sim.capt}
\title{Simulating SECR data}
\usage{
  sim.capt(fit, calls = 1, traps = NULL, mask = NULL,
    pars = NULL, method = "simple", detfn = "hn",
    cutoff = NULL, sound.speed = NULL)
}
\arguments{
  \item{fit}{a fitted \code{admbsecr} model object which
  provides the method, detection function, and parameter
  values from which to generate capture histories.}

  \item{calls}{a vector containing call frequencies from
  monitored individuals.}

  \item{traps}{a matrix containing the coordinates of trap
  locations. The object returned by
  \code{\link[secr]{read.traps}} is suitable.}

  \item{mask}{a mask object. The object returned by
  \code{\link[secr]{make.mask}} is suitable.}

  \item{pars}{a named vector of parameter values.}

  \item{method}{the \code{admbsecr} method which specifies
  the additional information to be simulated.}

  \item{detfn}{the admbsecr detection function to be used.}

  \item{cutoff}{the signal strength threshold of detection.
  Required if \code{method} is \code{"ss"} or
  \code{"sstoa"}.}

  \item{sound.speed}{the speed of sound in metres per
  second. Used for TOA analysis.}
}
\description{
  Simulates SECR data capture histories in the correct
  format for use with \code{\link[admbsecr]{admbsecr}}.
}
\details{
  If \code{fit} is provided then no other arguments are
  required. Otherwise, at least \code{traps}, \code{mask},
  \code{pars}, \code{method}, and \code{detfn} are needed.
}

