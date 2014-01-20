\name{sim.capt}
\alias{sim.capt}
\title{Simulating SECR data}
\usage{
  sim.capt(fit = NULL, traps = NULL, mask = NULL,
    infotypes = character(0), detfn = "hn", pars = NULL,
    ss.link = "identity", cutoff = NULL, sound.speed = 330)
}
\arguments{
  \item{fit}{A fitted \code{admbsecr} model object which
  provides the additional information types, detection
  function, and parameter values from which to generate
  capture histories.}

  \item{traps}{A matrix with two columns. The rows provide
  Cartesian coordiates for trap locations.}

  \item{mask}{A matrix with two columns. The rows provide
  Cartesian coordinates for the mask point locations.}

  \item{infotypes}{A character vector indicating the
  type(s) of additional information to be simulated.
  Elements can be a subset of \code{"ang"}, \code{"dist"},
  \code{"ss"}, \code{"toa"}, and \code{"mrds"} (NOTE:
  \code{"mrds"} not yet implemented).}

  \item{detfn}{A character string specifying the detection
  function to be used. Options are "hn" (halfnormal), "hr"
  (hazard rate), "th" (threshold), "lth" (log-link
  threshold), or "ss" (signal strength).}

  \item{pars}{A named list. Component names are parameter
  names, and each component is the value of the associated
  parameter. A value for the parameter \code{D}, animal
  density (or call density, if it an acoustic survey) must
  always be provided, along with values for parameters
  associated with the chosen detection function and
  additional information type(s).}

  \item{ss.link}{A character string, either
  \code{"indentity"} or \code{"log"}, which specifies the
  link function for the signal strength detection function.
  Only required when \code{detfn} is \code{"ss"}.}

  \item{cutoff}{The signal strength threshold, above which
  sounds are identified as detections. Only required when
  \code{detfn} is \code{"ss"}.}

  \item{sound.speed}{The speed of sound in metres per
  second, defaults to 330 (the speed of sound in air). Only
  used when \code{info} includes \code{"toa"}.}
}
\description{
  Simulates SECR capture histories and associated
  additional information in the correct format for use with
  \code{\link[admbsecr]{admbsecr}}.
}
\details{
  If \code{fit} is provided then no other arguments are
  required. Otherwise, at least \code{traps}, \code{mask},
  and \code{pars} are needed.
}

