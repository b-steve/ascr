\name{admbsecr}
\alias{admbsecr}
\title{Fitting SECR models in ADMB}
\usage{
admbsecr(capt, traps, mask, detfn = "hn", sv = NULL, bounds = NULL,
  fix = NULL, scalefactors = NULL, ss.link = "identity", cutoff = NULL,
  sound.speed = NULL, trace = FALSE, clean = TRUE, exe.type = "old")
}
\arguments{
  \item{capt}{A list with named components, containing the
  capture history and supplementary information.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap.}

  \item{mask}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a mask point.}

  \item{detfn}{A character string specifying the detection
  function to be used. Options are "hn" (halfnormal), "hr"
  (hazard rate), "th" (threshold), "lth" (log-link
  threshold), or "ss" (signal strength). If the latter is
  used, signal strength information must be provided in
  \code{capt}.}

  \item{sv}{A named list. Component names are parameter
  names, and each component is a start value for the
  associated parameter.}

  \item{bounds}{A named list. Component names are parameter
  names, and each components is a vector of length two,
  specifying the bounds for the associated parameter.}

  \item{fix}{A named list. Component names are parameter
  names to be fixed, and each component is the fixed value
  for the associated parameter.}

  \item{scalefactors}{A named list. Component names are
  parameter names, and each component is a scalefactor for
  the associated parameter. The default behaviour is to
  automatically select scalefactors based on parameter
  start values.}

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

  \item{trace}{Logical, if \code{TRUE} parameter values at
  each step of the optimisation algorithm are printed to
  the R session.}

  \item{clean}{Logical, if \code{TRUE} ADMB output files
  are removed.}

  \item{exe.type}{Character string, either "old" or "new",
  depending on which executable is to be used (for
  development purpouses only; please ignore).}
}
\description{
Fits an SECR model, with our without supplementary
information relevant to animal location. Parameter
estimation is done by maximum likelihood through and ADMB
executible.
}

