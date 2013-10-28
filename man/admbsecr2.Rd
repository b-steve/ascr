\name{admbsecr2}
\alias{admbsecr2}
\title{Fitting SECR models in ADMB}
\usage{
  admbsecr2(capt, traps, mask, detfn = "hn", sv = NULL,
    bounds = NULL, fix = NULL, scalefactors = NULL,
    trace = FALSE)
}
\arguments{
  \item{capt}{A list with named components, containing the
  capture history and supplementary information.}

  \item{traps}{A matrix with two columns. The rows provide
  Cartesian coordinates for trap locations.}

  \item{mask}{A matrix with two columns. The rows provide
  Cartesian coordinates for the mask point locations.}

  \item{detfn}{A character string specifying the detection
  function to be used. Options are "hn" (halfnormal), "hr"
  (hazard rate), "th" (threshold), "lth" (log-link
  threshold), "ss" (signal strength), or "logss" (log-link
  signal strength). If either of the latter two are used,
  signal strength information must be provided in
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
  parameter names, and each components is a scalefactor for
  the associated parameter. The default behaviour is to
  automatically select scalefactors based on parameter
  start values.}

  \item{trace}{logical, if \code{TRUE} parameter values at
  each step of the optimisation algorithm are printed to
  the R session.}
}
\description{
  Fits an SECR model, with our without supplementary
  information relevant to animal location. Parameter
  estimation is done by maximum likelihood through and ADMB
  executible.
}

