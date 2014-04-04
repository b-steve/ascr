\name{sim.capt}
\alias{sim.capt}
\title{Simulating SECR data}
\usage{
sim.capt(fit = NULL, traps = NULL, mask = NULL,
  infotypes = character(0), detfn = "hn", pars = NULL,
  ss.link = "identity", cutoff = NULL, call.freqs = NULL,
  freq.dist = "edf", sound.speed = 330, test.detfn = FALSE)
}
\arguments{
  \item{fit}{A fitted \code{admbsecr} model object which
  provides the additional information types, detection
  function, and parameter values from which to generate
  capture histories.}

  \item{infotypes}{A character vector indicating the
  type(s) of additional information to be simulated.
  Elements can be a subset of \code{"bearing"},
  \code{"dist"}, \code{"ss"}, \code{"toa"}, and
  \code{"mrds"} (NOTE: \code{"mrds"} not yet implemented).}

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

  \item{freq.dist}{A character string, either \code{"edf"}
  or \code{"norm"}, which specifies how the distribution
  function of the call frequencies should be estimated. If
  \code{"edf"}, then the distribution of call frequencies
  is estimated using the empirical distribution function.
  If \code{"norm"}, then a normal distribution is fitted to
  the call frequencies using the sample mean and variance.}

  \item{test.detfn}{Logical value, if \code{TRUE}, tests
  detection function to aid debugging.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap (or
  detector).}

  \item{mask}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a mask point.
  The function \link{create.mask} will return a suitable
  object.}

  \item{ss.link}{A character string, either
  \code{"indentity"} or \code{"log"}, which specifies the
  link function for the signal strength detection function.
  Only required when \code{detfn} is \code{"ss"} (i.e.,
  when there is signal strength information in
  \code{capt}).}

  \item{cutoff}{The signal strength threshold, above which
  sounds are identified as detections. Only required when
  \code{detfn} is \code{"ss"}.}

  \item{call.freqs}{A vector of call frequencies collected
  independently to an acoustic survey.}

  \item{sound.speed}{The speed of sound in metres per
  second, defaults to 330 (the speed of sound in air). Only
  used when \code{"toa"} is a component name of
  \code{capt}. \strong{Not yet implemented}.}
}
\value{
A list with named components, each corresponding to a data
type listed in \code{infotypes}. Each component is a matrix
where each row corresponds to each detected individual, and
each column corresponds to a trap (or detector). The
elements in the matrix indicate detection, and provide
simulated values of the additional information requested.
This object can be used as the \code{capt} argument for the
function \link{admbsecr}.
}
\description{
Simulates SECR capture histories and associated additional
information in the correct format for use with the function
\link{admbsecr}. If \code{fit} is provided then no other
arguments are required. Otherwise, at least \code{traps},
\code{mask}, and \code{pars} are needed.
}
\details{
See documentation for the function \link{admbsecr} for
information on the parameters corresponding to the
different detection functions, and to different types of
additional information.
}
\examples{
## Simulating based on model fits.
simple.capt <- sim.capt(simple.hn.fit)
bearing.capt <- sim.capt(bearing.hn.fit)
## Simulating from provided parameter values.
new.capt <- sim.capt(traps = example.traps, mask = example.mask, infotypes = c("bearing", "dist"), detfn = "hr",
                     pars = list(D = 2500, g0 = 0.9, sigma = 3, z = 2, kappa = 50, alpha = 10))
}

