\name{admbsecr}
\alias{admbsecr}
\title{Fitting SECR models in ADMB}
\usage{
admbsecr(capt, traps, mask, detfn = "hn", sv = NULL, bounds = NULL,
  fix = NULL, sf = NULL, ss.link = "identity", cutoff = NULL,
  call.freqs = NULL, sound.speed = 330, trace = FALSE, clean = TRUE,
  cbs = NULL, gbs = NULL, exe.type = "old")
}
\arguments{
  \item{capt}{A list with named components, containing the
  capture history and supplementary information. The
  function \link{create.capt} will return a suitable
  object. See 'Details' below.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap (or
  detector).}

  \item{mask}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a mask point.
  The function \link{create.mask} will return a suitable
  object.}

  \item{detfn}{A character string specifying the detection
  function to be used. Options are "hn" (halfnormal), "hr"
  (hazard rate), "th" (threshold), "lth" (log-link
  threshold), or "ss" (signal strength). If the latter is
  used, signal strength information must be provided in
  \code{capt}.}

  \item{sv}{A named list. Component names are parameter
  names, and each component is a start value for the
  associated parameter. See 'Details' for further
  information on the parameters to be fitted.}

  \item{bounds}{A named list. Component names are parameter
  names, and each components is a vector of length two,
  specifying the bounds for the associated parameter.}

  \item{fix}{A named list. Component names are parameter
  names to be fixed, and each component is the fixed value
  for the associated parameter.}

  \item{sf}{A named list. Component names are parameter
  names, and each component is a scalefactor for the
  associated parameter. The default behaviour is to
  automatically select scalefactors based on parameter
  start values. See the section on convergence below.}

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

  \item{trace}{Logical, if \code{TRUE} parameter values at
  each step of the optimisation algorithm are printed to
  the R console.}

  \item{clean}{Logical, if \code{TRUE} ADMB output files
  are removed. Otherwise, ADMB output file will remain in a
  directory, the location of which is reported after the
  model is fitted.}

  \item{cbs}{The CMPDIF_BUFFER_SIZE, set using the
  \code{-cbs} option of the executable created by ADMB.
  This can be increased to speed up optimisation if
  \code{cmpdiff.tmp} gets too large (please ignore, unless
  you are familiar with ADMB and know what you are doing).}

  \item{gbs}{The GRADSTACK_BUFFER_SIZE, set using the
  \code{-gbs} option of the executable created by ADMB.
  This can be increased to speed up optimisation if
  \code{gradfil1.tmp} gets too large (please ignore, unless
  you are familiar with ADMB and know what you are doing).}

  \item{exe.type}{Character string, either \code{"old"} or
  \code{"new"}, depending on which executable is to be used
  (for development purposes only; please ignore).}
}
\value{
A list of class \code{"admbsecr"}. Components contain
information such as estimated parameters and standard
errors. The best way to access such information, however,
is through the variety of helper functions provided by the
admbsecr package.
}
\description{
Fits an SECR model, with our without supplementary
information relevant to animal location. Parameter
estimation is done by maximum likelihood through an AD
Model Builder (ADMB) executable.
}
\details{
ADMB uses a quasi-Newton method to find maximum likelihood
estimates for the model parameters. Standard errors are
calculated by taking the inverse of the negative of the
Hessian. Alternatively, \link{boot.admbsecr} can be used to
carry out a parametric bootstrap procedure, from which
parameter uncertainty can also be inferred.

The class of model fitted by this function (and, indeed,
around which this package is based) was first proposed by
Borchers et al. (in press); this reference is a good
starting point for practitioners looking to implement these
methods.

If the data are from an acoustic survey where individuals
call more than once (i.e., the argument \code{call.freqs}
contains values that are not 1), then standard errors
calculated from the inverse of the negative Hessian are not
correct. They are therefore not provided in this case. The
method used by the function \link{boot.admbsecr} is
currently the only way to calculate these reliably (see
Stevenson et al., in prep., for details).
}
\section{The \code{capt} argument}{
  The \code{capt} argument is a list with named components.
  Each component must be an \eqn{n} by \eqn{k} matrix,
  where \eqn{n} is the number of detections made, and
  \eqn{k} is the number of traps (or detectors) deployed. A
  component named \code{bincapt} is compulsory.

  Further optional component names each refer to a type of
  information which is informative on animal location
  collected on each detection. Possible choices are:
  \code{bearing}, \code{dist}, \code{ss}, \code{toa}, and
  \code{mrds}.

  If the \eqn{i}th individual evaded the \eqn{j}th trap (or
  detector), then the \eqn{j}th element in the \eqn{i}th
  row should be 0 for all components. Otherwise, if the
  \eqn{i}th individual was trapped (or detected) by the
  \eqn{j}th trap (or detector), then: \itemize{ \item For
  the \code{bincapt} component, the element should be 1.
  \item For the \code{bearing} component, the element
  should be the estimated bearing from which the detector
  detected the individual.  \item For the \code{dist}
  component, the element should be the estimated distance
  between the individual and the detector at the time of
  the detection.  \item For the \code{ss} component, the
  element should be the measured signal strength of an
  acoustic signal detected by the detector (only possible
  when the detectors are microphones).  \item For the
  \code{toa} component, the element should be the measured
  time of arrival (in seconds) since the start of the
  survey (or some other reference time) of an acoustic
  signal detected by the detector (only possible when the
  detectors are microphones).  \item For the \code{mrds}
  component, the element should be the \emph{known} (not
  estimated) distance between the individual and the
  detector at the time of the detection. }
}

\section{Fitted parameters}{
  The parameter \code{D}, the density of individuals (or,
  in an acoustic survey, the density of calls) is always
  fitted. The effective survey area, \code{esa}, (see
  Borchers, 2012, for details) is always provided as a
  derived parameter, with a standard error calculated using
  the delta method.

  Further parameters to be fitted depend on the choice of
  the detection function (i.e., the \code{detfn} argument),
  and the types of additional information collected (i.e.,
  the components in the \code{capt}).

  Details of the detection functions are as follows:

  For \code{detfn = "hn"}: \itemize{ \item Estimated
  paramters are \code{g0} and \code{sigma}.  \item
  \eqn{g(d) = g_0\ exp(-d^2/(2\sigma^2))}{g(d) = g0 * exp(
  -d^2 / (2 * sigma^2 ))} }

  For \code{detfn = "hr"}: \itemize{ \item Estimated
  parameters are \code{g0}, \code{sigma}, and \code{z}.
  \item \eqn{g(d) = g_0\ (1 - exp(-(d/\sigma)^{-z}))}{g(d)
  = g0 * ( 1 - exp( -(d/sigma)^{-z} ) )} }

  For \code{detfn = "lth"}: \itemize{ \item Estimated
  parameters are \code{shape.1}
  \ifelse{latex}{(\eqn{\kappa})}{}, \code{shape.2}
  \ifelse{latex}{(\eqn{\nu})}{}, and \code{scale}
  \ifelse{latex}{(\eqn{\tau})}{}.  \item \eqn{g(d) = 0.5 -
  0.5\ erf(\kappa - exp(\nu - \tau d))}{g(d) = 0.5 - 0.5 *
  erf( shape.1 - exp( shape.2 - scale * d ) )} }

  For \code{detfn = "th"}: \itemize{ \item Estimated
  parameters are \code{shape}
  \ifelse{latex}{(\eqn{\kappa})}{} and \code{scale}
  \ifelse{latex}{(\eqn{\tau})}{}.  \item \eqn{g(d) = 0.5 -
  0.5\ erf(d/\tau - \kappa)}{g(d) = 0.5 - 0.5 * erf(
  d/scale - shape )} }

  For \code{detfn = "ss"}: \itemize{ \item The signal
  strength detection function is special in that it
  requires signal strength information to be collected in
  order for all parameters to be estimated.  \item
  Estimated parameters are \code{b0.ss}, \code{b1.ss}, and
  \code{sigma.ss}.  \item The expected signal strength is
  modelled as: \eqn{E(SS) = h^{-1}(\beta_0 -
  \beta_1d)}{E(SS) = h^{-1}(b0.ss - b1.ss*d)}, where
  \eqn{h} is specified by the argument \code{ss.link}. }

  Details of the parameters associated with different
  additional data types are as follows:

  For data type \code{"bearing"}, \code{kappa} is
  estimated. This is the concerntration parameter of the
  von-Mises distribution used for measurement error in
  estimated bearings.

  For data type \code{"dist"}, \code{alpha} is estimated.
  This is the shape parameter of the gamma distribution
  used for measurement error in estimated distances.

  For data type \code{"toa"}, \code{sigma.toa} is
  estimated. This is the standard deviation parameter of
  the normal distribution used for measurement error in
  recorded times of arrival.

  For data type \code{"mrds"}, no extra parameters are
  estimated. Animal location is assumed to be known.
}

\section{Convergence}{
  If maximum likelihood estimates could not be found during
  optimisation, then \code{admbsecr} will usually show a
  warning that the maximum gradient component is large, or
  possibly throw an error reporting that a \code{.par} file
  is missing.

  The best approach to fixing convergence issues is to
  re-run the \code{admbsecr} function with the argument
  \code{trace} set to \code{TRUE}. Parameter values will be
  printed out for each step of the optimisation algorithm.

  First, look for a large jump in a parameter to a value
  far from what is feasible. This issue can be fixed by
  using the \code{bounds} argument to restrict the
  parameter space over which ADMB searches for the maximum
  likelihood estimate.

  Alternatively, try a different set of start values using
  the argument \code{sv}; by default \code{admbsecr} will
  choose some start values, but these are not necessarily
  sensible. The start values that were used appear as the
  first line of text when \code{trace} is \code{TRUE}.

  Sometimes the algorithm appears to converge, but
  nevertheless perseveres reporting the same parameter
  values again and again for a while (prior to the
  calculation of the Hessian). This is because ADMB has
  failed to detect convergence as at least one of the
  gradient components is still larger than the convergence
  criterion (by default, 0.0001). It is possible to speed
  things up and help ADMB detect convergence earlier by
  either tightening parameter bounds (as above), or by
  setting appropriate scalefactors (using the argument
  \code{sf}). To do this, first identify which parameters
  have large gradient components from the "final
  statistics" section of the \code{trace} output. Next,
  find the default settings of the scalefactors by printing
  the object \code{fit$args$sf}, where \code{fit} is the
  original object returned by \code{admbsecr}. Finally,
  rerun \code{admbsecr} again, but this time set the
  argument \code{sf} manually. Set scalefactors for any
  parameters with small gradient components to the same as
  the defaults ascertained above, and increase those
  associated with large gradient components by a factor of
  10. If the problem persists, repeat this process (e.g.,
  if the same parameters still have large gradient
  components, increase the associated scalefactors by
  another factor of 10).
}
\examples{
\dontrun{
simple.capt <- example.capt["bincapt"]
simple.hn.fit <- admbsecr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, fix = list(g0 = 1))
simple.hr.fit <- admbsecr(capt = simple.capt, traps = example.traps,
                          mask = example.mask, detfn = "hr")
bearing.capt <- example.capt[c("bincapt", "bearing")]
bearing.hn.fit <- admbsecr(capt = bearing.capt, traps = example.traps,
                           mask = example.mask, fix = list(g0 = 1))
}
}
\references{
Borchers, D. L. (2012) A non-technical overview of
spatially explicit capture-recapture models. \emph{Journal
of Ornithology}, \strong{152}: 435--444.

Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas, L.,
and Marques, T. A. (in press) A unifying model for
capture-recapture and distance sampling surveys of wildlife
populations. \emph{Journal of the American Statistical
Association}.

Stevenson, B. C., Borchers, D. L., Altwegg, R., Measey, G.
J., Swift, R. J., and Gillespie, D. M. (in prep.) An
acoustic spatially explicit capture-recapture method for
estimating vocalizing amphibian density.
}
\seealso{
\link{boot.admbsecr} to calculate standard errors and
estimate bias using a parametric bootstrap.

\link{coef.admbsecr}, \link{stdEr.admbsecr}, and
\link{vcov.admbsecr} to extract estimated parameters,
standard errors, and the variance-covariance matrix,
respectively.

\link{confint.admbsecr} to calculate confidence intervals.

\link{summary.admbsecr} to get a summary of estimates and
standard errors.

\link{show.detfn} to plot the estimated detection function.

\link{locations} to plot estimated locations of particular
individuals or calls.
}

