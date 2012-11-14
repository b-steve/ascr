\name{admbsecr}
\alias{admbsecr}
\title{Fitting SECR models in ADMB}
\usage{
  admbsecr(capt, traps = NULL, mask, sv = "auto",
    bounds = NULL, ssqtoa = NULL, cutoff = NULL,
    admbwd = NULL, method = "simple", memory = NULL,
    profpars = NULL, clean = TRUE, verbose = FALSE,
    trace = FALSE, autogen = TRUE)
}
\arguments{
  \item{capt}{an array of dimension \code{(n, S, K)}, where
  \code{n} is the number of detected animals, \code{S} is
  number of individual sampling sessions, and \code{K} is
  the number of deployed traps. The object returned by
  \code{make.capthist()} is suitable if \code{method} is
  \code{"simple"}. Otherwise, the \code{1} values in this
  array must be changed to the value of the recorded
  supplementary information, which will depend on
  \code{method} (see 'Details'). When \code{method} is
  \code{"sstoa"} or \code{"mrds"}, this array must be of
  dimension \code{(n, S, K, 2)}. With \code{"sstoa"},
  \code{capt[, , , 1]} provides the signal strength
  information and \code{capt[, , , 2]} provides the time of
  arrival information. With \code{"mrds"}, \code{capt[, , ,
  1]} provides the binary capture history array and
  \code{capt[, , , 2]} provides the distances between all
  traps (regardless of capture) and detected animals.}

  \item{traps}{a matrix containing the coordinates of trap
  locations. The object returned by \code{read.traps()} is
  suitable.}

  \item{mask}{a mask object. The object returned by
  \code{make.mask()} is suitable.}

  \item{sv}{either \code{"auto"}, or a named vector. If
  \code{auto}, starting values for all parameters are
  automatically generated. If a vector, the elements are
  the starting values and the names indicate which
  parameter these correspond to. Starting values for all
  parameters need not be provided; they are automatically
  generated for any parameters that do not have a starting
  value explicitly provided. See 'Details' for list of
  parameters used by each method.}

  \item{bounds}{a list with optional components
  corresponding to parameters which are to have their
  default bounds overridden. Each component should be a
  vector of length two specifying the bounds, and the name
  of the component should be the name of the parameter to
  which these bounds apply. To remove default bounds from a
  parameter the component should be \code{NULL} rather than
  a vector of length two. Bounds for all parameters need
  not be provided; if there is no component corresponding
  to a parameter it keeps its default bounds. See 'Details'
  for list of parameters used by each method.}

  \item{ssqtoa}{an optional matrix. If calculated before
  call to \code{admbsecr}, providing this will prevent
  recalculation.}

  \item{cutoff}{The signal strength threshold of detection.
  Required if \code{method} is \code{"ss"} or
  \code{"sstoa"}.}

  \item{admbwd}{file path to the ADMB working directory.
  Only required if \code{autogen} is \code{TRUE}, in which
  case it points to the directory in which the \code{.tpl}
  file is located.}

  \item{method}{either \code{"simple"}, \code{"toa"},
  \code{"ang"}, \code{"ss"}, or \code{"sstoa"}. See
  'Details'.}

  \item{memory}{value of \code{arrmblsize} in ADMB.
  Increase this if ADMB reports a memory error.}

  \item{profpars}{character vector of names of parameters
  over which profile likelihood should occur.}

  \item{clean}{logical, if \code{TRUE} ADMB files are
  cleaned after fitting of the model.}

  \item{verbose}{logical, if \code{TRUE} ADMB details,
  along with error messages, are printed to the R session.}

  \item{trace}{logical, if \code{TRUE} parameter values at
  each step of the fitting algorithm are printed to the R
  session.}

  \item{autogen}{logical, if \code{TRUE}, the appropriate
  \code{.tpl} file is written to \code{admbwd} (or the
  current working directory if \code{admbwd} is
  \code{NULL}). If \code{FALSE}, the \code{.tpl} file
  should already be located in \code{admbwd} (or the
  current working directory if \code{admb} is \code{NULL}).
  Usually only set to \code{FALSE} for development
  purposes.}
}
\value{
  An object of class 'admb'.

  The following functions can be used to extract model
  components: \code{summary()}, \code{AIC()},
  \code{logLik()}, \code{deviance()}, \code{vcov()},
  \code{coef()}, \code{stdEr()}, and \code{confint()}.

  The latter takes arguments \code{level} and
  \code{method}, which specify the confidence level and
  calculation method respectively. The default method gives
  quadratic (Wald) intervals based on approximate standard
  errors; \code{"profile"} gives profile likelihood
  intervals, and can be used if the \code{admbsecr()}
  parameter \code{profpars} is non-null and provides names
  of model parameters that are to be profiled.
}
\description{
  Fits an SECR model, with our without supplementary
  information relevant to animal location. Parameter
  estimation is done by maximum likelihood in ADMB.
}
\details{
  ADMB is called to fit an SECR model through use of the
  R2admb package. Different methods are used depending on
  the additional information on animal location that is
  collected. These are: \itemize{ \item \code{"simple"}:
  Normal SECR with no additional information. Parameters to
  estimate are: \itemize{ \item D: Animal density.

  \item g0: Probability of detection at distance 0.

  \item sigma: 'Standard deviation' parameter for
  halfnormal detection function.  } \item \code{"toa"}:
  SECR with precise time of arrival (TOA) recorded.
  Parameters to estimate are: \itemize{ \item D: As above.

  \item g0: As above.

  \item sigma: As above.

  \item sigmatoa: Error term associated with TOA.  } \item
  \code{"ang"}: SECR with estimates of angle to animal
  recorded. Parameters to estimate are: \itemize{ \item D:
  As above.

  \item g0: As above.

  \item sigma: As above.

  \item kappa: Error term associated with angle estimation.
  } \item \code{"ss"}: SECR with received signal strengths
  at traps recorded. Parameters to estimate are: \itemize{
  \item D: As above.

  \item ssb0: Average signal strength at sound source.

  \item ssb1: Decrease in signal strength per unit distance
  due to sound propagation.

  \item sigmass: Error term associated with signal
  strength.  } \item \code{"sstoa"}: SECR with precise TOA
  and received signal strengths at traps recorded.
  Parameters to estimate are: \itemize{ \item D: As above.

  \item sigmatoa: As above.

  \item ssb0: As above.

  \item ssb1: As above.

  \item sigmass: As above.  } \item \code{"dist"}: SECR
  with estimated distances between animals and traps at
  which detections occurred. Parameters to estimate are:
  \itemize{ \item D: As above.

  \item g0: As above.

  \item sigma: As above.

  \item alpha: Shape parameter associated with the gamma
  distribution used to model estimated distances.  } }
}
\author{
  Ben Stevenson
}
\seealso{
  \code{\link[R2admb]{do_admb}},
  \code{\link[secr]{secr.fit}},
  \code{\link[secr]{make.capthist}},
  \code{\link[secr]{read.traps}}.
}

