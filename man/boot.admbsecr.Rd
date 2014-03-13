\name{boot.admbsecr}
\alias{boot.admbsecr}
\title{Bootstrapping SECR data}
\usage{
boot.admbsecr(fit, N, prog = TRUE)
}
\arguments{
  \item{fit}{A fitted \code{admbsecr} model object.}

  \item{N}{The number of bootstrap resamples.}

  \item{prog}{Logical, if \code{TRUE}, a progress bar is
  shown.}
}
\value{
A list of class \code{"admbsecr.boot"}. Components contain
information such as estimated parameters and standard
errors. The best way to access such information, however,
is through the variety of helper functions provided by the
admbsecr package. S3 methods \link{stdEr.admbsecr.boot} and
\link{vcov.admbsecr.boot} can be used to return values of
based on the bootstrap procedure.
}
\description{
Carries out a parametric bootstrap, based on a model fitted
using \link{admbsecr}.
}
\details{
For each bootstrap resample, a new population of
individuals is simulated within the mask area. Detections
of these individuals are simulated using the estimated
detection function. For detected individuals, additional
informatin is simulated from the estimated distribution of
measurement error. The original model is then re-fitted to
these simulated data, and parameter estimates for each
iteration saved in the component \code{boot} of the
returned list.

For fits based on acoustic surveys where the argument
\code{call.freqs} is provided to the \code{admbsecr}
function, the simulated data allocates multiple calls to
the same location based on an estimated distribution of the
call frequencies. Using a parametric bootstrap is currently
the only way parameter uncertainty can be estimated for
such models.
}
\examples{
\dontrun{
## In practice, N should be >> 100, but this leads to long computation time.
boot.fit <- boot.admbsecr(fit = simple.hn.fit, N = 100)
}
}

