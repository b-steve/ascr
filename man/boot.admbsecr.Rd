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
\description{
Carries out a parametric bootstrap, based on a model fitted
using \link[admbsecr]{admbsecr}. For fits based on acoustic
surveys where the argument \code{call.freqs} is provided to
the \code{admbsecr} function, the simulated data allocates
multiple calls to the same location based on an estimated
distribution of the call frequencies. Using a parametric
bootstrap is currently the only way parameter uncertainty
can be estimated for such models (see Details).
}

