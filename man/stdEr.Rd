\name{stdEr}
\alias{stdEr}
\title{Extract parameter standard errors.}
\usage{
  stdEr(fit, type = "fixed")
}
\arguments{
  \item{fit}{a fitted model from \code{admbsecr()}.}

  \item{type}{a character string, either \code{"fixed"}, or
  \code{"all"}. If \code{"fixed"} only model parameter
  standard errors are shown, otherwise a standard error
  (calculated using the delta method) is also provided for
  the effective sampling area.}
}
\description{
  Extracts standard errors from an admbsecr fit.
}

