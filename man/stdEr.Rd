\name{stdEr}
\alias{stdEr}
\title{Extract parameter standard errors.}
\usage{
  stdEr(fit, type = "fixed")
}
\arguments{
  \item{fit}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{type}{A character string. If \code{"fixed"} then
  only standard errors for estimated parameters are
  presented. If \code{"all"} then standard errors are also
  presented for derived paramters (e.g., the effective
  survey area).}
}
\description{
  Extracts standard errors from an admbsecr fit.
}

