\name{try.admbsecr}
\alias{try.admbsecr}
\title{Function for stable model fitting.}
\usage{
  try.admbsecr(sv = "auto", ...)
}
\arguments{
  \item{sv}{named vector of start values.}

  \item{...}{further arguments passed to \code{admbsecr}.}
}
\description{
  Tries to fit an admbsecr model, and reverts to automatic
  start values if it fails.
}

