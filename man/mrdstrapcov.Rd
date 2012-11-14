\name{mrdstrapcov}
\alias{mrdstrapcov}
\title{Fitting MRDS models.}
\usage{
  mrdstrapcov(capt, mask, traps, sv, admb.dir, clean,
    verbose, trace)
}
\arguments{
  \item{capt}{capture history array.}

  \item{mask}{mask point locations.}

  \item{traps}{trap locations.}

  \item{sv}{start values.}

  \item{admb.dir}{directory containing
  mrdstrapcovsecr.tpl.}

  \item{clean}{logical, if \code{TRUE} ADMB files are
  cleaned after fitting of the model.}

  \item{verbose}{logical, if \code{TRUE} ADMB details,
  along with error messages, are printed to the R session.}

  \item{trace}{logical, if \code{TRUE} parameter values at
  each step of the fitting algorithm are printed to the R
  session.}
}
\description{
  Fits a mark-recapture distance sampling (MRDS) model,
  with different detection function parameters for each
  trap. This is a temporary function; eventually
  \code{admbsecr()} will be flexible enough to do this.
}

