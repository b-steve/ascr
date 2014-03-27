\name{vcov.admbsecr.boot}
\alias{vcov.admbsecr.boot}
\title{Extract the variance-covariance matrix from a bootstrapped admbsecr
model object}
\usage{
\method{vcov}{admbsecr.boot}(object, pars = "fitted", ...)
}
\arguments{
  \item{object}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{pars}{A character containing a subset of
  \code{"all"}, \code{"derived"}, \code{"fitted"}, and
  \code{"linked"}; \code{"fitted"} corresponds to the
  parameters of interest, \code{"derived"} corresponds to
  quantities that are functions of these parameters (e.g.,
  the effective survey area or animal density from an
  acoustic survey), and \code{"linked"} corresponds to the
  parameters AD Model Builder has maximised the likelihood
  over.}

  \item{...}{Other parameters (for S3 generic
  compatibility).}
}
\description{
Extracts the variance-covariance matrix for parameters in a
model fitted using \link[admbsecr]{admbsecr}, with a
bootstrap procedure carried out using
\link[admbsecr]{boot.admbsecr}.
}

