\name{se.correct}
\alias{se.correct}
\title{SECR standard error correction}
\usage{
  se.correct(fit, size, calls = NULL)
}
\arguments{
  \item{fit}{an \code{admbsecr} model fit.}

  \item{size}{number of bootstrap samples.}

  \item{calls}{a vector containing the numbers of calls
  emitted from recorded individuals. If \code{NULL}, this
  defaults to what is provided by \code{fit}.}
}
\description{
  Corrects standard errors from an \code{admbsecr} fit.
}

