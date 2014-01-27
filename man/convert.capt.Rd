\name{convert.capt}
\alias{convert.capt}
\title{Convert capture history object}
\usage{
convert.capt(capt, traps, capthist = TRUE)
}
\arguments{
  \item{capt}{A capture history list.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap.}

  \item{capthist}{Logical, if \code{TRUE}, a
  \code{capthist} object is returned. Otherwise a data
  frame is returned, which is suitable for the
  \code{captures} argument to the
  \code{\link[secr]{make.capthist}} function.}
}
\description{
Converts an \code{admbsecr} capture history list to a
\code{secr} capthist object.
}

