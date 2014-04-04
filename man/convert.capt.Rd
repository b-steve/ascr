\name{convert.capt}
\alias{convert.capt}
\title{Convert capture history object}
\usage{
convert.capt(capt, traps, capthist = TRUE)
}
\arguments{
  \item{capthist}{Logical, if \code{TRUE}, a
  \code{capthist} object is returned. Otherwise a data
  frame is returned, which is suitable for the
  \code{captures} argument to the \link{make.capthist}
  function.}

  \item{capt}{A list with named components, containing the
  capture history and supplementary information. The
  function \link{create.capt} will return a suitable
  object. See 'Details' below.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap (or
  detector).}
}
\value{
An object of class \link{capthist}.
}
\description{
Converts an \code{admbsecr} capture history list to a
\code{secr} capthist object.
}
\details{
The returned object is suitable for use as the
\code{capthist} argument of the function \link{secr.fit}.
}
\examples{
capt <- convert.capt(capt = example.capt, traps = example.traps)
}

