\name{convert.traps}
\alias{convert.traps}
\title{Convert traps object}
\usage{
convert.traps(traps)
}
\arguments{
  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap (or
  detector).}
}
\value{
An object of class \code{traps} comprising a data frame of
x- and y-coordinates, the detector type ('single', 'multi',
'proximity', 'count', 'polygon' etc.), and possibly other
attributes.
}
\description{
Converts an \code{admbsecr} traps matrix to a \code{secr}
traps object.
}
\details{
The returned object is suitable for use as the \code{traps}
argument of the function \link{make.capthist}.
}
\examples{
traps <- convert.traps(traps = example.traps)
}

