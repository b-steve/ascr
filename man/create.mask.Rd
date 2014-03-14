\name{create.mask}
\alias{create.mask}
\title{Create mask object}
\usage{
create.mask(traps, buffer, ...)
}
\arguments{
  \item{buffer}{The minimum distance between trap locations
  and the edge of the generated mask.}

  \item{...}{Arguments to be passed to \link{make.mask}.}

  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap (or
  detector).}
}
\value{
An object of class \code{mask}.
}
\description{
Creates a mask object to use with the function
\code{\link{admbsecr}}.
}
\examples{
mask <- create.mask(traps = example.traps, buffer = 20)
}
\seealso{
\link{convert.mask} to convert a mask compatible with the
\link{secr} package.
}

