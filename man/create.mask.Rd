\name{create.mask}
\alias{create.mask}
\title{Create mask object}
\usage{
create.mask(traps, buffer, ...)
}
\arguments{
  \item{traps}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a trap.}

  \item{buffer}{The minimum distance between trap locations
  and the edge of the generated mask.}

  \item{...}{Arguments to be passed to
  \link[secr]{make.mask}.}
}
\description{
Creates a mask object to use with the function
\code{\link[admbsecr]{admbsecr}}.
}

