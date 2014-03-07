\name{convert.mask}
\alias{convert.mask}
\title{Convert mask object}
\usage{
convert.mask(mask)
}
\arguments{
  \item{mask}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a mask point.
  The function \link[admbsecr]{create.mask} will return a
  suitable object.}
}
\description{
Converts an \code{admbsecr} mask matrix to a \code{secr}
mask object.
}

