\name{contours}
\alias{contours}
\title{Plotting estimated animal location densities.}
\usage{
  contours(fit, which = "all", add = FALSE,
    partition = FALSE, col = "black", trapnos = FALSE)
}
\arguments{
  \item{fit}{a fitted model returned by \code{admbsecr()}.}

  \item{which}{which individuals' location densities are
  plotted.}

  \item{add}{logical, if \code{TRUE} the contours are added
  to an already existing plot.}

  \item{partition}{logical, if \code{TRUE} the
  contributions to the countour due to both the binary
  capture history data and the supplementary information
  are also plotted.}

  \item{col}{specifies the colour of the contours to be
  plotted.}

  \item{trapnos}{logical, if \code{TRUE} the trap
  identification numbers are plotted.}
}
\description{
  Plots a density of individual animals' locations from a
  fit returned by \code{admbsecr()}.
}

