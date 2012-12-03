\name{contours}
\alias{contours}
\alias{contours.ang}
\alias{contours.simple}
\alias{contours.ss}
\alias{contours.toa}
\title{Plotting estimated animal location densities.}
\usage{
  contours(fit, ...)

  \method{contours}{simple} (fit, which = "all",
    add = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE, ...)

  \method{contours}{toa} (fit, which = "all", add = FALSE,
    partition = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE, ...)

  \method{contours}{ss} (fit, which = "all", add = FALSE,
    heat = FALSE, col = "black", trapnos = FALSE, ...)

  \method{contours}{ang} (fit, which = "all", add = FALSE,
    partition = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE, ...)
}
\arguments{
  \item{fit}{a fitted model returned by \code{admbsecr()}.}

  \item{...}{additional arguments to be passed to
  \code{\link[graphics]{contour}}.}

  \item{which}{which individuals' location densities are
  plotted.}

  \item{add}{logical, if \code{TRUE} the contours are added
  to an already existing plot.}

  \item{heat}{logical, if \code{TRUE} a levelplot is used
  instead of contours.}

  \item{col}{specifies the colour of the contours to be
  plotted.}

  \item{trapnos}{logical, if \code{TRUE} the trap
  identification numbers are plotted.}

  \item{partition}{logical, if \code{TRUE} the
  contributions to the countour due to both the binary
  capture history data and the supplementary information
  are also plotted.}
}
\description{
  Plots a density of individual animals' locations from a
  fit returned by \code{admbsecr()}.
}

