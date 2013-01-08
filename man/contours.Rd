\name{contours}
\alias{contours}
\alias{contours.ang}
\alias{contours.disttc}
\alias{contours.simple}
\alias{contours.ss}
\alias{contours.toa}
\title{Plotting estimated animal location densities.}
\usage{
  contours(fit, ...)

  \method{contours}{simple} (fit, dets = "all",
    add = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all",
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{toa} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all",
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{ss} (fit, dets = "all", add = FALSE,
    heat = FALSE, col = "black", trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all",
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{ang} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE, col = "black",
    trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all",
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{disttc} (fit, dets = "all",
    add = FALSE, partition = FALSE, heat = FALSE,
    col = "black", trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all",
    xlim = NULL, ylim = NULL, ...)
}
\arguments{
  \item{fit}{a fitted model returned by
  \code{\link[admbsecr]{admbsecr}}.}

  \item{...}{additional arguments to be passed to
  \code{\link[graphics]{contour}}.}

  \item{dets}{which individuals' location densities are
  plotted.}

  \item{add}{logical, if \code{TRUE} the contours are added
  to an already existing plot.}

  \item{heat}{logical, if \code{TRUE} a levelplot is used
  instead of contours.}

  \item{col}{specifies the colour of the contours to be
  plotted.}

  \item{trapnos}{logical, if \code{TRUE} the trap
  identification numbers are plotted.}

  \item{showcapt}{logical, if \code{TRUE} circles are drawn
  around detectors on which the detection was made.}

  \item{xlim}{numeric vector of length 2, giving the x
  coordinate range.}

  \item{ylim}{numeric vector of length 2, giving the y
  coordinate range.}

  \item{partition}{if logical, \code{TRUE} indicates that
  the contributions to the countour due to both the binary
  capture history data and the supplementary information
  are also to be plotted. If character, \code{"all"} and
  \code{"none"} correspond to \code{"TRUE"} and
  \code{"FALSE"} respectively. \code{"simple"} indicates
  only to add contours due to the binary capture history
  data, and \code{"extra"} indicates only to add contours
  due to the supplementary information.}
}
\description{
  Plots a density of individual animals' locations from a
  fit returned by \code{admbsecr()}.
}

