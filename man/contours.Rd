\name{contours}
\alias{contours}
\alias{contours.ang}
\alias{contours.dist}
\alias{contours.disttc}
\alias{contours.mrds}
\alias{contours.simple}
\alias{contours.ss}
\alias{contours.sstoa}
\alias{contours.toa}
\title{Plotting estimated animal location densities.}
\usage{
  contours(fit, ...)

  \method{contours}{simple} (fit, dets = "all",
    add = FALSE, heat = FALSE, cols = "black", ltys = 1,
    trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{toa} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4), rgb(0, 0, 1, 0.4)),
    ltys = 1, trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{ss} (fit, dets = "all", add = FALSE,
    heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4), rgb(0, 0, 1, 0.4)),
    ltys = 1, trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{sstoa} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4), rgb(0, 0, 1, 0.4)),
    ltys = 1, trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{ang} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4), rgb(0, 0, 1, 0.4)),
    ltys = 1, trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    arrows = showcapt, xlim = NULL, ylim = NULL, ...)

  \method{contours}{disttc} (fit, dets = "all",
    add = FALSE, partition = FALSE, heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4)), ltys = 1,
    trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{dist} (fit, dets = "all", add = FALSE,
    partition = FALSE, heat = FALSE,
    cols = c("black", rgb(0, 1, 0, 0.4)), ltys = 1,
    trapnos = FALSE, problevels = NULL,
    showcapt = length(dets) == 1 && dets != "all" && !add,
    xlim = NULL, ylim = NULL, ...)

  \method{contours}{mrds} (fit, dets = "all", add = FALSE,
    trapnos = FALSE,
    showcapt = length(dets) == 1 && dets != "all" && !add,
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

  \item{cols}{a vector specifying the colours of the main,
  extra, and simple contours, in that order. Can be of
  length 1 if all three are to be the same.}

  \item{ltys}{vector specifying the line types for the
  main, simple, and extra contours, in that order. Can be
  of length 1 if all three are to be the same.}

  \item{trapnos}{logical, if \code{TRUE} the trap
  identification numbers are plotted.}

  \item{problevels}{a vector indicating which probabilities
  should be associated with the levels of the contours.}

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

  \item{arrows}{logical, indicating whether or not arrows
  are to be plotted to show estimated animal bearing from
  traps.}
}
\description{
  Plots a density of individual animals' locations from a
  fit returned by \code{admbsecr()}.
}

