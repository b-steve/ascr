\name{locations}
\alias{locations}
\title{Plotting estimated locations}
\usage{
locations(fit, id, infotypes = "combined", xlim = range(mask[, 1]),
  ylim = range(mask[, 2]), mask = fit$mask, cols = list(combined =
  "black", capt = "purple", ang = "green", dist = "brown", toa = "blue"),
  plot.arrows = any(c("ang", "all") \%in\% infotypes),
  plot.circles = any(c("dist", "all") \%in\% infotypes), legend = TRUE,
  add = FALSE)
}
\arguments{
  \item{fit}{A fitted model from
  \link[admbsecr]{admbsecr}.}

  \item{id}{A numeric vector with row numbers from
  \code{fit$capt}, indicating which individuals' locations
  are to be plotted.}

  \item{infotypes}{A character vector indicating the
  type(s) of information to be used when plotting the
  estimated density of location.  Elements can be a subset
  of \code{"capt"}, \code{"ang"}, \code{"dist"},
  \code{"ss"}, \code{"toa"}, \code{"combined"}, and
  \code{"all"}, where \code{"capt"} shows estimated
  location only using detection locations,
  \code{"combined"} combines all information types
  together, and \code{"all"} plots all possible contour
  types. When signal strength information is used in the
  model fit, \code{"capt"} and \code{"ss"} are equivalent
  as the signal strength information is built into the
  detection function.}

  \item{xlim}{A numeric vector of length 2, giving the x
  coordinate range.}

  \item{ylim}{A numeric vector of length 2, giving the y
  coordinate range.}

  \item{mask}{A matrix with two columns. Each row provides
  Cartesian coordinates for the location of a mask point.
  The function \link[admbsecr]{create.mask} will return a
  suitable object. The mask used to fit the model
  \code{fit} will be used by default; this argument is
  usually used when estimated location contours need to be
  plotted to a higher resolution than this.}

  \item{cols}{A list with named components corresponding to
  each contour type (i.e., a subset of \code{"capt"},
  \code{"ang"}, \code{"dist"}, \code{"toa"}, and
  \code{"combined"}). Each component provides the colour
  the associated contour type (e.g., using a character
  string such as \code{"red"}, or a call to the function
  \link[grDevices]{rgb}). By default, if only one contour
  is to be plotted, it will be plotted in black.}

  \item{plot.arrows}{Logical, if \code{TRUE} arrows
  indicating the estimated bearing to the individual are
  plotted from detectors at which detections were made.}

  \item{plot.circles}{Logical, if \code{TRUE} circles
  indicating the estimated distance to the individual are
  plotted around detectors at which detections were made.}

  \item{legend}{Logical, if \code{TRUE}, a legend will be
  added to the plot.}

  \item{add}{Logical, if \code{TRUE} contours will be added
  to an existing plot.}
}
\description{
Plots estimated densities of animal locations, which are
latent variables in SECR models.
}

