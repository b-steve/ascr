\name{show.detfn}
\alias{show.detfn}
\title{Show estimated detection function}
\usage{
  show.detfn(fit = NULL, detfn = NULL, pars = NULL,
    xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
    add = FALSE, ...)
}
\arguments{
  \item{fit}{a fitted model from \code{admbsecr()}.}

  \item{detfn}{the detection function to be plotted.}

  \item{pars}{parameter values for the chosen detection
  function.}

  \item{xlim}{x-axis limits.}

  \item{ylim}{y-axis limits.}

  \item{xlab}{x-axis label.}

  \item{ylab}{y-axis label.}

  \item{add}{logical, if \code{TRUE}, add to existing
  plot.}

  \item{...}{further arguments to be passed to plotting
  functions.}
}
\description{
  Plots the estimated detection function from an
  \code{admbsecr} fit.
}

