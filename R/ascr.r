## Main package-level documentation.

#' ascr: A package to fit acoustic spatial capture-recapture models
#'
#' Functions to fit acoustic spatial capture-recapture models.
#' 
#' Spatial capture-recapture (SCR) models estimate animal density from
#' records of when and were individuals were detected. The first SCR
#' methods for acoustic data were developed by Efford, Dawson, and
#' Borchers (2009).
#'
#' The \code{ascr} package fits SCR models to estimate animal density
#' from acoustic surveys. In particular, these models can incorporate
#' additional data like received signal strengths, times of arrival,
#' estimated distances, and estimated bearings, which are informative
#' about animals' locations. Other more general software packages
#' exist for SCR (e.g., the comprehensive \link{secr} R package), but
#' they cannot handle all data types listed above.
#'
#' Models that can be fitted in \code{ascr} include those described by
#' Efford, Dawson, and Borchers (2009), Borchers et al (2015), and
#' Stevenson et al (2015). The package also implements some unpublished
#' methods that are described in Stevenson (2016), a PhD thesis.
#'
#' @section Data structure:
#'
#' The easiest way to prepare your data for model fitting is via the
#' functions \link{create.capt} and \link{create.mask}. Also see the
#' vignette on data structure, which can be viewed by running
#' \code{vignette("data-structure")}.
#'
#' @section Model fitting:
#'
#' Model fitting is carried out by the \link{fit.ascr} function.
#' 
#' @section Output and plots:
#'
#' Useful information like parameter estimates, standard errors,
#' confidence intervals, and variance-covariance matrices can be
#' extracted from a fitted model object using the standard R
#' functions, like \code{coef}, \code{stdEr}, \code{summary},
#' \code{confint}, and \code{vcov} (see \link{coef.ascr},
#' \link{stdEr.ascr}, \link{summary.ascr}, \link{confint.ascr}, and
#' \link{vcov.ascr}).
#'
#' Plots of the estimated detection function and of the detection
#' probability surface can be obtained using \link{show.detfn} and
#' \link{show.detsurf}, respectively. For inhomogeneous density
#' models, \link{show.Dsurf} plots the estimated density surface. The
#' \link{locations} function plots estimated locations of individuals
#' or calls.
#'
#' @section Simulation:
#'
#' The \link{sim.capt} function can simulate SCR data, either from
#' user-specified parameters or estimated parameter values extracted
#' from a fitted model object.
#'
#' @references Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas,
#'     L., and Marques, T. A. (2015). A unifying model for
#'     capture-recapture and distance sampling surveys of wildlife
#'     populations. \emph{Journal of the American Statistical
#'     Association}, \strong{110}: 195--204.
#' @references Efford, M. G., Dawson, D. K., and Borchers,
#'     D. L. (2009). Population density estimated from locations of
#'     individuals on a passive detector array. \emph{Ecology},
#'     \strong{90}: 2676--2682.
#' @references Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift,
#'     R. J., Gillespie, D. M., and Measey, G. J. (2015). A general
#'     framework for animal density estimation from acoustic
#'     detections across a fixed microphone array. \emph{Methods in
#'     Ecology and Evolution}, \strong{6}: 38--48.
#' @references Stevenson, B. C. (2016) \emph{Methods in Spatially
#'     Explicit Capture-Recapture}. PhD thesis, University of St
#'     Andrews.
#'
#' @docType package
#' @name ascr
NULL

## Roxygen code for NAMESPACE and datasets.

## Package imports for roxygenise to pass to NAMESPACE.
#' @import graphics grDevices parallel plyr Rcpp R2admb stats testthat truncnorm
#' @importFrom CircStats dvm rvm
#' @importFrom optimx optimx
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom fields image.plot
#' @importFrom matrixStats colProds
#' @importFrom mgcv gam
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom secr make.capthist make.mask read.mask read.traps sim.popn
#' @importFrom utils example setTxtProgressBar txtProgressBar
#' @importFrom viridis viridis
#' @useDynLib ascr
NULL

## Data documentation.

#' \emph{Arthroleptella lightfooti} survey data
#'
#' Data from an acoustic survey of the Western Cape moss frog
#' \emph{Arthroleptella lightfooti}. These data are from a 25 s subset
#' of the original recording, taken on 16 May 2012 at Silvermine,
#' Table Mountain National Park, Cape Town, South Africa. Acoustic
#' signal strengths and times of arrival were measured, and this
#' information is contained in the capture history object.
#'
#' This object is a list which contains components:
#' \itemize{
#' \item \code{capt}: A capture history object.
#' \item \code{traps}: A traps object.
#' \item \code{mask}: A suitable mask object.
#' \item \code{cutoff}: The microphone cutoff value.
#' \item \code{freqs}: A vector of call frequencies measured
#'                     independently to the acoustic survey.
#' }
#'
#' @name lightfooti
#' @format A list.
#' @usage lightfooti
#' @docType data
#' @keywords datasets
NULL

#' Example data
#'
#' This object contains simulated data with all types of supplementary
#' information, corresponding trap locations, and a suitable mask
#' object. Also included are some example model fits, which were
#' generated from these data using the \link{fit.ascr} function.
#'
#' This object is a list which contains components:
#' \itemize{
#' \item \code{capt}: A capture history object.
#' \item \code{traps}: A traps object.
#' \item \code{mask}: A suitable mask object.
#' \item \code{cutoff}: The cutoff value used to simluate these data.
#' \item \code{fits}: Some example model fits.
#' }
#'
#' @name example.data
#' @format A list.
#' @usage example.data
#' @docType data
#' @keywords datasets
NULL

#' Multi-array example data
#'
#' This object contains simulated data From multiple detector
#' arrays.
#'
#' This object is a list which contains components:
#' \itemize{
#' 
#' \item \code{capt}: A list of capture history objects, one from each
#' detector array.
#'
#' \item \code{traps}: A list of traps objects.
#' \item \code{mask}: A list of suitable mask objects.
#' }
#'
#' @name multi.example.data
#' @format A list.
#' @usage multi.example.data
#' @docType data
#' @keywords datasets
NULL
