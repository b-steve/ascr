## Package imports for roxygenise to pass to NAMESPACE.
#' @import CircStats R2admb secr
NULL

#' Fitting SECR models in ADMB
#'
#' Fits an SECR model, with our without supplementary information relevant to animal
#' location. Parameter estimation is done by maximum likelihood in ADMB.
#'
#' ADMB is called to fit an SECR model through use of the R2admb package. Different
#' methods are used depending on the additional information on animal location that is
#' collected. These are:
#' \itemize{
#'    \item \code{"simple"}: Normal SECR with no additional information. Parameters to
#'      estimate are:
#'   \itemize{
#'          \item D:      Animal density.
#'
#'          \item g0:     Probability of detection at distance 0.
#'
#'          \item sigma:  'Standard deviation' parameter for halfnormal detection function.
#'   }
#'     \item \code{"toa"}: SECR with precise time of arrival (TOA) recorded. Parameters to
#'      estimate are:
#'   \itemize{
#'          \item D:        As above.
#'
#'          \item g0:       As above.
#'
#'          \item sigma:    As above.
#'
#'          \item sigmatoa: Error term associated with TOA.
#'   }
#'    \item \code{"ang"}: SECR with estimates of angle to animal recorded. Parameters to
#'      estimate are:
#'    \itemize{
#'          \item D:        As above.
#'
#'          \item g0:       As above.
#'
#'          \item sigma:    As above.
#'
#'          \item kappa:    Error term associated with angle estimation.
#'    }
#'    \item \code{"ss"}: SECR with received signal strengths at traps recorded. Parameters
#'      to estimate are:
#'    \itemize{
#'          \item D:        As above.
#'
#'          \item ssb0:     Average signal strength at sound source.
#'
#'          \item ssb1:     Decrease in signal strength per unit distance due to sound
#'                      propagation.
#'
#'          \item sigmass:  Error term associated with signal strength.
#'    }
#'    \item \code{"sstoa"}: SECR with precise TOA and received signal strengths at traps
#'      recorded. Parameters to estimate are:
#'    \itemize{
#'          \item D:        As above.
#'
#'          \item sigmatoa: As above.
#'
#'          \item ssb0:     As above.
#'
#'          \item ssb1:     As above.
#'
#'          \item sigmass:  As above.
#'    }
#' }
#' @param capt an array of dimension \code{(n, S, K)}, where \code{n} is the number of
#' detected animals, \code{S} is number of individual sampling sessions, and \code{K}
#' is the number of deployed traps. The object returned by  \code{make.capthist()} is
#' suitable if \code{method} is \code{"simple"}. Otherwise, the \code{1} values in
#' this array must be changed to the value of the recorded supplementary information,
#' which will depend on \code{method} (see 'Details'). When \code{method} is
#' \code{"sstoa"}, this array must be of dimension \code{(n, S, K, 2)}, where
#' \code{capt[, , , 1]} provides the signal strength information and
#' \code{capt[, , , 2]} provides the time of arrival information.
#' @param traps a matrix containing the coordinates of trap locations. The object
#' returned by \code{read.traps()} is suitable.
#' @param mask a mask object. The object returned by \code{make.mask()} is suitable.
#' @param sv either "auto", or a named vector of starting values for each of the model's
#' parameters. See 'Details' for list of parameters used by each method. If a named
#' vector, any of the elements can be set to \code{"auto"}, which will result in a
#' (hopefully) sensible starting value being automatically generated for the parameters
#' concerned. If \code{sv} itself is \code{"auto"}, all starting values are automatically
#' generated.
#' @param ssqtoa an optional matrix. If calculated before call to \code{admbsecr},
#' providing this will prevent recalculation.
#' @param cutoff The signal strength threshold of detection. Required if \code{method} is
#' \code{"ss"} or \code{"sstoa"}.
#' @param admbwd file path to the ADMB working directory. Only required if
#' \code{autogen} is \code{TRUE}, in which case it points to the directory in which the
#' \code{.tpl} file is located.
#' @param method either \code{"simple"}, \code{"toa"}, \code{"ang"}, \code{"ss"}, or
#' \code{"sstoa"}. See 'Details'.
#' @param memory value of \code{arrmblsize} in ADMB. Increase this if ADMB reports a
#' memory error.
#' @param profpars character vector of names of parameters over which profile likelihood
#' should occur (untested and probably very slow!).
#' @param clean logical, if \code{TRUE} ADMB files are cleaned after fitting of the model.
#' @param verbose logical, if \code{TRUE} ADMB details, along with error messages, are
#' printed to the R session.
#' @param trace logical, if \code{TRUE} parameter values at each step of the fitting
#' algorithm are printed to the R session.
#' @param autogen logical, if \code{TRUE}, the appropriate \code{.tpl} file is written
#' to \code{admbwd} (or the current working directory if \code{admbwd} is \code{NULL}).
#' If \code{FALSE}, the \code{.tpl} file already be located in \code{admbwd} (or the
#' current working directory if \code{admb} is \code{NULL}). Usually only set to
#' \code{FALSE} for development purposes.
#' @return An object of class 'admb'.
#' @author Ben Stevenson
#' @seealso \code{\link[R2admb]{do_admb}}, \code{\link[secr]{secr.fit}},
#' \code{\link[secr]{make.capthist}}, \code{\link[secr]{read.traps}}.
#' @export
admbsecr <- function(capt, traps = NULL, mask, sv = "auto", ssqtoa = NULL, cutoff = NULL,
                     admbwd = NULL, method = "simple", memory = NULL, profpars = NULL,
                     clean = TRUE, verbose = FALSE, trace = FALSE, autogen = TRUE){
  ## Warnings for incorrect input.
  if (length(method) != 1){
    stop("method must be of length 1")
  }
  if (method == "simple" & any(capt != 1 & capt != 0)){
    stop('capt must be binary when using the "simple" method')
  }
  if (method == "ss" & is.null(cutoff)){
    stop("cutoff must be supplied for signal strength analysis")
  }
  if (trace){
    verbose <- TRUE
  }
  trace <- as.numeric(trace)
  currwd <- getwd()
  ## If traps is NULL, see if it is provided as part of capt.
  if (is.null(traps)){
    traps <- traps(capt)
  }
  ## Moving to ADMB working directory.
  if (!is.null(admbwd)){
    setwd(admbwd)
  }
  if (autogen){
    prefix <- "secr"
    make.all.tpl.easy(memory = memory, methods = method)
  } else {
    prefix <- paste(method, "secr", sep = "")
  }
  ## If NAs are present in capture history object, change to zeros.
  capt[is.na(capt)] <- 0
  ## Extracting no. animals trapped (n) and traps (k) from capture history array.
  ## Only currently works with one capture session.
  n <- dim(capt)[1]
  k <- dim(capt)[3]
  ## Area covered by each mask location.
  A <- attr(mask, "area")
  bincapt <- capt
  bincapt[capt > 0] <- 1
  if (length(dim(bincapt)) == 4){
    bincapt <- bincapt[, , , 1, drop = FALSE]
  } else if (length(dim(bincapt)) > 4){
    stop("capt array cannot have more than 4 dimensions.")
  }
  ## Setting number of model parameters.
  npars <- c(3[method == "simple"], 4[method %in% c("toa", "ang", "ss")], 5[method == "sstoa"])
  ## Setting sensible start values if elements of sv are "auto".
  if (length(sv) == 1 & sv[1] == "auto"){
    sv <- rep("auto", npars)
  }
  if (any(sv == "auto")){
    ## Give sv vector names if it doesn't have them.
    if (is.null(names(sv))){
      names(sv) <- c("D", "g0"[!(method == "ss" | method == "sstoa")],
                     "sigma"[!(method == "ss" | method == "sstoa")],
                     "sigmatoa"[method == "toa" | method == "sstoa"],
                     "kappa"[method == "ang"],
                     "ssb0"[method == "ss" | method == "sstoa"],
                     "ssb1"[method == "ss" | method == "sstoa"],
                     "sigmass"[method == "ss" | method == "sstoa"])
    } else {
      ## Reordering sv vector if names are provided.
      sv <- sv[c("D", "g0"[!(method == "ss" | method == "sstoa")],
                 "sigma"[!(method == "ss" | method == "sstoa")],
                 "sigmatoa"[method == "toa" | method == "sstoa"],
                 "kappa"[method == "ang"],
                 "ssb0"[method == "ss" | method == "sstoa"],
                 "ssb1"[method == "ss" | method == "sstoa"],
                 "sigmass"[method = "ss" | method == "sstoa"])]
    }
    autofuns <- list("D" = autoD, "g0" = autog0, "sigma" = autosigma,
                     "sigmatoa" = autosigmatoa, "kappa" = autokappa,
                     "ssb0" = autossb0, "ssb1" = autossb1,
                     "sigmass" = autosigmass)
    ## Replacing "auto" elements of sv vector.
    for (i in rev(which(sv == "auto"))){
      sv[i] <- autofuns[[names(sv)[i]]](capt, bincapt, traps, mask, sv, cutoff, method)
    }
    sv <- as.numeric(sv)
  }
  ## Removing attributes from capt and mask objects as do_admb cannot handle them.
  bincapt <- matrix(as.vector(bincapt), nrow = n, ncol = k)
  capt <- array(as.vector(capt), dim = c(n, k, dim(capt)[4][length(dim(capt)) == 4]))
  mask <- as.matrix(mask)
  ## No. of mask locations.
  nm <- nrow(mask)
  traps <- as.matrix(traps)
  ## Distances between traps and mask locations.
  dist <- distances(traps, mask)
  ## Setting up parameters for do_admb.
  if (method == "simple"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt,
                 dist = dist, traps = traps, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3])
    bounds <- list(D = c(0, 10000000), g0 = c(0, 1), sigma = c(0, 100000))
  } else if (method == "toa"){
    if (is.null(ssqtoa)){
      ssqtoa <- apply(capt, 1, toa.ssq, dists = dist)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, toacapt = capt,
                 toassq = t(ssqtoa), dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], sigmatoa = sv[4])
    bounds <- list(D = c(0, 10000000), g0 = c(0, 1), sigma = c(0, 100000),
                   sigmatoa = c(0, 100000))
  } else if (method == "ang"){
    angs <- angles(traps, mask)
    data <- list(n = n, ntraps = k, nmask = nm, A = A, angcapt = capt,
                 ang = angs, dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], kappa = sv[4])
    bounds <- list(D = c(0, 10000000), g0 = c(0, 1), sigma = c(0, 100000),
                   kappa = c(0, 700))
  } else if (method == "ss"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, c = cutoff, sscapt = capt,
                 dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], ssb0 = sv[2], ssb1 = sv[3], sigmass = sv[4])
    bounds <- list(D = c(0, 10000000), sigmass = c(0, 100000), ssb1 = c(-100000, 0))
  } else if (method == "sstoa"){
    if (is.null(ssqtoa)){
      ssqtoa <- apply(capt[, , 1], 1, toa.ssq, dists = dist)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, c = cutoff, sscapt = capt[, , 1],
                 toacapt = capt[, , 2], toassq = t(ssqtoa), dist = dist, capt = bincapt,
                 trace = trace)
    params <- list(D = sv[1], sigmatoa = sv[2], ssb0 = sv[3], ssb1 = sv[4], sigmass = sv[5])
    bounds <- list(D = c(0, 10000000), sigmass = c(0, 100000), ssb1 = c(-100000, 0),
                   sigmatoa = c(0, 100000))
  } else {
    stop('method must be either "simple", "toa", "ang", "ss", or "sstoa"')
  }
  ## Fitting the model.
  if (!is.null(profpars)){
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = verbose,
                   profile = FALSE, profpars = profpars, safe = TRUE,
                   run.opts = run.control(checkdata = "write", checkparam = "write",
                     clean_files = clean))
  } else {
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = verbose,
                   safe = FALSE,
                   run.opts = run.control(checkdata = "write", checkparam = "write",
                     clean_files = clean))
  }
  if (autogen){
    file.remove("secr.tpl")
  }
  setwd(currwd)
  fit
}
