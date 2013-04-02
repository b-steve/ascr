## Package imports for roxygenise to pass to NAMESPACE.
#' @import CircStats R2admb secr
NULL

#' Fitting SECR models in ADMB
#'
#' THE FOLLOWING INFORMATION IS OUT OF DATE AND REQUIRES A REWRITE.
#'
#' Fits an SECR model, with our without supplementary information relevant to animal
#' location. Parameter estimation is done by maximum likelihood in ADMB.
#'
#' ADMB is called to fit an SECR model through use of the R2admb package. Different
#' methods are used depending on the additional information on animal location that is
#' collected. Different detection functions are used depending on the relationship between
#' detection probability and distance from detector.
#'
#' Note that the method \code{"ss"} is a special case in that it incorporates its own detection
#' function, and thus the half normal, hazard rate (etc) options cannot be specified. Instead,
#' either \code{"identity"} (the default) or \code{"log"} can be provided for the argument
#' \code{detfn}, which give the link function for the estimated received signal strengths.
#'
#' The parameter D, density of animals (in individuals per hectare) is always estimated.
#' The other parameters in the model depend on the method and the detection function used.
#'
#' Possible methods, along with their parameters, are as follows:
#'
#' \itemize{
#'    \item \code{"simple"}: Normal SECR with no additional information. No additional parameters.
#'     \item \code{"toa"}: SECR with precise time of arrival (TOA) recorded:
#'   \itemize{
#'          \item sigmatoa: Error term associated with the normal distribution used to model TOA.
#'   }
#'    \item \code{"ang"}: SECR with estimates of angle to animal recorded:
#'    \itemize{
#'          \item kappa:    Error term from a Von-Mises distribution, used to model estimated
#'                       angles.
#'    }
#'    \item \code{"ss"}: SECR with received signal strengths at traps recorded:
#'    \itemize{
#'          \item ssb0:     Signal strength at source.
#'
#'          \item ssb1:     Decrease in signal strength per unit distance due to sound
#'                      propagation.
#'
#'          \item sigmass:  Error term associated with the normal distribution used to model signal strength.
#'    }
#'    \item \code{"sstoa"}: SECR with precise TOA and received signal strengths at traps
#'      recorded:
#'    \itemize{
#'          \item sigmatoa: As above.
#'
#'          \item ssb0:     As above.
#'
#'          \item ssb1:     As above.
#'
#'          \item sigmass:  As above.
#'    }
#'    \item \code{"dist"}: SECR with estimated distances between animals
#'      and traps at which detections occurred:
#'    \itemize{
#'          \item alpha:    Shape parameter associated with the gamma distribution
#'                          used to model estimated distances.
#'    }
#' }
#' Possible detection functions, along with their parameters, are as follows:
#'
#' \itemize{
#'     \item \code{"hn"}: Half-normal detection function:
#'   \itemize{
#'          \item g0:       Probability of detection at distance 0.
#'
#'          \item sigma:    Scale parameter.
#'   }
#'    \item \code{"hr"}: Hazard rate detection function:
#'    \itemize{
#'          \item g0
#'
#'          \item sigma
#'
#'          \item z
#'    }
#'    \item \code{"th"}: Threshold detection function:
#'    \itemize{
#'          \item shape
#'
#'          \item scale
#'    }
#'    \item \code{"logth"}: Log-link threshold detection function:
#'    \itemize{
#'          \item shape1
#'
#'          \item shape2
#'
#'          \item scale
#'    }
#' }
#' @param capt an array of dimension \code{(n, S, K)}, where \code{n} is the number of
#' detected animals, \code{S} is number of individual sampling sessions, and \code{K}
#' is the number of deployed traps. The object returned by  \code{make.capthist()} is
#' suitable if \code{method} is \code{"simple"}. Otherwise, the \code{1} values in
#' this array must be changed to the value of the recorded supplementary information,
#' which will depend on \code{method} (see 'Details'). When \code{method} is
#' \code{"sstoa"} or \code{"mrds"}, this array must be of dimension \code{(n, S, K, 2)}.
#' With \code{"sstoa"},  \code{capt[, , , 1]} provides the signal strength information and
#' \code{capt[, , , 2]} provides the time of arrival information. With \code{"mrds"},
#' \code{capt[, , , 1]} provides the binary capture history array and \code{capt[, , , 2]}
#' provides the distances between all traps (regardless of capture) and detected animals.
#' @param traps a matrix containing the coordinates of trap locations. The object
#' returned by \code{\link[secr]{read.traps}} is suitable.
#' @param mask a mask object. The object returned by \code{\link[secr]{make.mask}} is
#' suitable.
#' @param sv either \code{"auto"}, or a named vector. If \code{auto}, starting values for
#' all parameters are automatically generated. If a vector, the elements are the starting
#' values and the names indicate which parameter these correspond to. Starting values for
#' all parameters need not be provided; they are automatically generated for any parameters
#' that do not have a starting value explicitly provided. See 'Details' for list of
#' parameters used by each method.
#' @param bounds a list with optional components corresponding to parameters which are to
#' have their default bounds overridden. Each component should be a vector of length two
#' specifying the bounds, and the name of the component should be the name of the
#' parameter to which these bounds apply. To remove default bounds from a parameter the
#' component should be \code{NULL} rather than a vector of length two. Bounds for all
#' parameters need not be provided; if there is no component corresponding to a
#' parameter it keeps its default bounds. See 'Details' for list of parameters used by
#' each method.
#' @param fix a list with optional components corresponding to parameters which are to
#' be fixed rather than estimated. Each component should be a vector of length one,
#' specifying the fixed value of the parameter, and the name of the component should be
#' the name of the paramter to which this value applies.
#' @param ssqtoa an optional matrix. If calculated before call to \code{admbsecr},
#' providing this will prevent recalculation.
#' @param cutoff The signal strength threshold of detection. Required if \code{method} is
#' \code{"ss"} or \code{"sstoa"}.
#' @param admbwd file path to the ADMB working directory. Only required if
#' \code{autogen} is \code{TRUE}, in which case it points to the directory in which the
#' \code{.tpl} file is located.
#' @param method either \code{"simple"}, \code{"toa"}, \code{"ang"}, \code{"ss"}, or
#' \code{"sstoa"}. See 'Details'.
#' @param detfn the detection function to be used. Either half normal (\code{"hn"}),
#' hazard rate (\code{"hr"}), threshold (\code{"th"}) or log-link threshold (\code{"logth"}.
#' If method is \code{"ss"}, this argument gives the link function for the expected received
#' signal strengths (either \code{"identity"}, the default, or \code{"log"}).
#' @param memory value of \code{arrmblsize} in ADMB. Increase this if ADMB reports a
#' memory error.
#' @param profpars character vector of names of parameters over which profile likelihood
#' should occur.
#' @param clean logical, if \code{TRUE} ADMB files are cleaned after fitting of the model.
#' @param verbose logical, if \code{TRUE} ADMB details, along with error messages, are
#' printed to the R session.
#' @param trace logical, if \code{TRUE} parameter values at each step of the fitting
#' algorithm are printed to the R session.
#' @param autogen logical, if \code{TRUE}, the appropriate \code{.tpl} file is written
#' to \code{admbwd} (or the current working directory if \code{admbwd} is \code{NULL}).
#' If \code{FALSE}, the \code{.tpl} file should already be located in \code{admbwd} (or
#' the current working directory if \code{admb} is \code{NULL}). Usually only set to
#' \code{FALSE} for development purposes.
#' @return An object of class 'admb'.
#'
#' The following functions can be used to extract model components:
#' \code{\link[base]{summary}}, \code{\link[R2admb:AIC.admb]{AIC}},
#' \code{\link[R2admb:AIC.admb]{logLik}}, \code{\link[R2admb:AIC.admb]{deviance}},
#' \code{\link[R2admb:AIC.admb]{vcov}}, \code{\link[R2admb:AIC.admb]{coef}},
#' \code{\link[R2admb:AIC.admb]{stdEr}}, and \code{\link[R2admb:AIC.admb]{confint}}.
#'
#' The latter takes arguments \code{level} and \code{method}, which specify the confidence
#' level and calculation method respectively. The default method gives quadratic (Wald)
#' intervals based on approximate standard errors; \code{"profile"} gives profile
#' likelihood intervals, and can be used if the \code{admbsecr()} parameter
#' \code{profpars} is non-null and provides names of model parameters that are to be
#' profiled.
#' @author Ben Stevenson
#' @export
admbsecr <- function(capt, traps = NULL, mask, sv = "auto", bounds = NULL, fix = NULL,
                     ssqtoa = NULL, cutoff = NULL, admbwd = NULL, method = "simple",
                     detfn = "hn" , memory = NULL, profpars = NULL, clean = TRUE,
                     verbose = FALSE, trace = FALSE, autogen = TRUE){
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
  if (!is.array(capt) | !(length(dim(capt)) == 3 | length(dim(capt)) == 4)){
    stop("capt must be a three or four-dimensional array.")
  }
  if (dim(capt)[2] != 1){
    stop("admbsecr only currently works for a single sampling session.")
  }
  if (method == "ss" | method == "sstoa"){
    if (missing(detfn)){
      detfn <- "identity"
    } else if (!(detfn == "identity" | detfn == "log"))
      stop("The \"ss\" and \"sstoa\" methods use their own detection function. \nThe 'detfn' argument can either be \"identity\" or \"log\" (see 'Details' in help file).")
  } else if (!(detfn == "hn" | detfn == "th" | detfn == "logth" | detfn == "hr")){
    stop("Detection function must be \"hn\", \"th\", \"logth\" or \"hr\"")
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
  if (diff(range(traps[, 1])) == 0 & diff(range(traps[, 2])) == 0 & any(sv == "auto")){
    stop("All traps are at the same location; please provide starting values.")
  }
  ## Moving to ADMB working directory.
  if (!is.null(admbwd)){
    setwd(admbwd)
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
  ## Detection function parameters.
  detnames <- c(c("g0", "sigma")[detfn == "hn" | detfn == "hr"],
                c("shape", "scale")[detfn == "th"],
                c("shape1", "shape2", "scale")[detfn == "logth"],
                "z"[detfn == "hr"])
  ## Parameter names.
  parnames <- c("D", detnames,
                c("ssb0", "ssb1", "sigmass")[method == "ss" | method == "sstoa"],
                "sigmatoa"[method == "toa" | method == "sstoa"],
                "kappa"[method == "ang"],
                "alpha"[method == "dist"])
  ## Setting number of model parameters.
  npars <- length(parnames)
  ## Setting up bounds.
  default.bounds <- list(D = c(0, 1e8),
                         g0 = c(0, 1),
                         sigma = c(0, 1e5),
                         shape = NULL,
                         shape1 = c(0, 1e5),
                         shape2 = NULL,
                         scale = c(-10, 0),
                         ssb0 = NULL,
                         ssb1 = c(-10, 0),
                         sigmass = c(0, 1e5),
                         z = NULL,
                         sigmatoa = c(0, 1e5),
                         kappa = c(0, 700),
                         alpha = c(0, 150))[parnames]
  if (!(is.list(bounds) | is.null(bounds))){
    stop("bounds must either be NULL or a list.")
  } else {
    bound.changes <- bounds
    bounds <- default.bounds
    for (i in names(default.bounds)){
      if (i %in% names(bound.changes)){
        bounds[[i]] <- bound.changes[[i]]
      } else {
        ## Removing NULL elements from list.
        bounds[[i]] <- bounds[[i]]
      }
    }
  }
  ## If sv is a list, turn it into a vector.
  if (is.list(sv)){
    sv <- c(sv, recursive = TRUE)
  }
  ## Setting sv to a vector full of "auto" if required.
  if (length(sv) == 1 & sv[1] == "auto"){
    sv <- rep("auto", npars)
    names(sv) <- parnames
  } else if (is.null(names(sv))){
    stop("sv is not a named vector.")
  } else if (length(unique(names(sv))) != length(names(sv))){
    stop("sv names are not all unique")
  } else {
    ## Warning if a listed parameter name is not used in this model.
    if (!all(names(sv) %in% parnames)){
      warning("One of the element names of sv is not a parameter used in this model.")
    }
    sv.old <- sv
    sv <- rep("auto", npars)
    names(sv) <- parnames
    for (i in parnames){
      if (any(names(sv.old) == i)){
        sv[i] <- sv.old[i]
      }
    }
    ## Reordering sv vector.
    sv <- sv[parnames]
  }
  ## Creating .tpl file.
  if (autogen){
    prefix <- "secr"
    make.all.tpl.easy(memory = memory, method = method,
                      detfn = detfn, parnames = parnames)
    bessel.exists <- file.access("bessel.cxx", mode = 0)
    if (bessel.exists == -1){
      make.bessel()
    }
  } else {
    prefix <- paste(method, "secr", sep = "")
  }
  ## Adding fixed parameters to "sv" in case they are required for
  ## determining further start values.
  for (i in names(fix)){
    sv[i] <- fix[[i]]
  }
  autofuns <- list("D" = autoD, "g0" = autog0, "sigma" = autosigma,
                   "shape" = autoshape, "shape1" = autoshape1, "shape2" = autoshape2,
                   "scale" = autoscale, "z" = autoz,
                   "ssb0" = autossb0, "ssb1" = autossb1,
                   "sigmass" = autosigmass, "sigmatoa" = autosigmatoa,
                   "kappa" = autokappa, "alpha" = autoalpha)
  ## Replacing "auto" elements of sv vector.
  for (i in rev(which(sv == "auto"))){
    sv[i] <- autofuns[[names(sv)[i]]](capt, bincapt, traps, mask, sv, cutoff, method, detfn)
  }
  sv <- as.numeric(sv)
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
                 dist = dist, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3])
  } else if (method == "toa"){
    if (is.null(ssqtoa)){
      ssqtoa <- apply(capt, 1, toa.ssq, dists = dist)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, toacapt = capt,
                 toassq = t(ssqtoa), dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], sigmatoa = sv[4])
  } else if (method == "ang"){
    angs <- angles(traps, mask)
    data <- list(n = n, ntraps = k, nmask = nm, A = A, angcapt = capt,
                 ang = angs, dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], kappa = sv[4])
  } else if (method == "ss"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, c = cutoff, sscapt = capt,
                 dist = dist, capt = bincapt, trace = trace)
    params <- list(D = sv[1], ssb0 = sv[2], ssb1 = sv[3], sigmass = sv[4])
  } else if (method == "sstoa"){
    if (is.null(ssqtoa)){
      ssqtoa <- apply(capt[, , 1], 1, toa.ssq, dists = dist)
    }
    data <- list(n = n, ntraps = k, nmask = nm, A = A, c = cutoff, sscapt = capt[, , 1],
                 toacapt = capt[, , 2], toassq = t(ssqtoa), dist = dist, capt = bincapt,
                 trace = trace)
    params <- list(D = sv[1], sigmatoa = sv[2], ssb0 = sv[3], ssb1 = sv[4], sigmass = sv[5])
  } else if (method == "dist"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, distcapt = capt, dist = dist,
                 capt = bincapt, trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3], alpha = sv[4])
  } else if (method == "mrds"){
    data <- list(n = n, ntraps = k, nmask = nm, A = A, capt = capt[, , 1],
                 dist = dist, indivdist = capt[, , 2], trace = trace)
    params <- list(D = sv[1], g0 = sv[2], sigma = sv[3])
  } else {
    stop('method must be either "simple", "toa", "ang", "ss", "sstoa", "dist", or "mrds"')
  }
  params <- list()
  for (i in 1:npars){
    params[[i]] <- sv[i]
  }
  names(params) <- parnames
  ## Removing fixed parameters from param list and adding them to the data instead.
  for (i in names(fix)){
    params[[i]] <- NULL
    bounds[[i]] <- NULL
    data[[i]] <- fix[[i]]
  }
  ## Fitting the model.
  if (!is.null(profpars)){
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = verbose,
                   profile = TRUE, profpars = profpars, safe = FALSE,
                   run.opts = run.control(checkdata = "write", checkparam = "write",
                     clean_files = clean))
  } else {
    fit <- do_admb(prefix, data = data, params = params, bounds = bounds, verbose = verbose,
                   safe = FALSE,
                   run.opts = run.control(checkdata = "write",
                     checkparam = "write", clean_files = clean))
  }
  if (autogen){
    file.remove("secr.tpl")
    if (bessel.exists == -1){
      file.remove("bessel.cxx")
    }
  }
  setwd(currwd)
  fit$data <- data
  fit$traps <- traps
  fit$mask <- mask
  fit$method <- method
  fit$detfn <- detfn
  fit$parnames <- parnames
  class(fit) <- c(class(fit), method, "admbsecr")
  fit
}
