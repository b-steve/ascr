## Package imports for roxygenise to pass to NAMESPACE.
#' @import CircStats Rcpp R2admb secr
#' @useDynLib admbsecr
NULL

#' Fitting SECR models in ADMB
#'
#' Fits an SECR model, with our without supplementary information
#' relevant to animal location. Parameter estimation is done by
#' maximum likelihood through and ADMB executible.
#'
#' @param capt A list with named components, containing the capture
#' history and supplementary information.
#' @param traps A matrix with two columns. The rows provide Cartesian
#' coordinates for trap locations.
#' @param mask A matrix with two columns. The rows provide Cartesian
#' coordinates for the mask point locations.
#' @param detfn A character string specifying the detection function
#' to be used. Options are "hn" (halfnormal), "hr" (hazard rate), "th"
#' (threshold), "lth" (log-link threshold), or "ss" (signal
#' strength). If the latter is used, signal strength information must
#' be provided in \code{capt}.
#' @param sv A named list. Component names are parameter names, and
#' each component is a start value for the associated parameter.
#' @param bounds A named list. Component names are parameter names,
#' and each components is a vector of length two, specifying the
#' bounds for the associated parameter.
#' @param fix A named list. Component names are parameter names to be
#' fixed, and each component is the fixed value for the associated
#' parameter.
#' @param scalefactors A named list. Component names are parameter
#' names, and each component is a scalefactor for the associated
#' parameter. The default behaviour is to automatically select
#' scalefactors based on parameter start values.
#' @param ss.link A character string, either \code{"indentity"} or
#' \code{"log"}, which specifies the link function for the signal
#' strength detection function. Only required when \code{detfn} is
#' \code{"ss"}.
#' @param cutoff The signal strength threshold, above which sounds are
#' identified as detections. Only required when \code{detfn} is
#' \code{"ss"}.
#' @param trace logical, if \code{TRUE} parameter values at each step
#' of the optimisation algorithm are printed to the R session.
#' @param clean logical, if \code{TRUE} ADMB output files are removed.
#' @export
#'
admbsecr <- function(capt, traps, mask, detfn = "hn", sv = NULL, bounds = NULL,
                     fix = NULL, scalefactors = NULL, ss.link = "identity",
                     cutoff = NULL, trace = FALSE, clean = TRUE){
  capt.bin <- capt$bincapt
  if (is.null(capt.bin)){
    stop("The binary capture history must be provided as a component of 'capt'.")
  }
  if (!is.list(sv) & !is.null(sv)){
    stop("The 'sv' argument must be 'NULL' or a list")
  }
  if (!is.list(bounds) & !is.null(bounds)){
    stop("The 'bounds' argument must be 'NULL' or a list")
  }
  if (!is.list(fix) & !is.null(fix)){
    stop("The 'fix' argument must be 'NULL' or a list")
  }
  n <- nrow(capt.bin)
  n.traps <- nrow(traps)
  n.mask <- nrow(mask)
  A <- attr(mask, "area")
  ## Removing attributes from mask.
  mask <- as.matrix(mask)
  ## TODO: Sort out how to determine supplementary parameter names.
  supp.types <- c("ang", "dist", "ss", "toa", "mrds")
  fit.types <- supp.types %in% names(capt)
  names(fit.types) <- supp.types
  ## Logical indicators for additional information types.
  fit.angs <- fit.types["ang"]
  fit.dists <- fit.types["dist"]
  fit.ss <- fit.types["ss"]
  fit.toas <- fit.types["toa"]
  fit.mrds <- fit.types["mrds"]
  ## Capture histories for additional information types (if they exist)
  capt.ang <- if (fit.angs) capt$ang else 0
  capt.dist <- if (fit.dists) capt$dist else 0
  capt.ss <- if (fit.ss) capt$ss else 0
  capt.toa <- if (fit.toas) capt$toa else 0
  mrds.dist <- if (fit.mrds) capt$mrds else 0
  suppar.names <- c("kappa", "alpha", "sigma.toa")[fit.types[c("ang", "dist", "toa")]]
  if (fit.ss){
    if (ss.link == "identity"){
      detfn <- "ss"
      linkfn.id <- 1
    } else if (ss.link == "log"){
      detfn <- "log.ss"
      linkfn.id <- 2
    } else {
      stop("ss.link must be either \"identity\" or \"log\"")
    }
  } else {
    ## Not sure what a linkfn.id of 3 means? Probably throws an error in ADMB.
    linkfn.id <- 3
  }
  detfns <- c("hn", "hr", "th", "lth", "ss", "log.ss")
  ## Sets detection function ID number for use in ADMB:
  ## 1 = Half normal
  ## 2 = Hazard rate
  ## 3 = Threshold
  ## 4 = Log-link threshold
  ## 5 = Identity-link signal strength
  ## 6 = Log-link signal strength.
  detfn.id <- which(detfn == detfns)
  detpar.names <- switch(detfn,
                         hn = c("g0", "sigma"),
                         hr = c("g0", "sigma", "z"),
                         th = c("shape", "scale"),
                         lth = c("shape.1", "shape.2", "scale"),
                         ss = c("b0.ss", "b1.ss", "sigma.ss"),
                         log.ss = c("b0.ss", "b1.ss", "sigma.ss"))
  par.names <- c("D", detpar.names, suppar.names)
  n.detpars <- length(detpar.names)
  n.suppars <- length(suppar.names)
  npars <- length(par.names)
  ## Sorting out start values. Start values are set to those provided,
  ## or else are determined automatically from functions in
  ## autofuns.r.
  sv.old <- sv
  sv <- vector("list", length = npars)
  names(sv) <- par.names
  sv[names(sv.old)] <- sv.old
  sv[names(fix)] <- fix
  auto.names <- par.names[sapply(sv, is.null)]
  sv.funs <- paste("auto", auto.names, "2", sep = "")
  for (i in seq(1, length(auto.names), length.out = length(auto.names))){
    sv[auto.names[i]] <- eval(call(sv.funs[i], capt, traps))
  }
  ## Sorting out phases.
  ## TODO: Add phases parameter so that these can be controlled by user.
  phases <- vector("list", length = npars)
  for (i in par.names){
    if (any(i == names(fix))){
      ## Phase of -1 in ADMB fixes parameter at starting value.
      phases[[i]] <- -1
    } else {
      phases[[i]] <- 0
    }
  }
  D.phase <- phases[["D"]]
  detpars.phase <- c(phases[detpar.names], recursive = TRUE)
  if (n.suppars > 0){
    suppars.phase <- c(phases[suppar.names], recursive = TRUE)
  } else {
    suppars.phase <- -1
  }
  ## Sorting out bounds.
  ## Below bounds are the defaults.
  default.bounds <- list(D = c(0, 1e8),
                         D.a = c(0, 1e8),
                         mu.C = c(0, 1e8),
                         sigma.C = c(0, 1e5),
                         g0 = c(0, 1),
                         sigma = c(0, 1e5),
                         shape = c(-1e8, 1e8),
                         shape.1 = c(0, 1e5),
                         shape.2 = c(-1e8, 1e8),
                         scale = c(0, 1e5),
                         b0.ss = c(0, 1e8),
                         b1.ss = c(0, 10),
                         sigma.ss = c(0, 1e5),
                         z = c(0, 1e5),
                         sigma.toa = c(0, 1e5),
                         kappa = c(0, 700),
                         alpha = c(0, 10000))[par.names]
  bound.changes <- bounds
  bounds <- default.bounds
  for (i in names(default.bounds)){
    if (i %in% names(bound.changes)){
      bounds[[i]] <- bound.changes[[i]]
    }
  }
  D.bounds <- bounds[["D"]]
  D.lb <- D.bounds[1]
  D.ub <- D.bounds[2]
  detpar.bounds <- bounds[detpar.names]
  detpars.lb <- sapply(detpar.bounds, function(x) x[1])
  detpars.ub <- sapply(detpar.bounds, function(x) x[2])
  if (n.suppars > 0){
    suppar.bounds <- bounds[suppar.names]
    suppars.lb <- sapply(suppar.bounds, function(x) x[1])
    suppars.ub <- sapply(suppar.bounds, function(x) x[2])
  } else {
    suppars.lb <- 0
    suppars.ub <- 0
  }
  ## Sorting out scalefactors.
  ## TODO: Sort these out in a better way.
  if (is.null(scalefactors)){
    sv.vec <- c(sv, recursive = TRUE)
    ## Currently, by default, the scalefactors are the inverse
    ## fraction of each starting value to the largest starting
    ## value. Not sure how sensible this is.
    sf <- max(sv.vec)/sv.vec
  } else {
    sf <- numeric(npars)
    names(sf) <- par.names
    for (i in par.names){
      sf[i] <- ifelse(i %in% names(scalefactors), scalefactors[i], 1)
    }
  }
  D.sf <- sf[["D"]]
  detpars.sf <- c(sf[detpar.names], recursive = TRUE)
  if (n.suppars > 0){
    suppars.sf <- c(sf[suppar.names], recursive = TRUE)
  } else {
    suppars.sf <- 1
  }
  ## Setting small number so that numerical under/overflow in ADMB
  ## does not affect estimation.
  dbl.min <- 1e-150
  ## Some stuff being set as defaults for testing.
  n.freqs <- 1
  call.freqs <- 1
  ## Calculating distances and angles.
  ## TODO: Try calculating these in the PROCEDURE_SECTION instead.
  dists <- distances(traps, mask)
  if (fit.angs){
    angs <- bearings(traps, mask)
  } else {
    angs <- 0
  }
  if (fit.toas){
    toa.ssq <- make_toa_ssq(capt$toa, dists)
  } else {
    toa.ssq <- 0
  }
  if (is.null(cutoff)){
    cutoff <- 0
  }
  ## Kludge to fix number of parameters for no supplementary
  ## information.
  if (n.suppars == 0){
    n.suppars <- max(c(n.suppars, 1))
    sv$dummy <- 0
  }
  ## Stuff for the .dat file.
  data.list <- list(D_lb = D.lb, D_ub = D.ub, D_phase = D.phase, D_sf
                    = D.sf, n_detpars = n.detpars, detpars_lb =
                    detpars.lb, detpars_ub = detpars.ub, detpars_phase
                    = detpars.phase, detpars_sf = detpars.sf,
                    n_suppars = n.suppars, suppars_lb = suppars.lb,
                    suppars_ub = suppars.ub, suppars_phase =
                    suppars.phase, suppars_sf = suppars.sf, detfn_id =
                    detfn.id, trace = as.numeric(trace), DBL_MIN =
                    dbl.min, n = n, n_traps = n.traps, n_mask =
                    n.mask, A = A, n_freqs = n.freqs, call_freqs =
                    call.freqs, capt_bin = capt.bin, fit_angs =
                    as.numeric(fit.angs), capt_ang = capt.ang,
                    fit_dists = as.numeric(fit.dists), capt_dist =
                    capt.dist, fit_ss = as.numeric(fit.ss), cutoff =
                    cutoff, linkfn_id = linkfn.id, capt_ss = capt.ss,
                    fit_toas = as.numeric(fit.toas), capt_toa =
                    capt.toa, fit_mrds = as.numeric(fit.mrds),
                    mrds_dist = mrds.dist, dists = dists, angs = angs,
                    toa_ssq = toa.ssq)
  ## Idea of running executable as below taken from glmmADMB.
  ## Working out correct command to run from command line.
  if (.Platform$OS == "windows"){
    os.type <- "windows"
    cmd <- "secr -ind secr.dat -ainp secr.pin"
  } else if (.Platform$OS == "unix"){
    os.type <- "linux"
    cmd <- "./secr -ind secr.dat -ainp secr.pin"
  } else {
    stop("OS not supported yet.")
  }
  ## Finding executable folder (possible permission problems?).
  exe.dir <- paste(system.file(package = "admbsecr"), "ADMB", "bin", os.type, sep = "/")
  curr.dir <- getwd()
  ## Moving to executable location.
  setwd(exe.dir)
  curr.files <- list.files()
  ## Creating .pin and .dat files.
  write_pin("secr", sv)
  write_dat("secr", data.list)
  ## Running ADMB executable.
  system(cmd, ignore.stdout = !trace)
  ## Reading in model results.
  out <- read.admbsecr("secr")
  ## Cleaning up files.
  all.files <- list.files()
  new.files <- all.files[!all.files %in% curr.files]
  if (clean){
    file.remove(new.files)
  }
  ## Moving back to original working directory.
  setwd(curr.dir)
  out
}

