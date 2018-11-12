#' Fitting acoustic SCR models
#'
#' Fits an acoustic SCR model. Parameter estimation is done by maximum
#' likelihood through an AD Model Builder (ADMB) executable.
#'
#' ADMB uses a quasi-Newton method to find maximum likelihood
#' estimates for the model parameters. Standard errors are calculated
#' by taking the inverse of the negative of the
#' Hessian. Alternatively, \link{boot.ascr} can be used to carry
#' out a parametric bootstrap procedure, from which parameter
#' uncertainty can also be inferred.
#'
#' If the data are from an acoustic survey where individuals call more
#' than once (i.e., the argument \code{cue.rates} contains values
#' that are not 1), then standard errors calculated from the inverse
#' of the negative Hessian are not correct. They are therefore not
#' provided in this case. The method used by the function
#' \link{boot.ascr} is currently the only way to calculate these
#' reliably (see Stevenson et al., 2015, for details)
#'
#' @section The \code{capt} argument:
#'
#' The \code{capt} argument is a list with named components. Each
#' component must be an \eqn{n} by \eqn{k} matrix, where \eqn{n} is
#' the number of detections made, and \eqn{k} is the number of traps
#' (or detectors) deployed. A component named \code{bincapt} is
#' compulsory.
#'
#' Further optional component names each refer to a type of additional
#' information that is informative on animal location collected from
#' each detection. Possible choices are: \code{bearing}, \code{dist},
#' \code{ss}, \code{toa}, and \code{mrds}.
#'
#' If the \eqn{i}th individual evaded the \eqn{j}th trap (or
#' detector), then the \eqn{j}th element in the \eqn{i}th row should
#' be 0 for all components. Otherwise, if the \eqn{i}th individual was
#' trapped (or detected) by the \eqn{j}th trap (or detector), then:
#' \itemize{
#'   \item For the \code{bincapt} component, the element should be 1.
#'   \item For the \code{bearing} component, the element should be the
#'         estimated bearing from which the detector detected the
#'         individual.
#'   \item For the \code{dist} component, the element should be the
#'         estimated distance between the individual and the detector
#'         at the time of the detection.
#'   \item For the \code{ss} component, the element should be the
#'         measured signal strength of an acoustic signal detected by
#'         the detector (only possible when the detectors are
#'         microphones).
#'   \item For the \code{toa} component, the element should be the
#'         measured time of arrival (in seconds) since the start of
#'         the survey (or some other reference time) of an acoustic
#'         signal detected by the detector (only possible when the
#'         detectors are microphones).
#'   \item For the \code{mrds} component, the element should be the
#'         \emph{known} (not estimated) distance between the individual
#'         and the detector at the time of the detection.
#' }
#'
#' @section The \code{ss.opts} argument:
#'
#' This argument allows the user to select options for the signal
#' strength detection function (for more details, see the section
#' below on fitted parameters). It is therefore only required if
#' signal strength information appears in the \code{capt} argument,
#' and is ignored (with a warning) otherwise.
#'
#' The argument \code{ss.opts} is a list with up to seven components:
#' \itemize{
#'   \item \code{cutoff}: Compulsory. The signal strength threshold,
#'         above which sounds are identified as detections.
#'   \item \code{lower.cutoff}: Optional. Used for models where only
#'         the first detected call is used in the capture history. The lower
#'         cutoff is the signal strength value above which calls can be
#'         assumed to have been detected with certainty.
#'   \item \code{het.source}: Optional. Logical, if \code{TRUE} a model with
#'         heterogeneity in source signal strengths is used. If unspecified,
#'         it will default to \code{FALSE}.
#'   \item \code{het.source.method}: Optional. A character string, either
#'         \code{"GH"} or \code{"rect"}. If "GH", integration over source strengths
#'         uses Gauss-Hermite quadrature. If "rect", the rectangle method is used.
#'   \item \code{n.het.source.quadpoints}: Optional. An integer, giving the number
#'         of quadrature points used for numerical integration over source strengths.
#'         Defaults to 15. A larger number of quadrature points leads to more accurate
#'         results, but will increase computation time.
#'   \item \code{directional}: Optional. Logical, if \code{TRUE} a
#'         directional signal strength model is used; see the section below on
#'         fitted parameters. If unspecified, it will default to \code{FALSE},
#'         unless the \code{b2.ss} parameter is provided in \code{sv} or
#'         \code{fix}, in which case it will default to \code{TRUE}.
#'   \item \code{n.dir.quadpoints}: Optional. An integer, giving the number of
#'         quadrature points used for numerical integration over the possible
#'         call directions. Defaults to 8, but needs to be larger when calls are
#'         more directional (i.e., b2.ss parameter is large). A larger number of
#'         quadrature points leads to more accurate results, but will increase computation
#'         time.
#'   \item \code{ss.link}: Optional. A character string, either
#'         \code{"identity"}, \code{"log"}, or \code{"spherical"}, which
#'         specifies the relationship between the expected received signal
#'         strength and distance from the microphone. See details on the
#'         signal strength detection function in the section 'Fitted
#'         parameters' below. Defaults to \code{"identity"}.
#'
#'
#' }
#'
#' @section The \code{ihd.opts} argument:
#'
#' This argument allows the user to select options for the fitting on
#' inhomogeneous density surfaces.
#'
#' The argument \code{ihd.opts} is a list with up to three components:
#' \itemize{
#'    \item \code{model}: Compulsory. An equation for the relationship between
#'          covariates and the log of the density surface.
#'    \item \code{covariates}: Compulsory. A list of data frames, one
#'          for each session. Each data frame provides covariate
#'          values at each mask point.
#'    \item \code{scale}: Optional. If \code{TRUE}, the default,
#'          covariates are scaled by subtracting the mean and dividing
#'          by the standard deviaton. This does not affect model
#'          inference and improves optimisation stability, but makes
#'          it more difficult to interpret estimated coefficients.
#' }
#' @section The \code{optim.opts} argument:
#'
#' This argument allows the user to select options for the
#' maximisation of the likelihood.
#'
#' The argument \code{optim.opts} is a list with up to four components:
#' \itemize{
#'
#'   \item \code{cbs}: Optional. The CMPDIF_BUFFER_SIZE, set using the
#' \code{-cbs} option of the executable created by ADMB. This can be
#' increased to speed up optimisation if \code{cmpdiff.tmp} gets too
#' large (please ignore, unless you are familiar with ADMB and know
#' what you are doing).
#'
#'   \item \code{gbs}: Optional. The GRADSTACK_BUFFER_SIZE, set using
#' the \code{-gbs} option of the executable created by ADMB. This can
#' be increased to speed up optimisation if \code{gradfil1.tmp} gets
#' too large (please ignore, unless you are familiar with ADMB and
#' know what you are doing).
#'
#'   \item \code{exe.type}: Optional. Character string, either
#' \code{"old"} or \code{"new"}, depending on which executable is to
#' be used (for development purposes only; please ignore).
#'
#'   \item \code{neld.mead}: Optional. A logical value specifying
#' whether or not to use Nelder-Mead optimisation. Defaults to
#' \code{FALSE}, which is recommended.
#'
#' }
#'
#' @section Fitted parameters:
#'
#' The parameter \code{D}, the density of individuals (or, in an
#' acoustic survey, the density of calls) is always fitted. The
#' effective survey area, \code{esa}, (see Borchers, 2012, for
#' details) is always provided as a derived parameter, with a standard
#' error calculated using the delta method.
#'
#' Further parameters to be fitted depend on the choice of the
#' detection function (i.e., the \code{detfn} argument), and the types
#' of additional information collected (i.e., the components in the
#' \code{capt}).
#'
#' Details of the detection functions are as follows:
#'
#' For \code{detfn = "hn"}:
#' \itemize{
#'    \item Estimated paramters are \code{g0} and \code{sigma}.
#'    \item \eqn{g(d) = g_0\ exp(-d^2/(2\sigma^2))}{g(d) = g0 * exp( -d^2 / (2 * sigma^2 ))}
#' }
#'
#' For \code{detfn = "hr"}:
#' \itemize{
#'    \item Estimated parameters are \code{g0}, \code{sigma}, and
#'          \code{z}.
#'    \item \eqn{g(d) = g_0\ (1 - exp(-(d/\sigma)^{-z}))}{g(d) = g0 * ( 1 - exp( -(d/sigma)^{-z} ) )}
#' }
#'
#' For \code{detfn = "lth"}:
#' \itemize{
#'   \item Estimated parameters are \code{shape.1}
#'         \ifelse{latex}{(\eqn{\kappa})}{}, \code{shape.2}
#'         \ifelse{latex}{(\eqn{\nu})}{}, and \code{scale}
#'         \ifelse{latex}{(\eqn{\tau})}{}.
#'   \item \eqn{g(d) = 0.5 - 0.5\ erf(\kappa - exp(\nu - \tau d))}{g(d) = 0.5 - 0.5 * erf( shape.1 - exp( shape.2 - scale * d ) )}
#' }
#'
#' For \code{detfn = "th"}:
#' \itemize{
#'   \item Estimated parameters are \code{shape}
#'         \ifelse{latex}{(\eqn{\kappa})}{} and \code{scale}
#'         \ifelse{latex}{(\eqn{\tau})}{}.
#'   \item \eqn{g(d) = 0.5 - 0.5\ erf(d/\tau - \kappa)}{g(d) = 0.5 - 0.5 * erf( d/scale - shape )}
#' }
#'
#' For \code{detfn = "ss"} in a non-directional model:
#' \itemize{
#'   \item The signal strength detection function is special in that
#'         it requires signal strength information to be collected in
#'         order for all parameters to be estimated.
#'   \item Estimated parameters are \code{b0.ss}, \code{b1.ss}, and
#'         \code{sigma.ss}.
#'   \item The expected signal strength is modelled as:
#'         \eqn{E(SS) = h^{-1}(\beta_0 - \beta_1d)}{E(SS) = h^{-1}(b0.ss - b1.ss*d)},
#'         where \eqn{h} is specified by the argument \code{ss.link}.
#' }
#'
#' For \code{detfn = "ss"} in a directional model:
#' \itemize{
#'   \item Estimated paramters are \code{b0.ss}, \code{b1.ss}, \code{b2.ss} and
#'         \code{sigma.ss}.
#'   \item The expected signal strength is modelled differently depending on the value of \code{ss.link} in \code{ss.opts}:
#'   \itemize{
#'     \item For \code{ss.link = "identity"} (the default):
#'     \itemize{
#'       \item \eqn{E(SS) = \beta_0 - (\beta_1 - (\beta_2(\cos(\theta) - 1)))d)}{E(SS) = h^{-1}( b0.ss - ( b1.ss - ( b2.ss * ( cos( theta ) - 1 ) ) ) * d }
#'     }
#'     \item For \code{ss.link = "log"}:
#'     \itemize{
#'       \item \eqn{E(SS) = log(\beta_0 - (\beta_1 - (\beta_2(\cos(\theta) - 1)))d)}{E(SS) = h^{-1}( b0.ss - ( b1.ss - ( b2.ss * ( cos( theta ) - 1 ) ) ) * d ) }
#'     }
#'     \item For \code{ss.link = "spherical"}:
#'     \itemize{
#'       \item \eqn{E(SS) = \beta_0 - 10\log_{10}(d^2) - ( \beta_1 - ( \beta_2(\cos(\theta ) - 1)))(d - 1)}{E(SS) = \beta_0 - 10 * \log_{10}(d^2) - ( b1.ss - ( b2.ss( \cos( \theta ) - 1 ) ) ) * ( d - 1 )}
#'     }
#'   }
#'   \item In all cases \eqn{\theta}{theta} is the difference between
#'   the bearing the animal is facing when it makes a call, and the
#'   bearing from the animal to the detector.
#'
#' }
#'
#' Details of the parameters associated with different additional data
#' types are as follows:
#'
#' For data type \code{"bearing"}, \code{kappa} is estimated. This is
#' the concerntration parameter of the von-Mises distribution used for
#' measurement error in estimated bearings.
#'
#' For data type \code{"dist"}, \code{alpha} is estimated. This is the
#' shape parameter of the gamma distribution used for measurement
#' error in estimated distances.
#'
#' For data type \code{"toa"}, \code{sigma.toa} is estimated. This is
#' the standard deviation parameter of the normal distribution used
#' for measurement error in recorded times of arrival.
#'
#' For data type \code{"mrds"}, no extra parameters are
#' estimated. Animal location is assumed to be known.
#'
#' @section Convergence:
#'
#' If maximum likelihood estimates could not be found during
#' optimisation, then \code{fit.ascr()} will usually show a warning that
#' the maximum gradient component is large, or possibly throw an error
#' reporting that a \code{.par} file is missing.
#'
#' The best approach to fixing convergence issues is to re-run the
#' \code{fit.ascr()} function with the argument \code{trace} set to
#' \code{TRUE}. Parameter values will be printed out for each step of
#' the optimisation algorithm.
#'
#' First, look for a large jump in a parameter to a value far from
#' what is feasible. This issue can be fixed by using the
#' \code{bounds} argument to restrict the parameter space over which
#' ADMB searches for the maximum likelihood estimate.
#'
#' Alternatively, try a different set of start values using the
#' argument \code{sv}; by default \code{fit.ascr()} will choose some
#' start values, but these are not necessarily sensible. The start
#' values that were used appear as the first line of text when
#' \code{trace} is \code{TRUE}.
#'
#' Sometimes the algorithm appears to converge, but nevertheless
#' perseveres reporting the same parameter values again and again for
#' a while (prior to the calculation of the Hessian). This is because
#' ADMB has failed to detect convergence as at least one of the
#' gradient components is still larger than the convergence criterion
#' (by default, 0.0001). It is possible to speed things up and help
#' ADMB detect convergence earlier by tightening parameter bounds (as
#' above), by setting parameter phases, or by setting appropriate
#' scalefactors.
#'
#' To improve convergence using parameter phases use the \code{phases}
#' argument. By default all parameters are given a phase of 1, unless
#' it is changed. First, all parameters with a phase of 1 will be
#' maximised over, while all others are fixed at their original
#' values. Following that, parameters with a phase of 2 are
#' introduced, with all parameters with later phases remaining
#' fixed. This process continues until all parameters are maximised
#' over. Maximising paramters in phases can greatly improve the
#' stability of optimisation.
#'
#' To improve convergence using scalefactors, first identify which
#' parameters have large gradient components from the "final
#' statistics" section of the \code{trace} output. Next, find the
#' default settings of the scalefactors by printing the object
#' \code{fit$args$sf}, where \code{fit} is the original object
#' returned by \code{fit.ascr()}. Finally, rerun \code{fit.ascr()} again,
#' but this time set the argument \code{sf} manually. Set scalefactors
#' for any parameters with small gradient components to the same as
#' the defaults ascertained above, and increase those associated with
#' large gradient components by a factor of 10. If the problem
#' persists, repeat this process (e.g., if the same parameters still
#' have large gradient components, increase the associated
#' scalefactors by another factor of 10).
#'
#' @section Local integration:
#'
#' For SECR models, the likelihood is calculated by integrating over
#' the unobserved animal activity centres (see Borchers & Efford,
#' 2008). Here, the integral is approximated numerically by taking a
#' finite sum over the mask points. The integrand is negligible in
#' size for mask points far from detectors that detected a particular
#' individual, and so to increase computational efficiency the region
#' over which this sum takes place can be reduced.
#'
#' Setting \code{local} to \code{TRUE} will only carry out this sum
#' across mask points that are within the mask buffer distance of
#' \emph{all} detectors that made a detection. So long as the buffer
#' suitably represents a distance beyond which detection is
#' practically impossible, the effect this has on parameter estimates
#' is negligible, but processing time can be substantially reduced.
#'
#' Note that this increases the parameter estimates' sensitivity to
#' the buffer. A buffer that is too small will lead to inaccurate
#' results.
#'
#' @references Borchers, D. L., and Efford, M. G. (2008) Spatially
#'     explicit maximum likelihood methods for capture-recapture
#'     studies. \emph{Biometrics}, \strong{64}: 377--385.
#'
#' @references Borchers, D. L. (2012) A non-technical overview of
#'     spatially explicit capture-recapture models. \emph{Journal of
#'     Ornithology}, \strong{152}: 435--444.
#'
#' @references Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas,
#'     L., and Marques, T. A. (2015) A unifying model for
#'     capture-recapture and distance sampling surveys of wildlife
#'     populations. \emph{Journal of the American Statistical
#'     Association}, \strong{110}: 195--204.
#'
#' @references Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift,
#'     R. J., Gillespie, D. M., and Measey, G. J. (2015) A general
#'     framework for animal density estimation from acoustic
#'     detections across a fixed microphone array. \emph{Methods in
#'     Ecology and Evolution}, \strong{6}: 38--48.
#'
#' @return A list of class \code{"ascr"}. Components contain
#'     information such as estimated parameters and standard
#'     errors. The best way to access such information, however, is
#'     through the variety of helper functions provided by the
#'     ascr package.
#'
#' @param capt A list with named components, containing the capture
#'     history and supplementary information. The function
#'     \link{create.capt} will return a suitable object. See 'Details'
#'     below. Alternatively, this can be a list of such lists if
#'     detections from multiple detector arrays are being used to fit
#'     a single model.
#' @param traps A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a trap (or
#'     detector). Alternatively, this can be a list of such matrices
#'     if detections from multiple detector arrays are being used to
#'     fit a single model.
#' @param mask A matrix with two columns. Each row provides Cartesian
#'     coordinates for the location of a mask point. The function
#'     \link{create.mask} will return a suitable
#'     object. Alternatively, this can be a list of such matrices if
#'     detections from multiple detector arrays are being used to fit
#'     a single model.
#' @param detfn A character string specifying the detection function
#'     to be used. One of "hn" (halfnormal), "hr" (hazard rate), "th"
#'     (threshold), "lth" (log-link threshold), or "ss" (signal
#'     strength). If the latter is used, signal strength information
#'     must be provided in \code{capt}.
#' @param sv A named list. Component names are parameter names, and
#'     each component is a start value for the associated
#'     parameter. See 'Details' for further information on the
#'     parameters to be fitted.
#' @param bounds A named list. Component names are parameter names,
#'     and each components is a vector of length two, specifying the
#'     bounds for the associated parameter.
#' @param fix A named list. Component names are parameter names to be
#'     fixed, and each component is the fixed value for the associated
#'     parameter.
#' @param phases A named list. Component names are parameter names,
#'     and each component is a phase for the associated parameter. See
#'     the section on convergence below for information on parameter
#'     phases.
#' @param sf A named list. Component names are parameter names, and
#'     each component is a scalefactor for the associated
#'     parameter. The default behaviour is to automatically select
#'     scalefactors based on parameter start values. See the section
#'     on convergence below.
#' @param ss.opts Options for models using the signal strength
#'     detection function. See 'Details' below.
#' @param cue.rates A vector of call rates collected independently of
#'     the main acoustic survey. This must be measured in calls per
#'     unit time, where the time units are equivalent to those used by
#'     \code{survey.length}.
#' @param survey.length The length of a cue-based survey. If provided,
#'     the estimated density \code{Dc} is measured in cues per unit
#'     time (using the same units as \code{survey.length}). For
#'     multi-session data, this must be a vector, giving the survey
#'     lengths for each session.
#' @param sound.speed The speed of sound in metres per second,
#'     defaults to 330 (the speed of sound in air). Only used when
#'     \code{"toa"} is a component name of \code{capt}.
#' @param local Logical, if \code{TRUE} integration over unobserved
#'     animal activity centres is only carried out in a region local
#'     to detectors that detected individuals. See 'Details'.
#' @param ihd.opts Options for inhomogeneous density. See 'Details'
#'     below.
#' @param hess Logical, if \code{TRUE} the Hessian is estimated,
#'     allowing for calculation of standard errors, the
#'     variance-covariance matrix, and the correlation matrix, at the
#'     expense of a little processing time. If \code{FALSE}, the
#'     Hessian is not estimated. Note that if individuals are
#'     detectable more than once (e.g., by calling more than once on
#'     an acoustic survey) then parameter uncertainty is not properly
#'     represented by these calculations.
#' @param trace Logical, if \code{TRUE} parameter values at each step
#'     of the optimisation algorithm are printed to the R console.
#' @param clean Logical, if \code{TRUE} ADMB output files are
#'     removed. Otherwise, ADMB output file will remain in a
#'     directory, the location of which is reported after the model is
#'     fitted.
#' @param optim.opts Optimisation options. See 'Details' for further
#'     information.
#' @param ... Other arguments (mostly for back-compatibility).
#'
#' @seealso \link{boot.ascr} to calculate standard errors and
#'     estimate bias using a parametric bootstrap.
#' @seealso \link{coef.ascr}, \link{stdEr.ascr}, and
#'     \link{vcov.ascr} to extract estimated parameters, standard
#'     errors, and the variance-covariance matrix, respectively.
#' @seealso \link{confint.ascr} to calculate confidence intervals.
#' @seealso \link{summary.ascr} to get a summary of estimates and
#'     standard errors.
#' @seealso \link{show.detfn} to plot the estimated detection
#'     function.
#' @seealso \link{locations} to plot estimated locations of particular
#'     individuals or calls.
#'
#' @examples
#' \dontrun{
#' simple.capt <- example$capt["bincapt"]
#' simple.hn.fit <- fit.ascr(capt = simple.capt, traps = example$traps,
#'                           mask = example$mask, fix = list(g0 = 1))
#' simple.hr.fit <- fit.ascr(capt = simple.capt, traps = example$traps,
#'                           mask = example$mask, detfn = "hr")
#' bearing.capt <- example$capt[c("bincapt", "bearing")]
#' bearing.hn.fit <- fit.ascr(capt = bearing.capt, traps = example$traps,
#'                            mask = example$mask, fix = list(g0 = 1))
#' }
#'
#' @export
fit.ascr <- function(capt, traps, mask, detfn = "hn", sv = NULL, bounds = NULL,
                     fix = NULL, phases = NULL, sf = NULL, ss.opts = NULL,
                     cue.rates = NULL, survey.length = NULL, sound.speed = 330,
                     ihd.opts = NULL, local = FALSE, hess = NULL, trace = FALSE,
                     clean = TRUE, optim.opts = NULL, ...){
    arg.names <- names(as.list(environment()))
    extra.args <- list(...)
    ## Sorting out multi-session stuff.
    if (is.data.frame(traps)){
        traps <- as.matrix(traps)
    }
    if (is.data.frame(mask)){
        mask <- as.matrix(mask)
    }
    multi.session <- ifelse(is.list(traps) & is.list(mask), TRUE, FALSE)
    ## If only a single session, just make multi-session objects with a single component.
    if (!multi.session){
        if (is.list(traps) | is.list(mask)){
            stop("For multi-session models, both `traps' and `mask' must be lists.")
        }
        capt <- list(capt)
        traps <- list(traps)
        mask <- list(mask)
    }
    if (!is.list(capt[[1]])){
        capt <- list(capt)
    }
    if (length(traps) != length(mask)){
        stop("For multi-session models, both `traps' and `mask' must have the same number of components.")
    }
    n.sessions <- length(traps)
    if (length(capt) != n.sessions){
        stop("For multi-session models, `capt' must have a component for each session.")
    }
    capt.names <- names(capt[[1]])
    for (i in 1:length(capt)){
        if (!all(names(capt[[i]]) == capt.names)){
            stop("For multi-session models, all components of `capt' must have the same data types.")
        }
    }
    if (any(names(extra.args) == "call.freqs")){
        if (!missing(cue.rates)){
            stop("The argument `cue.rates' has replaced `call.freqs'; use only the former.")
        }
        warning("The argument `call.freqs' is deprecated; please rename to `cue.rates' instead.")
        cue.rates <- extra.args[["call.freqs"]]
    }
    ## Determining supplementary parameter names.
    supp.types <- c("bearing", "dist", "ss", "toa", "mrds")
    fit.types <- supp.types %in% capt.names
    names(fit.types) <- supp.types
    ## Logical indicators for additional information types.
    fit.bearings <- fit.types["bearing"]
    fit.dists <- fit.types["dist"]
    fit.ss <- fit.types["ss"]
    fit.toas <- fit.types["toa"]
    fit.mrds <- fit.types["mrds"]
    ## Warning from cue.rates without survey.length.
    if (is.null(survey.length)){
        survey.length <- rep(1, n.sessions)
        if (!is.null(cue.rates)){
            stop("The use of `cue.rates' without `survey.length' is no longer supported. Please provide `survey.length', and ensure `cue.rates' is measured in the same time units.")
        }
    } else {
        if (length(survey.length) != n.sessions){
            stop("The argument `survey.length' must have a value for each session.")
        }
    }
    ## Sorting out cues per survey.
    if (!is.null(cue.rates)){
        cue.freqs <- cue.rates*survey.length
    }
    ## Storing objects from ss.opts.
    cutoff <- ss.opts$cutoff
    ss.link <- ss.opts$ss.link
    directional <- ss.opts$directional
    het.source <- ss.opts$het.source
    het.source.method <- ss.opts$het.source.method
    n.dir.quadpoints <- ss.opts$n.dir.quadpoints
    n.het.source.quadpoints <- ss.opts$n.het.source.quadpoints
    lower.cutoff <- ss.opts$lower.cutoff
    ## Sorting objects from optim.opts.
    cbs <- optim.opts$cbs
    gbs <- optim.opts$gbs
    exe.type <- optim.opts$exe.type
    neld.mead <- optim.opts$neld.mead
    if (is.null(exe.type)){
        exe.type <- "old"
    }
    if (is.null(neld.mead)){
        neld.mead <- FALSE
        neld.mead.force <- FALSE
    } else {
        neld.mead.force <- TRUE
    }
    ## Setting up first.calls indicator.
    first.calls <- FALSE
    if (fit.ss){
        if (missing(ss.opts)){
            ## Error if ss.opts not provided for signal strength model.
            stop("Argument 'ss.opts' is missing.")
        }
        ## Error for unspecified cutoff.
        if (is.null(cutoff)){
            stop("The 'cutoff' component of 'ss.opts' must be specified.")
        }
        if (!is.null(lower.cutoff)){
            first.calls <- TRUE
            if (!(lower.cutoff < cutoff)){
                stop("The 'lower.cutoff' component of 'ss.opts' must be lower than the 'cutoff' component.")
            }
        } else {
            ss.opts["lower.cutoff"] <- list(NULL)
            lower.cutoff <- NULL
        }
        ## Removing detections below the cutoff.
        n.removed <- 0
        for (i in 1:n.sessions){
            rem <- capt[[i]]$ss < cutoff
            capt[[i]] <- lapply(capt[[i]], function(x, rem){
                x[rem] <- 0
                x
            }, rem = rem)
            keep <- apply(capt[[i]]$bincapt, 1, sum) > 0
            capt[[i]] <- lapply(capt[[i]], function(x, keep) x[keep, ], keep = keep)
            n.removed <- n.removed + sum(!keep)
        }
        if (trace & n.removed > 0){
            message(n.removed, " capture history entries have no received signal strengths above the cutoff and have therefore been removed.\n", sep = "")
        }
    }
    capt.bin <- vector(mode = "list", length = n.sessions)
    for (i in 1:n.sessions){
        ## Checking for bincapt.
        if (!any(names(capt[[i]]) == "bincapt")){
            stop("The binary capture history must be provided as a component of 'capt'.")
        }
        capt.bin[[i]] <- capt[[i]]$bincapt
        ## Checking for correct number of trap locations.
        if (ncol(capt.bin[[i]]) != nrow(traps[[i]])){
            stop("There must be a trap location for each column in the components of 'capt'.")
        }
        ## Checking that each component of 'capt' is a matrix.
        if (any(!laply(capt[[i]], is.matrix))){
            stop("At least one component of 'capt' is not a matrix.")
        }
        ## Checking for agreement in matrix dimensions.
        if (length(capt[[i]]) > 1){
            all.dims <- laply(capt[[i]], dim)
            if (any(aaply(all.dims, 2, function(x) diff(range(x))) != 0)){
                stop("Components of 'capt' object within a session have different dimensions.")
            }
        }
    }
    ## Various checks for other arguments.
    if (!is.list(sv) & !is.null(sv)){
        stop("The 'sv' argument must be 'NULL' or a list.")
    }
    if (!is.list(bounds) & !is.null(bounds)){
        stop("The 'bounds' argument must be 'NULL' or a list.")
    }
    if (!is.list(phases) & !is.null(phases)){
        stop("The 'phases' argument must be 'NULL' or a list.")
    }
    if (is.list(bounds)){
        if (any(laply(bounds, length) != 2)){
            stop("Each component of 'bounds' must be a vector of length 2.")
        }
    }
    if (!is.list(fix) & !is.null(fix)){
        stop("The 'fix' argument must be 'NULL' or a list.")
    }
    n <- sapply(capt.bin, nrow)
    n.traps <- sapply(traps, nrow)
    n.mask <- sapply(mask, nrow)
    A <- sapply(mask, function(x) attr(x, "area"))
    buffer <- sapply(mask, function(x) attr(x, "buffer"))
    ## Removing attributes from mask.
    for (i in 1:n.sessions){
        mask[[i]] <- as.matrix(mask[[i]])
        traps[[i]] <- as.matrix(traps[[i]])
        attr(mask[[i]], "area") <- A[i]
        attr(mask[[i]], "buffer") <- buffer[i]
    }
    ## Sorting out inhomogeneous density stuff.
    if (!is.null(ihd.opts)){
        if (ihd.opts$model != ~1){
            fit.ihd <- TRUE
        } else {
            fit.ihd <- FALSE
        }
    } else {
        ihd.opts$model <- ~ 1
        fit.ihd <- FALSE
    }
    if (is.data.frame(ihd.opts$covariates)){
        covariates <- list()
        covariates[[1]] <- data.frame(mask[[1]], ihd.opts$covariates)
    } else if (is.list(ihd.opts$covariates)){
        covariates <- list()
        for (i in 1:n.sessions){
            covariates[[i]] <- data.frame(mask[[i]], ihd.opts$covariates[[i]])
        }
    } else if (is.null(ihd.opts$covariates)){
        covariates <- mask
    }
    if (is.null(ihd.opts$scale)){
        cov.scale <- TRUE
    } else {
        cov.scale <- ihd.opts$scale
    }
    D.mask <- list()
    mm.ihd <- list()
    for (i in 1:n.sessions){
        if (cov.scale){
            covariates[[i]] <- as.data.frame(apply(covariates[[i]], 2, function(x) (x - mean(x))/sd(x)))
        }
        ## Extracting the formula.
        model.formula <- ihd.opts$model
        ## Need a response variable for gam() to work.
        model.formula <- as.formula(paste("rep(0, nrow(covariates[[i]]))", paste(as.character(model.formula), collapse="")))
        fgam <- gam(model.formula, data = covariates[[i]], fit = FALSE)
        mm.ihd[[i]] <- fgam$X
        colnames(mm.ihd[[i]]) <- fgam$term.names
    }
    D.betapars.names <- paste("D.", colnames(mm.ihd[[1]]), sep = "")
    ## Sorting out signal strength options.
    if (fit.ss){
        ## Warning for unexpected component names.
        if (!all(names(ss.opts) %in% c("cutoff", "het.source", "het.source.method", "n.het.source.quadpoints", "directional", "n.dir.quadpoints", "ss.link", "lower.cutoff"))){
            warning("Components of 'ss.opts' may only consist of \"cutoff\", \"het.source\", \"het.source.method\", \"n.het.source.quadpoints\", \"directional\",  \"n.dir.quadpoints\", \"ss.link\", and \"lower.cutoff\"; others are being ignored.")
        }
        ## Setting default values for ss.link, het.source and directional.
        if (is.null(ss.link)){
            ss.opts$ss.link <- "identity"
            ss.link <- "identity"
        } else if (!(ss.link %in% c("identity", "log", "spherical"))){
            stop("Component 'ss.link' in 'ss.opts' must be \"identity\", \"log\", or \"spherical\".")
        }
        if (first.calls & ss.link != "identity"){
            stop("First-call models are only implemented for ss.link = \"identity\".")
        }
        ## By default, directional calling model is only used if b2.ss appears in sv or fix.
        if (is.null(directional)){
            if (is.null(sv$b2.ss) & is.null(fix$b2.ss)){
                ss.opts$directional <- FALSE
                directional <- FALSE
            } else {
                ss.opts$directional <- TRUE
                directional <- TRUE
            }
        }
        ## Fixing b2.ss to 0 if a directional calling model is not being used.
        if (!directional){
            warn.directional <- FALSE
            if (!is.null(sv$b2.ss)){
                if (sv$b2.ss != 0){
                    warn.directional <- TRUE
                }
                sv$b2.ss <- NULL
            }
            if (!is.null(fix$b2.ss)){
                if (fix$b2.ss != 0){
                    warn.directional <- TRUE
                }
                fix$b2.ss <- NULL
            }
            if (warn.directional){
                warning("As the 'directional' component of 'ss.opts' is FALSE, the values of parameter b2.ss in 'sv' and 'fix' are being ignored")
            }
            fix$b2.ss <- 0
        }
        ## By default, heterogeneity source strength model is only
        ## used if sigma.b0.ss appears in sv or fix.
        if (is.null(het.source)){
            if (is.null(sv$sigma.b0.ss) & is.null(fix$sigma.b0.ss)){
                ss.opts$het.source <- FALSE
                het.source <- FALSE
            } else {
                ss.opts$het.source <- TRUE
                het.source <- TRUE
                ss.opts$het.source.method <- "GH"
                het.source.method <- "GH"
            }
        }
        if (het.source){
            if (is.null(het.source.method)){
                ss.opts$het.source.method <- "GH"
                het.source.method <- "GH"
            }
        } else {
            ## Fixing sigma.b0.ss to 0 if a heterogeneous source
            ## strength model is not being used.
            warn.het <- FALSE
            if (!is.null(sv$sigma.b0.ss)){
                if (sv$sigma.b0.ss != 0){
                    warn.het <- TRUE
                }
                sv$sigma.b0.ss <- NULL
            }
            if (!is.null(fix$sigma.b0.ss)){
                if (fix$sigma.b0.ss != 0){
                    warn.het <- TRUE
                }
                fix$sigma.b0.ss <- NULL
            }
            if (warn.het){
                warning("As the 'het.source' component of 'ss.opts' is FALSE, the values of the parameter sigma.b0.ss in 'sv' and 'fix' are being ignored")      
            }
            fix$sigma.b0.ss <- 0
        }
    } else {
        if (!is.null(ss.opts)){
            warning("Argument 'ss.opts' is being ignored as a signal strength model is not being fitted.")
        }
        ss.opts <- NULL
        lower.cutoff <- NULL
    }
    ## Setting fit.dir.
    if (fit.ss){
        fit.dir <- TRUE
        if ("b2.ss" %in% names(fix)){
            if (fix[["b2.ss"]] == 0){
                fit.dir <- FALSE
            }
        }
    } else {
        fit.dir <- FALSE
    }
    ## Setting fit.het.source.
    if (fit.ss){
        fit.het.source <- TRUE
        if ("sigma.b0.ss" %in% names(fix)){
            if (fix[["sigma.b0.ss"]] == 0){
                fit.het.source <- FALSE
            }
        }
    } else {
        fit.het.source <- FALSE
    }
    ## Setting het.source.gh.
    if (fit.het.source){
        het.source.gh <- het.source.method == "GH"
    } else {
        het.source.gh <- FALSE
    }
    if (fit.het.source & first.calls){
        stop("Models with both first calls and heterogeneity in source signal strengths are not yet implemented.")
    }
    if (fit.dir & first.calls){
        stop("Models with both first calls and directional calling are not yet implemented.")
    }
    ## Supplementary parameter names.
    suppar.names <- c("kappa", "alpha", "sigma.toa")[fit.types[c("bearing", "dist", "toa")]]
    ## Generating ordered binary capture history.
    capt.bin.order <- vector(mode = "list", length = n.sessions)
    capt.bin.unique <- vector(mode = "list", length = n.sessions)
    capt.bin.freqs <- vector(mode = "list", length = n.sessions)
    n.unique <- numeric(n.sessions)
    capt.bearing <- vector(mode = "list", length = n.sessions)
    capt.dist <- vector(mode = "list", length = n.sessions)
    capt.ss <- vector(mode = "list", length = n.sessions)
    capt.toa <- vector(mode = "list", length = n.sessions)
    mrds.dist <- vector(mode = "list", length = n.sessions)
    capt.ord <- vector(mode = "list", length = n.sessions)
    for (i in 1:n.sessions){
        capt.bin.order[[i]] <- do.call(order, as.data.frame(capt.bin[[i]]))
        capt.bin.unique[[i]] <- capt.bin[[i]][capt.bin.order[[i]], , drop = FALSE]
        capt.bin.freqs[[i]] <- as.vector(table(apply(capt.bin.unique[[i]], 1, paste, collapse = "")))
        names(capt.bin.freqs[[i]]) <- NULL
        capt.bin.unique[[i]] <- capt.bin.unique[[i]][!duplicated(as.data.frame(capt.bin.unique[[i]])), , drop = FALSE]
        n.unique[i] <- nrow(capt.bin.unique[[i]])
        ## Reordering all capture history components.
        capt.ord[[i]] <- capt[[i]]
        for (j in 1:length(capt[[i]])){
            capt.ord[[i]][[j]] <- capt[[i]][[j]][capt.bin.order[[i]], , drop = FALSE]
        }
        ## Capture histories for additional information types (if they exist)
        capt.bearing[[i]] <- if (fit.bearings) capt.ord[[i]]$bearing else 0
        capt.dist[[i]] <- if (fit.dists) capt.ord[[i]]$dist else 0
        capt.ss[[i]] <- if (fit.ss) capt.ord[[i]]$ss else 0
        capt.toa[[i]] <- if (fit.toas) capt.ord[[i]]$toa else 0
        mrds.dist[[i]] <- if (fit.mrds) capt.ord[[i]]$mrds else 0
        ## Data check for bearings.
        if (any(capt.bearing[[i]] < 0 | capt.bearing[[i]] > 2*pi)){
            warning("Some estimated bearings are not in the interval [0, 2*pi)")
        }
    }
    if (fit.ss){
        if (!missing(detfn) & detfn != "ss"){
            warning("Argument 'detfn' is being ignored as signal strength information is provided in 'capt'. A signal strength detection function has been fitted instead.")
        }
        if (ss.link == "identity"){
            detfn <- "ss"
            linkfn.id <- 1
        } else if (ss.link == "log"){
            detfn <- "log.ss"
            linkfn.id <- 2
        } else if (ss.link == "spherical"){
            detfn <- "spherical.ss"
            linkfn.id <- 3
        }
    } else {
        ## Not sure what a linkfn.id of 4 means? Probably throws an error in ADMB.
        linkfn.id <- 4
    }
    detfns <- c("hn", "hr", "th", "lth", "ss", "log.ss", "spherical.ss")
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
                           ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"),
                           log.ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"),
                           spherical.ss = c("b0.ss", "b1.ss", "b2.ss", "sigma.b0.ss", "sigma.ss"))
    par.names <- c(D.betapars.names, detpar.names, suppar.names)
    n.detpars <- length(detpar.names)
    n.suppars <- length(suppar.names)
    n.D.betapars <- length(D.betapars.names)
    any.suppars <- n.suppars > 0
    n.pars <- length(par.names)
    ## Checking par.names against names of sv, fix, bounds, and sf.
    for (i in c("sv", "fix", "bounds", "phases", "sf")){
        obj <- get(i)
        if (!is.null(obj)){
            if (i == "sv"){
                obj.fitted <- names(obj) %in% c("D", par.names)
            } else {
                obj.fitted <- names(obj) %in% c("D", par.names)
            }
            if(!all(obj.fitted)){
                msg <- paste("Some parameters listed in '", i, "' are not being used. These are being removed.",
                             sep = "")
                warning(msg)
                assign(i, obj[obj.fitted])
            }
        }
    }
    ## Sets link function ID number for use in ADMB:
    ## 1 = identity
    ## 2 = log
    ## 3 = logit
        parlinks.list <- list(g0 = 3,
                              sigma = 2,
                              shape = 1,
                              shape.1 = 2,
                              shape.2 = 1,
                              scale = 2,
                              b0.ss = 2,
                              b1.ss = 2,
                              b2.ss = 2,
                              sigma.b0.ss = 2,
                              sigma.ss = 2,
                              z = 2,
                              sigma.toa = 2,
                              kappa = 2,
                              alpha = 2)
    for (i in D.betapars.names){
        parlinks.list[i] <- 1
    }
    links <- parlinks.list[par.names]
    link.list <- list(identity, log.link, logit.link)
    unlink.list <- list(identity, exp, inv.logit)
    par.links <- llply(links, function(x, link.list) link.list[[x]], link.list)
    par.unlinks <- llply(links, function(x, unlink.list) unlink.list[[x]], unlink.list)
    ## Sorting out start values. Start values are set to those provided,
    ## or else are determined automatically from functions in
    ## autofuns.r.
    sv.link <- vector("list", length = n.pars)
    names(sv.link) <- par.names
    if (any(names(sv) == "D")){
        sv[["D"]] <- log(sv[["D"]])
        names(sv)[names(sv) == "D"] <- "D.(Intercept)"
    }
    autoD.generic <- function(args) 0
    for (i in D.betapars.names){
        if (i != "D.(Intercept)"){
            assign(paste0("auto", i), autoD.generic)
        }
    }
    sv.link[names(sv)] <- sv
    sv.link[names(fix)] <- fix
    auto.names <- par.names[sapply(sv.link, is.null)]
    sv.funs <- paste("auto", auto.names, sep = "")
    ## Use the first session with detections to generate start values.
    sess.to.use <- which(n > 0)[1]
    same.traplocs <-  all(distances(traps[[sess.to.use]], traps[[sess.to.use]]) == 0)
    ## Done in reverse so that D is calculated last (requires detfn
    ## parameters). D not moved to front as it should appear as the
    ## first parameter in any output.  Note that start value
    ## parameters are automatically generated only from the first
    ## session's data.
    for (i in rev(seq(1, length(auto.names), length.out = length(auto.names)))){
        sv.link[auto.names[i]] <- eval(call(sv.funs[i],
                                       list(capt = capt[[sess.to.use]], detfn = detfn,
                                            detpar.names = detpar.names,
                                            mask = mask[[sess.to.use]], traps = traps[[sess.to.use]],
                                            sv = sv.link, ss.opts = ss.opts,
                                            A = A[1], survey.length = survey.length[1],
                                            same.traplocs = same.traplocs)))
    }
    ## Converting start values to link scale.
    sv <- sv.link
    for (i in names(sv.link)){
        sv.link[[i]] <- link.list[[links[[i]]]](sv.link[[i]])
    }
    ## Sorting out phases.
    phases.save <- phases
    phases <- vector("list", length = n.pars)
    names(phases) <- par.names
    for (i in par.names){
        if (any(i == names(fix))){
            ## Phase of -1 in ADMB fixes parameter at starting value.
            phases[[i]] <- -1
        } else {
            if (any(i == names(phases.save))){
                phases[[i]] <- phases.save[[i]]
            } else {
                phases[[i]] <- 1
            }
        }
    }
    D.betapars.phase <- c(phases[D.betapars.names], recursive = TRUE)
    detpars.phase <- c(phases[detpar.names], recursive = TRUE)
    if (any.suppars){
        suppars.phase <- c(phases[suppar.names], recursive = TRUE)
    } else {
        suppars.phase <- -1
    }
    ## Sorting out bounds.
    ## Below bounds are the defaults.
    default.bounds.list <- list(g0 = c(0, 1),
                                sigma = c(0, 1e8),
                                shape = c(-100, 100),
                                shape.1 = c(0, 1e8),
                                shape.2 = c(-100, 100),
                                scale = c(0, 1e8),
                                b0.ss = c(0, 1e8),
                                b1.ss = c(0, 1e8),
                                b2.ss = c(0, 1e8),
                                sigma.b0.ss = c(0, 1e8),
                                sigma.ss = c(0, 1e8),
                                z = c(0, 1e8),
                                sigma.toa = c(0, 1e8),
                                kappa = c(0, 700),
                                alpha = c(0, 1e8))
    for (i in D.betapars.names){
        default.bounds.list[[i]] <- c(-log(1e20), log(1e20))
    }
    default.bounds <- default.bounds.list[par.names]
    if (any(names(bounds) == "D")){
        if (bounds[["D"]][1] == 0){
            bounds[["D"]][1] <- 1e-20
        }
        bounds[["D"]] <- log(bounds[["D"]])
        names(bounds)[names(bounds) == "D"] <- "D.(Intercept)"
    }
    bound.changes <- bounds
    bounds <- default.bounds
    for (i in names(default.bounds)){
        if (i %in% names(bound.changes)){
            bounds[[i]] <- bound.changes[[i]]
        }
    }
    ## Converting bounds to link scale.
    bounds.link <- bounds
    for (i in names(bounds)){
        bounds.link[[i]] <- link.list[[links[[i]]]](bounds[[i]])
    }
    D.betapars.bounds <- bounds.link[D.betapars.names]
    D.betapars.lb <- sapply(D.betapars.bounds, function(x) x[1])
    D.betapars.ub <- sapply(D.betapars.bounds, function(x) x[2])
    detpar.bounds <- bounds.link[detpar.names]
    detpars.lb <- sapply(detpar.bounds, function(x) x[1])
    detpars.ub <- sapply(detpar.bounds, function(x) x[2])
    if (any.suppars){
        suppar.bounds <- bounds.link[suppar.names]
        suppars.lb <- sapply(suppar.bounds, function(x) x[1])
        suppars.ub <- sapply(suppar.bounds, function(x) x[2])
    } else {
        suppars.lb <- 0
        suppars.ub <- 0
    }
    ## Sorting out scalefactors.
    if (is.null(sf)){
        sv.vec <- c(sv.link, recursive = TRUE)
        ## Currently, by default, the scalefactors are the inverse
        ## fraction of each starting value to the largest starting
        ## value. Not sure how sensible this is.
        ##bound.ranges <- laply(bounds.link, function(x) diff(range(x)))
        ##sf <- max(bound.ranges)/bound.ranges
        sf <- abs(max(sv.vec)/sv.vec)
        names(sf) <- par.names
    } else if (is.list(sf)){
        sf <- numeric(n.pars)
        names(sf) <- par.names
        for (i in par.names){
            sf[i] <- ifelse(i %in% names(sf), sf[[i]], 1)
        }
    } else if (is.vector(sf) & length(sf) == 1){
        sf <- rep(sf, length(par.names))
        names(sf) <- par.names
    }
    ## Replacing infinite scalefactors.
    sf[!is.finite(sf)] <- 1
    D.betapars.sf <- c(sf[D.betapars.names], recursive = TRUE)
    detpars.sf <- c(sf[detpar.names], recursive = TRUE)
    if (any.suppars){
        suppars.sf <- c(sf[suppar.names], recursive = TRUE)
    } else {
        suppars.sf <- 1
    }
    ## Creating link objects to pass to ADMB.
    detpars.link <- c(links[detpar.names], recursive = TRUE)
    if (any.suppars){
        suppars.link <- c(links[suppar.names], recursive = TRUE)
    } else {
        suppars.link <- 1
    }
    ## Setting small number so that numerical under/overflow in ADMB
    ## does not affect estimation.
    dbl.min <- 1e-150
    ## Calculating distances and angles.
    dists <- vector(mode = "list", length = n.sessions)
    bearings <- vector(mode = "list", length = n.sessions)
    toa.ssq <- vector(mode = "list", length = n.sessions)
    for (i in 1:n.sessions){
        dists[[i]] <- distances(traps[[i]], mask[[i]])
        if (fit.bearings | fit.dir){
            bearings[[i]] <- bearings(traps[[i]], mask[[i]])
        } else {
            bearings[[i]] <- 0
        }
        if (fit.toas){
            toa.ssq[[i]] <- make_toa_ssq(capt.ord[[i]]$toa, dists[[i]], sound.speed)
        } else {
            toa.ssq <- 0
        }
    }
    if (is.null(cutoff)){
        cutoff <- 0
    }
    ## Kludge to fix number of parameters for no supplementary
    ## information.
    if (!any.suppars){
        n.suppars <- max(c(n.suppars, 1))
        sv.link$dummy <- 0
        suppar.names <- "dummy"
    }
    ## Reordering sv.link. Matters when there is no supplementary information.
    sv.link <- sv.link[c(D.betapars.names, detpar.names, suppar.names)]
    ## ... And resetting suppar.names to NULL if required.
    if (!any.suppars){
        suppar.names <- NULL
    }
    ## Sorting out which mask points are local to each detection.
    all.which.local <- vector(mode = "list", length = n.sessions)
    all.n.local <- vector(mode = "list", length = n.sessions)
    if (local){
        for (i in 1:n.sessions){
            all.which.local[[i]] <- find_local(capt.bin.unique[[i]], dists[[i]], buffer[i])
            all.n.local[[i]] <- laply(all.which.local[[i]], length)
            all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
        }
    } else {
        for (i in 1:n.sessions){
            all.n.local[[i]] <- rep(1, n.unique[i])
            all.which.local[[i]] <- rep(0, n.unique[i])
        }
    }
    ## Sorting out number of quadrature points for directional calling.
    if (fit.dir){
        if (is.null(n.dir.quadpoints)){
            n.dir.quadpoints <- 8
        }
    } else {
        n.dir.quadpoints <- 1
    }
    ## Sorting out number of quadrature points for heteroegeous source strength.
    if (fit.het.source){
        if (is.null(n.het.source.quadpoints)){
            n.het.source.quadpoints <- 15
        }
    } else {
        n.het.source.quadpoints <- 1
    }
    ## Getting nodes and weights for Gauss-Hermite quadrature.
    if (het.source.gh){
        GHd <- gaussHermiteData(n.het.source.quadpoints)
        het.source.nodes <- GHd$x
        het.source.weights <- GHd$w
    } else {
        het.source.nodes <- 0
        het.source.weights <- 0
    }
    ## Error thrown for model with both heterogeneity in source strength and directional calling.
    if (fit.het.source & fit.dir){
        stop("Fitting of models with both heterogeneity in source signal strength and directional calling is not yet implemented.")
    }
    ## Error thrown for model with heterogeneity in source strength and a log-link function.
    if (fit.het.source){
        if (ss.link != "identity"){
            stop("Fitting of signal strength models with heterogeneity in source signal strength is only implemented with an identity link function.")
        }
    }
    if (first.calls){
        vectorise <- function(x) x[[1]]
    } else {
        vectorise <- list.to.vector
    }
    ## Stuff for the .dat file.
    data.list <- list(n_sessions = n.sessions,
                      survey_length = survey.length,
                      n_unique_per_sess = n.unique,
                      local = as.numeric(local),
                      n_local_per_unique = c(all.n.local, recursive = TRUE),
                      which_local_per_unique = c(all.which.local, recursive = TRUE),
                      n_D_betapars = n.D.betapars,
                      D_betapars_lb = D.betapars.lb,
                      D_betapars_ub = D.betapars.ub,
                      D_betapars_phase = D.betapars.phase,
                      D_betapars_sf = D.betapars.sf,
                      n_detpars = n.detpars,
                      detpars_lb = detpars.lb,
                      detpars_ub = detpars.ub,
                      detpars_phase = detpars.phase,
                      detpars_sf = detpars.sf,
                      detpars_linkfns = detpars.link,
                      n_suppars = n.suppars, suppars_lb = suppars.lb,
                      suppars_ub = suppars.ub,
                      suppars_phase = suppars.phase,
                      suppars_sf = suppars.sf,
                      suppars_linkfns = suppars.link,
                      detfn_id = detfn.id, trace = as.numeric(trace),
                      dbl_min = dbl.min, n_per_sess = n,
                      n_traps_per_sess = n.traps,
                      n_mask_per_sess = n.mask, A_per_sess = A,
                      capt_bin_unique = vectorise(capt.bin.unique),
                      capt_bin_freqs = c(capt.bin.freqs, recursive = TRUE),
                      fit_angs = as.numeric(fit.bearings),
                      fit_dir = as.numeric(fit.dir),
                      n_dir_quadpoints = n.dir.quadpoints,
                      fit_het_source = as.numeric(fit.het.source),
                      het_source_gh = as.numeric(het.source.gh),
                      n_het_source_quadpoints = n.het.source.quadpoints,
                      het_source_nodes = het.source.nodes,
                      het_source_weights = het.source.weights,
                      capt_ang = vectorise(capt.bearing),
                      fit_dists = as.numeric(fit.dists),
                      capt_dist = vectorise(capt.dist),
                      fit_ss = as.numeric(fit.ss), cutoff = cutoff,
                      first_calls = as.numeric(first.calls),
                      lower_cutoff = ifelse(is.null(lower.cutoff), 0, lower.cutoff),
                      linkfn_id = linkfn.id,
                      capt_ss = vectorise(capt.ss),
                      fit_toas = as.numeric(fit.toas),
                      capt_toa = vectorise(capt.toa),
                      fit_mrds = as.numeric(fit.mrds),
                      mrds_dist = vectorise(mrds.dist),
                      dists = vectorise(dists),
                      angs = vectorise(bearings),
                      toa_ssq = vectorise(toa.ssq),
                      mm_ihd = vectorise(mm.ihd))
    ## Determining whether or not standard errors should be calculated.
    if (!is.null(cue.rates)){
        fit.freqs <- any(cue.freqs != 1)
    } else {
        fit.freqs <- FALSE
    }
    ## Setting hess.
    if (is.null(hess)){
        hess <- !fit.freqs
    }
    ## Using optimx() for first call fits.
    if (first.calls){
        if (n.sessions > 1){
            stop("First-call models not yet implemented for multisession data.")
        }
        ## All possible capture histories.
        n.combins <- 2^n.traps
        combins <- matrix(NA, nrow = n.combins, ncol = n.traps)
        for (i in 1:n.traps){
            combins[, i] <- rep(rep(c(0, 1), each = 2^(n.traps - i)), times = 2^(i - 1))
        }
        data.list$combins <- combins
        fit <- optimx(c(sv.link[c("D.(Intercept)", "b0.ss", "b1.ss", "sigma.ss")], recursive = TRUE),
                      secr_nll, dat = data.list, get_esa = FALSE, method = "nmkb",
                      hessian = hess)
        out <- vector("list", 15)
        names(out) <- c("fn", "coefficients", "coeflist", "se", "loglik", "maxgrad", 
                        "cor", "vcov", "npar", "npar_re", "npar_sdrpt", "npar_rep",
                        "npar_total", "hes", "eratio")
        out$fn <- "optimx"
        c <- as.vector(coef(fit))
        n.opars <- length(c)
        coeflist <- list()
        for (i in 1:n.opars){
            coeflist[[i]] <- c[i]
        }
        names(coeflist) <- paste(colnames(coef(fit)), "_link", sep = "")
        out$coeflist <- coeflist
        ## Delta method for unlinked paramters.
        if (hess){
            vcov.link <- solve(attr(fit, "details")[1, "nhatend"][[1]])
            jacobian <- matrix(0, nrow = 2*n.opars, ncol = n.opars)
            jacobian[1:n.opars, ] <- diag(n.opars)
            jacobian[(n.opars + 1):(2*n.opars), ] <- diag(exp(c))
            vcov.all <- jacobian %*% vcov.link %*% t(jacobian)
            cor.all <- cov2cor(vcov.all)
            vcov.all <- rbind(vcov.all, NA)
            vcov.all <- cbind(vcov.all, NA)
            cor.all <- rbind(cor.all, NA)
            cor.all <- cbind(cor.all, NA)
            se.all <- c(sqrt(diag(vcov.all)), NA)
        } else {
            vcov.all <- matrix(NA, nrow = 2*n.opars + 1, ncol = 2*n.opars + 1)
            se.all <- rep(NA, 2*length(coeflist) + 2)
            cor.all <- matrix(NA, nrow = 2*n.opars + 1, ncol = 2*n.opars + 1)
        }
        rownames(vcov.all) <- colnames(vcov.all) <- rownames(cor.all) <-
            colnames(cor.all) <- names(se.all) <- c(paste("pars_link", 1:n.opars, sep = "."),
                                                    paste("par_ests", 1:n.opars, sep = "."),
                                                    paste("esa", 1:n.sessions, sep = "."))
        out$se <- se.all
        out$se[length(out$se)] <- out$se[n.opars + 1]
        out$loglik <- -fit$value
        out$maxgrad <- c(fit$kkt1, fit$kkt2)
        out$cor <- cor.all
        out$vcov <- vcov.all
        out$npar <- n.opars
        out$npar_re <- 0
        out$npar_sdrpt <- 0
        out$npar_rep <- 0
        out$npar_total <- out$npar + out$npar_re + out$npar_sdrpt + out$npar_rep
        out$eratio <- as.logical(NA)
        esa <- secr_nll(coef(fit), data.list, TRUE)
    } else {
        ## Idea of running executable as below taken from glmmADMB.
        ## Working out correct command to run from command line.
        if (exe.type == "new"){
            exe.name <- "secr_new"
        } else if (exe.type == "old"){
            exe.name <- "secr"
        } else if (exe.type == "test"){
            exe.name <- "secr_test"
        } else {
            stop("Argument 'exe.type' must be \"old\" or \"new\".")
        }
        prefix.name <- exe.name
        if (.Platform$OS == "windows"){
            os.type <- "windows"
            exe.name <- paste(prefix.name, ".exe", sep = "")
        } else if (.Platform$OS == "unix"){
            if (Sys.info()["sysname"] == "Linux"){
                os.type <- "linux"
            } else if (Sys.info()["sysname"] == "Darwin"){
                os.type <- "mac"
            } else {
                stop("Unknown OS type.")
            }
        } else {
            stop("Unknown OS type.")
        }
        ## Finding executable folder (possible permission problems?).
        exe.dir <- paste(system.file(package = "ascr"), "ADMB", "bin", os.type, sep = "/")
        exe.loc <- paste(exe.dir, exe.name, sep = "/")
        ## Creating command to run using system().
        curr.dir <- getwd()
        ## Creating temporary directory.
        temp.dir <- tempfile("ascr", curr.dir)
        dir.create(temp.dir)
        ## Cleaning up files on function exit.
        if (clean){
            on.exit({setwd(curr.dir)
                unlink(temp.dir, recursive = TRUE)})
        } else {
            on.exit({setwd(curr.dir)
                message("ADMB files found in:", "\n", temp.dir, "\n")})
        }
        setwd(temp.dir)
        ## Creating .pin and .dat files.
        write_pin("secr", sv.link)
        write_dat("secr", data.list)
        ## Creating link to executable.
        if (os.type == "windows"){
            file.copy(exe.loc, exe.name)
            dll1.name <- "libstdc++-6.dll"
            dll2.name <- "libgcc_s_dw2-1.dll"
            dll1.loc <- paste(exe.dir, dll1.name, sep = "/")
            dll2.loc <- paste(exe.dir, dll2.name, sep = "/")
            file.copy(dll1.loc, dll1.name)
            file.copy(dll2.loc, dll2.name)
        } else {
            file.symlink(exe.loc, exe.name)
        }
        ## Sorting out -cbs and -gbs.
        if (!is.null(cbs)){
            cbs.cmd <- paste(" -cbs", format(cbs, scientific = FALSE))
        } else {
            cbs.cmd <- NULL
        }
        if (!is.null(gbs)){
            gbs.cmd <- paste(" -gbs", format(gbs, scientific = FALSE))
        } else {
            gbs.cmd <- NULL
        }
        ## Running ADMB executable.
        cmd <- paste("./"[os.type != "windows"], exe.name,
                     " -ind secr.dat -ainp secr.pin", " -neldmead"[neld.mead],
                     " -nohess"[!hess], cbs.cmd, gbs.cmd, sep = "")
        if (trace){
            message(cmd, "\n")
        }
        if (os.type == "windows"){
            system(cmd, ignore.stdout = !trace, show.output.on.console = trace)
        } else {
            system(cmd, ignore.stdout = !trace)
        }
        ## Reading in model results.
        if (exe.type == "test"){
            prefix.name <- strsplit(list.files(), "\\.")[[which(substr(list.files(),
                                                                       nchar(list.files()) - 3,
                                                                       nchar(list.files())) == ".par")]][1]
        }
        out <- suppressWarnings(try(read.ascr(prefix.name), silent = TRUE))
        if (fit.ihd){
            out$se <- out$se[names(out$se) != "D"]
        }
        ## Getting ESAs from .rep file for better accuracy. Also getting mask densities.
        rep.pars <- read_rep("secr")$est
        esa <- rep.pars[substr(names(rep.pars), 1, 3) == "esa"]
        names(esa) <- NULL
        D.mask <- vector(mode = "list", length = n.sessions)
        for (i in 1:n.sessions){
            D.mask[[i]] <- rep.pars[substr(names(rep.pars), 1, 8 + nchar(i)) == paste("D_mask[", i, "]", sep = "")]
        }
        out$D.mask <- D.mask
        out$mm.ihd <- mm.ihd
        setwd(curr.dir)
        if (class(out)[1] == "try-error"){
            stop("Parameters not found. There was either a problem with the model fit, or the executable did not run properly.")
        }
        ## Warning for non-convergence.
        if (out$maxgrad < -0.1){
            warning("Failed convergence -- maximum gradient component is too large.")
        }
        ## Moving back to original working directory.
        setwd(curr.dir)
        ## Removing fixed coefficients from list.
        if (!hess){
            out$coeflist[c(D.betapars.phase, detpars.phase, suppars.phase, D.betapars.phase) == -1] <- NULL
        }
        if (fit.ihd){
            out$coeflist <- out$coeflist[names(out$coeflist) != "D"]
        }
    }
    ## Creating coefficients vector.
    est.pars <- c(D.betapars.names, detpar.names, suppar.names)[c(D.betapars.phase, detpars.phase, suppars.phase[any.suppars]) > -1]
    n.est.pars <- length(est.pars)
    out$coefficients <- numeric(2*n.est.pars + n.sessions + !fit.ihd)
    names(out$coefficients) <- c(paste(est.pars, "_link", sep = ""), est.pars, paste("esa.", 1:n.sessions, sep = ""), "D"[!fit.ihd])
    if (hess){
        names(out$se) <- names(out$coefficients)
    }
    for (i in 1:n.est.pars){
        out$coefficients[i] <- out$coeflist[[i]]
    }
    if (!fit.ihd){
        out$coefficients["D"] <- exp(out$coefficients[D.betapars.phase[1]])
    }
    for (i in 1:n.est.pars){
        out$coefficients[n.est.pars + i] <-
            unlink.list[[links[[est.pars[i]]]]](out$coeflist[[i]])
    }
    ## Adding extra components to list.
    if (detfn == "log.ss") detfn <- "ss"
    ## Putting in updated argument names.
    args <- vector(mode = "list", length = length(arg.names))
    names(args) <- arg.names
    for (i in arg.names){
        if (!is.null(get(i))){
            args[[i]] <- get(i)
        }
    }
    out$args <- args
    out$n.sessions <- n.sessions
    out$fit.types <- fit.types
    out$infotypes <- names(fit.types)[fit.types]
    out$detpars <- detpar.names
    out$suppars <- suppar.names
    out$D.betapars <- D.betapars.names
    out$phases <- phases
    out$par.links <- par.links
    out$par.unlinks <- par.unlinks
    out$fit.ihd <- fit.ihd
    ## Logical value for random effects in the detection function.
    out$re.detfn <- FALSE
    if (detfn == "ss"){
        if (get.par(out, "b2.ss") != 0 | get.par(out, "sigma.b0.ss") != 0){
            out$re.detfn <- TRUE
        }
    }
    ## Putting in esa estimate.
    out$coefficients[2*n.est.pars + (1:n.sessions)] <- esa
    ## Putting in call frequency information and correct parameter names.
    if (fit.freqs){
        mu.rates <- mean(cue.rates)
        Dc <- get.par(out, "D")
        Da <- Dc/mu.rates
        names.vec <- c(names(out[["coefficients"]]), "Da", "Dc", "mu.rates")
        coefs.updated <- c(out[["coefficients"]], Da, Dc, mu.rates)
        names(coefs.updated) <- names.vec
        out[["coefficients"]] <- coefs.updated
        ## Removing ses, cor, vcov matrices.
        cor.updated <- matrix(NA, nrow = length(names.vec),
                              ncol = length(names.vec))
        dimnames(cor.updated) <- list(names.vec, names.vec)
        vcov.updated <- matrix(NA, nrow = length(names.vec),
                               ncol = length(names.vec))
        dimnames(vcov.updated) <- list(names.vec, names.vec)
        if (hess){
            ses.updated <- c(out[["se"]], rep(NA, 3))
            max.ind <- length(names.vec) - 3
            cor.updated[1:max.ind, 1:max.ind] <- out[["cor"]]
            vcov.updated[1:max.ind, 1:max.ind] <- out[["vcov"]]
        } else {
            ses.updated <- rep(NA, length(names.vec))
        }
        names(ses.updated) <- names.vec
        out[["se"]] <- ses.updated
        out[["cor"]] <- cor.updated
        out[["vcov"]] <- vcov.updated
        if (trace){
            if (!hess){
                message("NOTE: Standard errors not calculated; use boot.ascr().", "\n")
            } else {
                message("NOTE: Standard errors are probably not correct; use boot.ascr().", "\n")
            }
        }
    } else {
        if (hess){
            ## Putting correct parameter names into se, cor, vcov.
            replace <- substr(rownames(out$vcov), 1, 8) == "par_ests"
            rownames(out$vcov)[replace] <-
                colnames(out$vcov)[replace] <- rownames(out$cor)[replace] <-
                colnames(out$cor)[replace] <- est.pars
            replace <- 1:length(est.pars)
            rownames(out$vcov)[replace] <-
                colnames(out$vcov)[replace] <- rownames(out$cor)[replace] <-
                colnames(out$cor)[replace] <- paste(est.pars, "_link", sep = "")
        } else {
            ## Filling se, cor, vcov with NAs.
            names.vec <- names(out[["coefficients"]])
            ses.updated <- rep(NA, length(names.vec))
            names(ses.updated) <- names.vec
            cor.updated <- matrix(NA, nrow = length(names.vec),
                                  ncol = length(names.vec))
            dimnames(cor.updated) <- list(names.vec, names.vec)
            vcov.updated <- matrix(NA, nrow = length(names.vec),
                                   ncol = length(names.vec))
            dimnames(vcov.updated) <- list(names.vec, names.vec)
            out[["se"]] <- ses.updated
            out[["cor"]] <- cor.updated
            out[["vcov"]] <- vcov.updated
        }
    }
    out$fit.freqs <- fit.freqs
    out$first.calls <- first.calls
    if (out$fn == "secr"){
        if (out$maxgrad < -0.01){
            warning("Maximum gradient component is large.")
        }
    }
    if (out$fn == "optimx"){
        if (any(!out$maxgrad)){
            warning("Nelder-Mead algorithm has not converged.")
        }
    }
    class(out) <- c("ascr", "admb")
    out
    }

## Aliasing old admbsecr() function name.
#' @rdname fit.ascr
#' @export
admbsecr <- fit.ascr

#' Parallelising ascr fits using a cluster
#'
#' Fits SECR models on different cores within a cluster.
#'
#' @param ... Lists with components comprising arguments for a call to
#' \link{fit.ascr}. Component names must be the argument names.
#' @param arg.list Alternatively, a list with components comprising
#' the lists of arguments, as above.
#' @inheritParams boot.ascr
#'
#' @return A list, where components are objects returned by
#' \link{fit.ascr}. There is one component for each list of arguments
#' provide in \code{...}.
#'
#' @examples
#' \dontrun{
#' ## Running the examples in the fit.ascr() documentation in parallel.
#' simple.capt <- example$capt["bincapt"]
#' simple.hn.args <- list(capt = simple.capt, traps = example$traps,
#'                        mask = example$mask, fix = list(g0 = 1))
#' simple.hr.args <- list(capt = simple.capt, traps = example$traps,
#'                        mask = example$mask, detfn = "hr")
#' bearing.capt <- example$capt[c("bincapt", "bearing")]
#' bearing.hn.args <- list(capt = bearing.capt, traps = example$traps,
#'                         mask = example$mask, fix = list(g0 = 1))
#' ## This will only run if you have 4 cores available, you may need
#' ## to alter n.cores as appropriate.
#' fits <- par.fit.ascr(n.cores = 4, simple.hn.args, simple.hr.args,
#'                      bearing.hn.args)
#' }
#'
#' @export
par.fit.ascr <- function(n.cores, ..., arg.list = NULL){
    if (n.cores > detectCores()){
        stop("The argument n.cores is greater than the number of available cores.")
    }
    if (is.null(arg.list)){
        arg.list <- list(...)
    }
    n.fits <- length(arg.list)
    FUN <- function(i, arg.list){
        out <- try(do.call(fit.ascr, arg.list[[i]]), silent = TRUE)
    }
    cluster <- makeCluster(n.cores)
    clusterEvalQ(cluster, {
        library(ascr)
    })
    out <- parLapplyLB(cluster, 1:n.fits, FUN, arg.list = arg.list)
    stopCluster(cluster)
    out
}

## Aliasing old par.admbsecr() function name.
#' @rdname par.fit.ascr
#' @export
par.admbsecr <- par.fit.ascr

## Roxygen code for NAMESPACE and datasets.

## Package imports for roxygenise to pass to NAMESPACE.
#' @import graphics grDevices parallel plyr Rcpp R2admb stats testthat truncnorm
#' @importFrom CircStats dvm rvm
#' @importFrom optimx optimx
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom matrixStats colProds
#' @importFrom mgcv gam
#' @importFrom mvtnorm rmvnorm
#' @importFrom secr make.capthist make.mask read.mask read.traps sim.popn
#' @importFrom utils example setTxtProgressBar txtProgressBar
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
#' @name example
#' @format A list.
#' @usage example
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
#' @name multi.example
#' @format A list.
#' @usage multi.example
#' @docType data
#' @keywords datasets
NULL
