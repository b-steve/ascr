#' Create mask object
#'
#' Creates a mask object to use with the function
#' \code{\link{fit.ascr}}.
#'
#' @param buffer The minimum distance between trap locations and the
#'     edge of the generated mask.
#' @param ... Arguments to be passed to \link{make.mask}.
#' @inheritParams fit.ascr
#'
#' @return An object of class \code{mask}. The object has two useful attributes: \code{area} which gives the area of each grid cell in hectares and \code{buffer} which gives the value of \code{buffer} supplied to this function.
#'
#' @seealso \link{convert.mask} to convert a mask compatible with the
#'     \link{secr} package.
#'
#' @examples
#' mask <- create.mask(traps = example$traps, buffer = 20)
#'
#' # calculate the area of the mask (in hectares)
#' mask_area <- attr(mask, "area") * nrow(mask)
#'
#' @export
create.mask <- function(traps, buffer, ...){
    traps <- convert.traps(traps)
    mask <- make.mask(traps, buffer = buffer, type = "trapbuffer", ...)
    A <- attr(mask, "area")
    mask <- as.matrix(mask)
    attr(mask, "area") <- A
    attr(mask, "buffer") <- buffer
    mask
}

#' Creating capture history object.
#'
#' Creates a capture history object to use with the function
#' \code{\link{fit.ascr}}.
#'
#' The \code{captures} argument to this function is intended to be of
#' a similar format to the \code{captures} argument to
#' \link{make.capthist} in the \link{secr} package. That is, users can
#' use the same \code{captures} data frame with \code{create.capt} and
#' \code{make.capthist}, which generate capture histories for use with
#' the \code{ascr} and \link{secr} packages respectively.
#'
#' As such, the first, second, and fourth columns should provide the
#' session, the identification number of the detected animal or call,
#' and the trap number of the trap which made the detection (where the
#' trap number is the row number of the corresponding trap in the
#' matrix of trap locations), respectively. Note that, for the
#' \link{secr} package, the third column of the \code{captures} data
#' frame provides the 'occassion' of the detection for
#' \link{make.capthist}.  However, the ascr package does not presently
#' have the capabilities to deal with multi-occassion data, so the
#' third column is ignored by \code{create.capt}.
#'
#' Additional optional columns can specify the additional information
#' collected over the course of the survey:
#' \itemize{
#'   \item A column named \code{bearing} containing estimated bearings
#'         from which the detector detected the individual.
#'   \item A column named \code{dist} containing the estimated
#'         distance between the individual detected and the detector.
#'   \item A column named \code{ss} containing the measured signal
#'         strengh of the detected acoustic signal (only possible when
#'         detectors are microphones).
#'   \item A column named \code{toa} containing the measured time of
#'         arrival (in seconds) since the start of the survey (or some
#'         other reference time) of the detected acoustic signal (only
#'         possible when the detectors are microphones).
#'   \item A column named \code{mrds} containing the \emph{known} (not
#'         estimated) distance between the individual detected and the
#'         detector.
#' }
#'
#' @param captures A data frame of capture records, see 'Details' for
#'     the correct format.
#' @param n.traps The total number of traps. If \code{NULL} then the
#'     number of traps is assumed to be the largest value in the
#'     \code{trap} column of the \code{captures} argument.
#' @param n.sessions The total number of sessions. If \code{NULL} then
#'     the number of sessions is assumed to be the largest value in
#'     the \code{session} column of the \code{captures} argument.
#'
#' @export
create.capt <- function(captures, n.traps = NULL, n.sessions = NULL){
    session.full <- captures[, 1]
    id.full <- captures[, 2]
    trap.full <- captures[, 4]
    if (is.null(n.sessions)){
        n.sessions <- max(session.full)
    } else if (any(session.full > n.sessions)){
        stop("Session ID in arguments 'captures' exceeds 'n.sessions'.")
    }
    if (is.null(n.traps)){
        n.traps <- max(trap.full)
    } else if (any(trap.full > n.traps)){
        stop("Trap ID in argument 'captures' exceeds 'n.traps'.")
    }
    all.types <- c("bearing", "dist", "ss", "toa", "mrds")
    info.types <- all.types[all.types %in% colnames(captures)]
    if (n.sessions > 1){
        out.list <- vector(mode = "list", length = n.sessions)
    }
    captures.full <- captures
    for (s in 1:n.sessions){
        captures <- captures.full[session.full == s, ]
        session <- captures[, 1]
        id <- captures[, 2]
        trap <- captures[, 4]
        out <- vector(mode = "list", length = length(info.types) + 1)
        names(out) <- c("bincapt", info.types)
        n <- length(unique(id))
        for (i in 1:length(out)){
            out[[i]] <- matrix(0, nrow = n, ncol = n.traps)
        }
        rnames <- character(n)
        for (i in 1:n){
            u.id <- unique(id)[i]
            trig <- trap[id == u.id]
            if (length(trig) != length(unique(trig))){
                msg <- paste("Ignoring that individual", u.id, "was detected by some traps more than once.")
                warning(msg)
            }
            out[["bincapt"]][i, trig] <- 1
            for (j in info.types){
                for (k in trig){
                    out[[j]][i, k] <- captures[id == u.id & trap == k, j][1]
                }
            }
            rnames[i] <- u.id
        }
        for (i in 1:length(out)){
            rownames(out[[i]]) <- rnames
        }
        if (n.sessions > 1){
            out.list[[s]] <- out
        }
    }
    if (n.sessions > 1){
        out <- out.list
    }
    out
}

#' Convert traps object
#'
#' Converts an \code{ascr} traps matrix to a \code{secr} traps
#' object.
#'
#' The returned object is suitable for use as the \code{traps}
#' argument of the function \link{make.capthist}.
#'
#' @param ss Logical, set to \code{TRUE} if a signal strength
#'     detection function is to be used.
#' @inheritParams fit.ascr
#'
#' @return An object of class \code{traps} comprising a data frame of
#'     x- and y-coordinates, the detector type ('single', 'multi',
#'     'proximity', 'count', 'polygon' etc.), and possibly other
#'     attributes.
#'
#' @examples
#' traps <- convert.traps(traps = example$traps)
#'
#' @export
convert.traps <- function(traps, ss = FALSE){
    n.traps <- nrow(traps)
    colnames(traps) <- c("x", "y")
    traps.df <- data.frame(names = 1:n.traps, traps)
    detector <- ifelse(ss, "signal", "proximity")
    read.traps(data = traps.df, detector = detector)
}

#' Convert mask object
#'
#' Converts an \code{ascr} mask matrix to a \code{secr} mask
#' object.
#'
#' The returned object is suitable for use as the \code{mask}
#' argument of the function \link{secr.fit}.
#'
#' @inheritParams fit.ascr
#'
#' @return An object of class \code{mask}.
#'
#' @examples mask <- convert.mask(mask = example$mask)
#'
#' @export
convert.mask <- function(mask){
    read.mask(data = as.data.frame(mask))
}

#' Capture history conversion.
#'
#' These functions convert a capture history object between the
#' structures required for the \code{ascr} and \code{secr}
#' packages.
#'
#' @param capt A \code{secr} capture history object for
#'     \code{convert.capt.to.ascr}, or an \code{ascr} capture
#'     history object for \code{convert.capt.to.secr}.
#' @param capthist Logical, if \code{TRUE}, a \code{capthist} object
#'     is returned. Otherwise a data frame is returned, which is
#'     suitable for the \code{captures} argument to the
#'     \link{make.capthist} function.
#' @param cutoff The signal strength threshold for detection, if
#'     required.
#' @inheritParams fit.ascr
#'
#' @return A capture history object appropriate for analysis using
#'     either the \code{ascr} or the \code{secr} package.
#' @examples capt <- convert.capt.to.secr(capt = example$capt, traps = example$traps, cutoff = example$cutoff)
#'
#' @name convert.capt
NULL

#' @rdname convert.capt
#' @export
convert.capt.to.ascr <- function(capt){
    bincapt <- capt[, 1, ]
    nr <- nrow(bincapt)
    nc <- ncol(bincapt)
    out <- list(bincapt = bincapt)
    if (!is.null(attr(capt, "signalframe"))){
        ss.capt <- matrix(attr(capt, "signalframe")[, 1],
                          nrow = nr, ncol = nc)
        out$ss <- ss.capt
    } else {
        
    }
    out
}

## Aliasing old convert.capt.to.admbsecr() function name.
#' @rdname convert.capt
#' @export
convert.capt.to.admbsecr <- convert.capt.to.ascr

#' @rdname convert.capt
#' @export
convert.capt.to.secr <- function(capt, traps, capthist = TRUE, cutoff = NULL){
    n <- nrow(capt$bincapt)
    n.dets <- sum(capt$bincapt)
    session <- rep(1, n.dets)
    n.indiv.dets <- apply(capt$bincapt, 1, sum)
    ID <- rep(1:n, times = n.indiv.dets)
    occasion <- rep(1, n.dets)
    trap <- c(apply(capt$bincapt, 1, function(x) which(x == 1)),
              recursive = TRUE)
    names(trap) <- NULL
    out <- data.frame(session = session, ID = ID, occasion = occasion,
                      trap = trap)
    fit.ss <- !is.null(capt$ss)
    if (fit.ss){
        out <- data.frame(out, ss = t(capt$ss)[t(capt$bincapt) == 1])
    }
    for (i in names(capt)[!(names(capt) %in% c("bincapt", "ss"))]){
        out[, i] <- t(capt[[i]])[t(capt$bincapt) == 1]
    }
    if (capthist){
        traps <- convert.traps(traps, ss = fit.ss)
        out <- make.capthist(out, traps, fmt = "trapID", noccasions = 1,
                             cutval = cutoff)
    }
    out
}

#' Create a capture history object from a PAMGuard output file
#'
#' Converts a PAMGuard output file to a capture history object
#' suitable for use with the \link{fit.ascr} function. This uses
#' \link{make.acoustic.captures} to allocate call identities to
#' detections.
#'
#' @param dets Detection output dataframe from PAMGuard.
#' @param mics A matrix containing the coordinates of microphone
#'     locations.
#' @param time.range A vector of length two, providing a range of
#'     times for which a subset should be taken to create the capture
#'     history.
#' @param sound.speed The speed of sound in metres per second.
#' @param new.allocation Logical, if \code{TRUE}, an improved
#'     call-allocation method is used. The old version is retained so
#'     that older analyses can be replicated.
#' @export
convert.pamguard <- function(dets, mics, time.range = NULL,
                             sound.speed = 330, new.allocation = TRUE){
    toa.info <- dets$startSeconds
    mic.id <- log2(dets$channelMap) + 1
    ss.info <- dets$amplitude
    n <- nrow(dets)
    clicks <- data.frame(session = rep(1, n), ID = 1:n,
                         occasion = rep(1, n), trap = mic.id,
                         ss = ss.info, toa = toa.info)
    if (!is.null(time.range)){
        keep <- toa.info > time.range[1] & toa.info < time.range[2]
        if (!any(keep)){
            stop("No calls were detected within the specified time.range.")
        }
        clicks <- clicks[keep, ]
    }
    ord <- order(clicks$toa)
    clicks <- clicks[ord, ]
    ## Old and new way to allocate IDs below.
    if (new.allocation){
        captures <- clicks
        captures[, 2] <- allocate.calls(mics, clicks, sound.speed)
    } else {
        captures <- make.acoustic.captures(mics, clicks, sound.speed)
    }
    create.capt(captures)
}



