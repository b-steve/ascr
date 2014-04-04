#' Create mask object
#'
#' Creates a mask object to use with the function
#' \code{\link{admbsecr}}.
#'
#' @param buffer The minimum distance between trap locations and the
#' edge of the generated mask.
#' @param ... Arguments to be passed to \link{make.mask}.
#' @inheritParams admbsecr
#'
#' @return An object of class \code{mask}.
#'
#' @seealso \link{convert.mask} to convert a mask compatible with the
#' \link{secr} package.
#'
#' @examples
#' mask <- create.mask(traps = example.traps, buffer = 20)
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
#' Creates a capture history object to use with the function \code{\link{admbsecr}}.
#'
#' @param captures A data frame.
#' @param n.traps The total number of traps.
#' 
#' @export
create.capt <- function(captures, n.traps = NULL){
    ids <- captures[, 2]
    traps <- captures[, 4]
    if (is.null(n.traps)){
        n.traps <- max(traps)
    } else if (any(traps > n.traps)){
        stop("Trap ID in argument 'captures' exceeds 'n.traps'.")
    }
    all.types <- c("bearing", "dist", "ss", "toa", "mrds")
    info.types <- all.types[all.types %in% colnames(captures)]
    out <- vector(mode = "list", length = length(info.types) + 1)
    names(out) <- c("bincapt", info.types)
    n <- length(unique(ids))
    for (i in 1:length(out)){
        out[[i]] <- matrix(0, nrow = n, ncol = n.traps)
    }
    for (i in 1:n){
        id <- unique(ids)[i]
        trig <- traps[ids == id]
        if (length(trig) != length(unique(trig))){
            msg <- paste("Ignoring that individual", id, "was detected by some traps more than once.")
            warning(msg)
        }
        out[["bincapt"]][i, trig] <- 1
        for (j in info.types){
            for (k in trig){
                out[[j]][i, k] <- captures[ids == id & traps == k, j][1]
            }
        }
    }
    out
}

#' Convert traps object
#'
#' Converts an \code{admbsecr} traps matrix to a \code{secr} traps
#' object.
#'
#' The returned object is suitable for use as the \code{traps}
#' argument of the function \link{make.capthist}.
#'
#' @inheritParams admbsecr
#'
#' @return An object of class \code{traps} comprising a data frame of
#' x- and y-coordinates, the detector type ('single', 'multi',
#' 'proximity', 'count', 'polygon' etc.), and possibly other
#' attributes.
#'
#' @examples
#' traps <- convert.traps(traps = example.traps)
#'
#' @export
convert.traps <- function(traps){
    n.traps <- nrow(traps)
    colnames(traps) <- c("x", "y")
    traps.df <- data.frame(names = 1:n.traps, traps)
    read.traps(data = traps.df, detector = "proximity")
}

#' Convert mask object
#'
#' Converts an \code{admbsecr} mask matrix to a \code{secr} mask
#' object.
#'
#' The returned object is suitable for use as the \code{mask}
#' argument of the function \link{secr.fit}.
#'
#' @inheritParams admbsecr
#'
#' @return An object of class \code{mask}.
#'
#' @examples
#' mask <- convert.mask(mask = example.mask)
#'
#' @export
convert.mask <- function(mask){
    read.mask(data = as.data.frame(mask))
}

#' Convert capture history object
#'
#' Converts an \code{admbsecr} capture history list to a \code{secr}
#' capthist object.
#'
#' The returned object is suitable for use as the \code{capthist}
#' argument of the function \link{secr.fit}.
#'
#' @param capthist Logical, if \code{TRUE}, a \code{capthist} object
#' is returned. Otherwise a data frame is returned, which is suitable
#' for the \code{captures} argument to the \link{make.capthist}
#' function.
#' @inheritParams admbsecr
#'
#' @return An object of class \link{capthist}.
#'
#' @examples
#' capt <- convert.capt(capt = example.capt, traps = example.traps)
#'
#' @export
convert.capt <- function(capt, traps, capthist = TRUE){
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
    if (!is.null(capt$ss)){
        out <- data.frame(out, ss = t(capt$ss)[t(capt$bincapt) == 1])
    }
    for (i in names(capt)[!(names(capt) %in% c("bincapt", "ss"))]){
        out[, i] <- t(capt[[i]])[t(capt$bincapt) == 1]
    }
    if (capthist){
        traps <- convert.traps(traps)
        out <- make.capthist(out, traps, fmt = "trapID", noccasions = 1)
    }
    out
}
