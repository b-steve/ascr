#' Create mask object
#'
#' Creates a mask object to use with the function
#' \code{\link[admbsecr]{admbsecr}}.
#'
#' @param traps A matrix with two columns. Each row provides Cartesian
#' coordinates for the location of a trap.
#' @param buffer The minimum distance between trap locations and the
#' edge of the generated mask.
#' @param ... Arguments to be passed to \code{link[secr]{make.mask}}.
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

#' Convert traps object
#'
#' Converts an \code{admbsecr} traps matrix to a \code{secr} traps
#' object.
#'
#' @param traps A matrix with two columns. Each row provides Cartesian
#' coordinates for the location of a trap.
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
#' @param mask A matrix with two columns. Each row provides Cartesian
#' coordinates for the location of a mask point.
#' @export
convert.mask <- function(mask){
    read.mask(data = as.data.frame(mask))
}

#' Convert capture history object
#'
#' Converts an \code{admbsecr} capture history list to a \code{secr}
#' capthist object.
#'
#' @param capt A capture history list.
#' @param traps A matrix with two columns. Each row provides Cartesian
#' coordinates for the location of a trap.
#' @param capthist Logical, if \code{TRUE}, a \code{capthist} object
#' is returned. Otherwise a data frame is returned, which is suitable
#' for the \code{captures} argument to the
#' \code{\link[secr]{make.capthist}} function.
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
    if (capthist){
        traps <- convert.traps(traps)
        out <- make.capthist(out, traps, fmt = "trapID", noccasions = 1)
    }
    out
}
