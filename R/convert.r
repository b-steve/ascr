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
#' mask <- create.mask(traps = example.data$traps, buffer = 20)
#'
#' # calculate the area of the mask (in hectares)
#' mask_area <- attr(mask, "area") * nrow(mask)
#'
#' @export
create.mask <- function(traps, buffer, ...){
    ## Changing data frame to matrix to avoid is.list() issues.
    if (is.data.frame(traps)){
        traps <- as.matrix(traps)
    }
    if (is.list(traps)){
        ## Calling create.mask() on each session separately.
        n.sessions <- length(traps)
        mask <- vector(mode = "list", length = n.sessions)
        for (i in 1:n.sessions){
            mask[[i]] <- create.mask(traps[[i]], buffer = buffer, ...)
        }
    } else {
        ## Creating mask for a single session.
        traps <- convert.traps(traps)
        mask <- secr::make.mask(traps, buffer = buffer, type = "trapbuffer", ...)
        A <- attr(mask, "area")
        mask <- as.matrix(mask)
        attr(mask, "area") <- A
        attr(mask, "buffer") <- buffer
    }
    mask
}


#' Creating capture history object.
#'
#' @param captures 
#' @param traps 
#' @param n.traps 
#' @param n.sessions 
#' @param mrds.locs 
#' @param use.name 
#'
#' @return
#' @export
create.capt <- function(captures, traps = NULL, n.traps = NULL, n.sessions = NULL,
                        mrds.locs = NULL, use.name = FALSE){
    
    #captures must be provided as a matrix or data frame with at least 4 columns
    #traps could be list, matrix or data frame
    #traps, or (n.traps & n.sessions), one of them must be provided
    stopifnot(is.data.frame(captures) | is.matrix(captures), ncol(captures) >= 4,
              any(is.null(traps), is.list(traps), is.matrix(traps), is.data.frame(traps)),
              any(!is.null(traps), !is.null(n.traps) & !is.null(n.sessions)))
    
    if ((!missing(n.traps) | !missing(n.sessions)) & is.null(traps)){
        warning("Arguments 'n.traps' and 'n.sessions' are deprecated.
    Please provide the the 'traps' argument instead. Future versions of ascr will require
            the 'traps' argument to be provided to the 'create.capt()' function.")
    }
    
    if(is.matrix(captures)) captures = as.data.frame(captures, stringsAsFactors = FALSE)
    
    #rename the 1, 2, 4 column, and make it sure that "session" and "trap" are numeric
    colnames(captures)[c(1,2,4)] = c("session", "ID", "trap")
    captures$session = as.numeric(captures$session)
    stopifnot(all(captures$session > 0), all(captures$session %% 1 == 0))
    stopifnot(all(captures$trap > 0), all(captures$trap %% 1 == 0))
    captures = captures[order(captures$session, captures$ID, captures$trap),]
    
    #------------------------------------------------------------------------------------------------------
    #deal with n.sessions
    
    #generate n.sessions based on "captures" data
    tem.n.sessions.capt = max(captures$session)
    
    if(!is.null(traps)){
        #generate n.sessions based on "traps" data
        tem.n.sessions.trap = ifelse(is(traps, 'list'), length(traps), 1)
        
        #length of traps is 1 is a special case
        if(tem.n.sessions.trap == 1){
            #if n.sessions is not provided, we assume it should be just 1
            if(is.null(n.sessions)) n.sessions = 1
            #if n.sessions is not null, we just use it, no code needed
        } else {
            #if length of traps is not 1, we should use it as n.sessions
            if(!is.null(n.sessions)) warning("'traps' is provided, argument 'n.sessions' will be overwriten")
            n.sessions = tem.n.sessions.trap
        }
    }
    
    #if "traps" is null, the constraint as very beginning of this function has confirmed that
    #n.sessions will be provided, so just use it, no code needed
    
    #the maximum of session in captures data should be <= n.sessions
    if(tem.n.sessions.capt > n.sessions) {
        stop('The number of sessions based on "n.sessions" or "traps" can not match the captures data')
    }
    
    
    if(!all(captures$session %in% seq(n.sessions))){
        stop('"session" must be assigned as successive positive integers begins from one.')
    }
    
    
    #-----------------------------------------------------------------------------------------------
    
    #deal with n.traps
    
    #generate n.traps based on captures data frame
    tem.n.traps.capt = numeric(n.sessions)
    tem.n.traps.capt[unique(captures$session)] = aggregate(captures$trap,
                                                           list(session = captures$session), max)$x
    
    #check if "trap" is labeled by natural numbers
    for(i in 1:n.sessions){
        tem = subset(captures, session == i)
        if(nrow(tem) > 0){
            if(any(!tem$trap %in% 1:tem.n.traps.capt[i])) {
                stop('labels of traps must be successive natural numbers')
            }
        }
    }
    
    
    if(!is.null(traps)){
        #if 'traps' is provided deal with the special case that its length is 1
        if(any(is(traps, 'list') & length(traps) == 1, is.matrix(traps), is.data.frame(traps))){
            tem.traps = vector('list', n.sessions)
            if(is(traps, 'list')){
                for(i in 1:n.sessions) tem.traps[[i]] = traps[[1]]
            } else {
                for(i in 1:n.sessions) tem.traps[[i]] = traps
            }
            traps = tem.traps
        }
        
        if(!is.null(n.traps)) warning("'traps' is provided, argument 'n.traps' will be overwriten")
        n.traps = sapply(traps, nrow)
    } else {
        #length(n.traps) == 1 is a special case
        if(length(n.traps) == 1){
            n.traps = rep(n.traps, n.sessions)
        } else {
            msg = paste0("The length of 'n.traps' should match the number of sessions: ", n.sessions)
            if(length(n.traps) != n.sessions) stop(msg)
        }
    }
    
    #n.traps generated from capture data should always smaller than it from "traps"
    if(any(tem.n.traps.capt > n.traps)){
        msg = paste0('In session(s): (', paste(which(tem.n.traps.capt > n.traps), collapse = ","),
                     '), the number of traps based on "n.traps" or "traps" cannot match the captures data')
        stop(msg)
    }
    
    stopifnot(all(n.traps) > 0)
    
    #-----------------------------------------------------------------------------------------------------
    
    #deal with mrds.locs. "ID" will be dealt with together as they are closely related
    #here we only require mrds.locs for the sessions that we have detection,
    #so not using n.sessions for checking
    is.mrds <- !is.null(mrds.locs)
    
    if (is.mrds){
        if(is(mrds.locs, "list")) {
            if(n.sessions != length(mrds.locs)) {
                msg = paste0("mrds.locs must be a list with a length of n.sessions: ", n.sessions)
                stop(msg)
            }
        } else {
            stopifnot(n.sessions == 1, is.matrix(mrds.locs) | is.data.frame(mrds.locs))
            mrds.locs = list(mrds.locs)
        }
        
    }
    
    
    #we added two new components, "animal_ID" column and "use named ID"
    is.animalID = "animal_ID" %in% colnames(captures)
    
    #deal with new components with all possible combinations
    if(is.animalID | use.name == TRUE){
        
        if(is.animalID){
            captures$key = paste(captures$animal_ID, captures$ID, sep = "_")
        } else {
            captures$key =as.character(captures$ID)
        }
        
        dupli = duplicated(captures[,c("session", "key", "trap")])
        if(any(dupli)){
            warning("Ignoring the duplicated observations")
        }
        captures = captures[!dupli,]
        
        n.keys = numeric(n.sessions)
        n.keys[unique(captures$session)] = aggregate(captures$key, list(session = captures$session),
                                                     function(x) length(unique(x)))$x
        
        
        if(is.mrds) {
            tem.n.IDs = lapply(mrds.locs, nrow)
            tem.n.IDs = lapply(tem.n.IDs, function(x) ifelse(is.null(x), 0, x))
            tem.n.IDs = do.call('c', tem.n.IDs)
            stopifnot(all(tem.n.IDs == n.keys))
            
            for(i in unique(captures$session)){
                if(is.animalID){
                    if(!all(c("animal_ID", "ID") %in% colnames(mrds.locs[[i]]))) {
                        stop("Please provide identification information 'animal_ID' and 'ID' in mrds.locs argument")
                    }
                } else {
                    if(!("ID" %in% colnames(mrds.locs[[i]]))){
                        stop("Please provide identifiation information 'ID' in mrds.locs argument")
                    }
                }
            }
            
        }
        
    } else {
        dupli = duplicated(captures[,c("session", "ID", "trap")])
        if(any(dupli)){
            warning("Ignoring the duplicated observations")
        }
        captures = captures[!dupli,]
        n.IDs = numeric(n.sessions)
        n.IDs[unique(captures$session)] = aggregate(captures$ID, list(session = captures$session),
                                                    function(x) length(unique(x)))$x
        
        if(is.mrds) {
            tem.n.IDs = lapply(mrds.locs, nrow)
            tem.n.IDs = lapply(tem.n.IDs, function(x) ifelse(is.null(x), 0, x))
            tem.n.IDs = do.call('c', tem.n.IDs)
            stopifnot(all(tem.n.IDs == n.IDs))
            for(i in 1:n.sessions){
                tem = subset(captures, session == i)
                if(any(!tem$ID %in% 1:n.IDs[i])) {
                    stop('labels of IDs with "mrds.locs" provided must be successive natural numbers')
                }
            }
            
        }
    }
    
    
    #-----------------------------------------------------------------------------------------------------
    #all checks have been done, below we generate output list
    
    all.types <- c("bearing", "dist", "ss", "toa")
    info.types <- all.types[all.types %in% colnames(captures)]
    out.list <- vector(mode = "list", length = n.sessions)
    
    
    
    
    for (s in 1:n.sessions){
        tem <- captures[captures$session == s, ]
        if(is.animalID){
            tem.key = data.frame(animal_ID = tem$animal_ID, ID = tem$ID)
            tem.key = tem.key[!duplicated(tem.key)]
            id <- paste(tem.key$animal_ID, tem.key$ID, sep = '_')
            n <- n.keys[s]
        } else {
            id <- tem$ID
            n <- ifelse(use.name == TRUE, n.keys[s], n.IDs[s])
        }
        
        out <- vector(mode = "list", length = length(info.types) + 1 + as.numeric(is.mrds))
        #if animal_ID is used, we use data frame for 'bincapt', and we add two columns after
        #n.traps[s] columns, so that the order of the first n.traps[s] columns still make sense
        
        
        for (i in 1:length(out)){
            out[[i]] <- matrix(0, nrow = n, ncol = n.traps[s])
        }
        
        
        if(is.animalID){
            out[[1]] <- as.data.frame(out[[1]], stringsAsFactors = FALSE)
            out[[1]]$animal_ID = NA
            out[[1]]$ID = NA
        }
        names(out) <- c("bincapt", info.types, "mrds"[is.mrds])
        
        
        if (nrow(tem) > 0){
            trap <- tem$trap
            rnames <- character(n)
            for (i in 1:n){
                u.id <- unique(id)[i]
                trig <- trap[id == u.id]
                out[["bincapt"]][i, trig] <- 1
                for (j in info.types){
                    for (k in trig){
                        out[[j]][i, k] <- tem[id == u.id & trap == k, j][1]
                    }
                }
                if (is.mrds){
                    out[["mrds"]] <- mrds.locs[[s]]
                }
                rnames[i] <- u.id
                
                if(is.animalID){
                    
                    out[["bincapt"]]$animal_ID = tem.key$animal_ID
                    out[["bincapt"]]$ID = tem.key$ID
                    
                }
                
            }
            for (i in 1:length(out)){
                rownames(out[[i]]) <- rnames
            }
        }
        out.list[[s]] = out
    }
    
    if(n.sessions == 1) {
        return(out.list[[1]])
    } else {
        return(out.list)
    }
    
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
#' traps <- convert.traps(traps = example.data$traps)
#'
#' @export
convert.traps <- function(traps, ss = FALSE){
    if (is.list(traps)){
        stop("The convert.traps() function will only convert single-session trap objects.")
    }
    n.traps <- nrow(traps)
    colnames(traps) <- c("x", "y")
    traps.df <- data.frame(names = 1:n.traps, traps)
    detector <- ifelse(ss, "signal", "proximity")
    secr::read.traps(data = traps.df, detector = detector)
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
#' @examples mask <- convert.mask(mask = example.data$mask)
#'
#' @export
convert.mask <- function(mask){
    if (is.list(mask)){
        stop("The convert.mask() function will only convert single-session mask objects.")
    }
    secr::read.mask(data = as.data.frame(mask))
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
#' @examples capt <- convert.capt.to.secr(capt = example.data$capt, traps = example.data$traps, cutoff = example.data$cutoff)
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
    if (!any(names(capt) == "bincapt")){
        if (length(capt) > 1){
            stop("The convert.capt.to.secr() function will only convert single-session capture history objects.")
        }
        capt <- capt[[1]]
    }
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
        out <- secr::make.capthist(out, traps, fmt = "trapID", noccasions = 1,
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
    mics <- as.matrix(mics)
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
    create.capt(captures, traps = mics)
}


#' Create a capture history object from a Raven output file
#'
#' Converts a Raven output file to a capture history object
#' suitable for use with the \link{fit.ascr} function. This uses
#' \link{make.acoustic.captures} to allocate call identities to
#' detections.
#'
#' @param dets Detection output dataframe from Raven.
#' @inheritParams convert.pamguard
#' @export
convert.raven <- function(dets, mics, time.range = NULL, sound.speed = 330,
                          new.allocation = TRUE){
    mics <- as.matrix(mics)
    toa.info <- dets[, 4]
    mic.id <- log2(dets[, 3]) + 1
    ss.info <- dets[, 8]
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
    create.capt(captures, traps = mics)
}


#' Assigning ID numbers to sounds
#'
#' Identifies recaptures and assigns ID numbers to sounds recorded for
#' an SECR model.
#'
#' Detected sounds are assumed to come from the same animal if times
#' of arrival at different microphones are closer together than the
#' time it would take for sound to travel between these microphones.
#'
#' @param mics a matrix containing the coordinates of trap locations.
#' @param dets a data frame containing (at least): (i) \code{$toa},
#'     the precise time of arrival of the received sound, and (ii)
#'     \code{$trap} the trap at which the sound was recorded.
#' @param sound.speed the speed of sound in metres per second.
#' @return A data frame. Specifically, the \code{dets} dataframe, now
#'     with a new variable, \code{ID}.
#' @author David Borchers
#'
#' @export
make.acoustic.captures <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    dists <- distances(mics, mics)
    dt <- dists/sound.speed
    K <- dim(mics)[1]
    captures <- dets
    ct <- rep(-Inf, K)
    ID <- 1
    ct[dets$trap[1]] <- dets$toa[1]
    new <- FALSE
    ndets <- length(dets$toa)
    for (i in 2:ndets){
        if (ct[dets$trap[i]] > -Inf){
            nd <- length(which(ct > -Inf))
            captures$ID[(i - nd):(i - 1)] <- ID
            ct <- rep(-Inf, K)
            ct[dets$trap[i]] <- dets$toa[i]
            ID <- ID + 1
            if(i == ndets) captures$ID[i] <- ID
        }
        else {
            ct[dets$trap[i]] <- dets$toa[i]
            ctset <- which(ct > -Inf)
            dts <- dt[ctset, dets$trap[i]]
            cts <- -(ct[ctset] - dets$toa[i])
            if (any((cts - dts) > 0)) new <- TRUE
            if (new) {
                nd <- length(which(ct > -Inf)) - 1
                captures$ID[(i - nd):(i - 1)] <- ID
                ct <- rep(-Inf, K)
                ct[dets$trap[i]] <- dets$toa[i]
                ID <- ID + 1
                new <- FALSE
                if (i == ndets) captures$ID[i] <- ID
            } else if(i == ndets){
                nd <- length(which(ct > -Inf))
                captures$ID[(i - nd + 1):i] <- ID
            }
        }
    }
    captures
}



allocate.calls <- function(mics, dets, sound.speed){
    mics <- as.matrix(mics)
    trap.dists <- distances(mics, mics)
    n.dets <- nrow(dets)
    ## Allocating pairwise plausibility of common cue sources.
    dist.mat <- detection_dists(trap.dists, dets$trap)
    timediff.mat <- detection_timediffs(dets$toa, dets$trap)
    maxtime.mat <- dist.mat/sound.speed
    match.mat <- timediff.mat <= maxtime.mat
    ## Finding blocks of multiple cues with possible common sources.
    incomplete.blocks <- find_incomplete_blocks(match.mat)
    n.blocks <- max(incomplete.blocks)
    complete.block <- logical(n.blocks)
    final.mat <- matrix(FALSE, nrow = n.dets, ncol = n.dets)
    ## Allocating possible common cues to sources.
    reqss.mat <- dist.mat/timediff.mat
    for (i in 1:max(incomplete.blocks)){
        ## Grabbing a block.
        block <- match.mat[incomplete.blocks == i, incomplete.blocks == i]
        reqss <- reqss.mat[incomplete.blocks == i, incomplete.blocks == i]
        ## Working out if there is any possible ambiguity.
        is.complete <- all(block)
        ## If ambiguity, resolve it.
        if (!is.complete){
            block <- blockify(block, reqss)
        }
        final.mat[incomplete.blocks == i, incomplete.blocks == i] <- block
    }
    find_incomplete_blocks(final.mat)
}