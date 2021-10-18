create_capt <- function(captures, traps = NULL, n.traps = NULL, n.sessions = NULL,
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
