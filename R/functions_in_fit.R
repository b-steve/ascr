capture.fun = function(capt, animal.model){
  all.types <- c("bearing", "dist", "ss", "toa")
  
  if(!animal.model){
            if("bincapt" %in% names(capt)){
              n.sessions = 1
              n.traps = ncol(capt$bincapt)
              n.IDs = nrow(capt$bincapt)
              bucket_info <- all.types[all.types %in% names(capt)]
              is.mrds = "mrds" %in% names(capt)
              
            } else {
              n.sessions = length(capt)
              n.traps = sapply(capt, function(x) ncol(x$bincapt))
              
              n.IDs = sapply(capt, function(x) nrow(x$bincapt))
              bucket_info <- all.types[all.types %in% names(capt[[1]])]
              is.mrds = "mrds" %in% names(capt[[1]])
              
            }
            
            tem.data.capt = vector('list', n.sessions)
            for(i in 1:n.sessions){
              if(n.sessions == 1){
                tem = capt
              } else {
                tem = capt[[i]]
              }
              
              number.row = n.IDs[i] * n.traps[i]
              tem.df = data.frame(session = rep(i, number.row), ID = numeric(number.row))
              for(j in c("trap", "bincapt", all.types, "mrds_x", "mrds_y")) tem.df[[j]] = numeric(number.row)
              
              if(number.row > 0){
                tem.df$ID = rep(1:nrow(tem$bincapt), n.traps[i])
                tem.df$trap = rep(1:n.traps[i], each = n.IDs[i])
                for(k in c('bincapt', bucket_info)) tem.df[[k]] = as.vector(tem[[k]])
                
                if(is.mrds){
                  tem$mrds = as.data.frame(tem$mrds, stringsAsFactors = FALSE)
                  colnames(tem$mrds) = c('mrds_x', 'mrds_y')
                  tem$mrds$ID = 1:nrow(tem$mrds)
                  tem.df = merge(tem.df, tem$mrds, by = "ID")
                }
                
              }
              tem.data.capt[[i]] = tem.df
            }
            data.capt = do.call("rbind", tem.data.capt)
            #n.animals is useless in this case, just use 0 to occupy the variable name
            n.animals = numeric(n.sessions)
  #####################################################################################################          
  } else {
  #####################################################################################################
            #according to the data cleaning process in "create.capt()", all these indices will be
            #successive natural numbers, so take its max is equivalent to take length(unique(xxx))
            
            base_cols = c('session', 'animal_ID', 'ID', 'trap', 'bincapt')
            all_cols = c(base_cols, 'mrds_x', 'mrds_y', all.types)
            
            if(any(!base_cols %in% colnames(capt))){
              stop(paste0('columns: ', paste(base_cols, collapse = ", "), ", must be included."))
            }
    
            #for animal included model, since there is no compatible concern, the "create.capt()" has
            #already sorted out the data structure, here we just need to extract some basic information
            
            n.sessions = max(capt$session)
            n.traps = agg_sort(capt, 'trap', 'session', max)$x
            
            n.animals = agg_sort(capt, 'animal_ID', 'session', function(x) length(unique(x)))$x

            no.det = agg_sort(capt, 'animal_ID', 'session', function(x) any(is.na(x)))$x

            
            n.animals = ifelse(no.det, 0, n.animals)
            
            #n.IDs becomes the number of length(unique(animal_ID, ID)) for each session
            n.IDs = agg_sort(capt, c('animal_ID', 'ID'), 'session', function(x) length(unique(x)))$x
            n.IDs = ifelse(no.det, 0, n.IDs)
            
            is.mrds = "mrds_x" %in% colnames(capt)
            bucket_info = all.types[all.types %in% colnames(capt)]
            
            data.capt = capt[, colnames(capt) %in% all_cols]
            for(i in c('mrds_x', 'mrds_y', all.types)){
              if(!i %in% colnames(data.capt)) data.capt[[i]] = 0
            }
  }
  
  if(animal.model){
    data.capt$animal_ID = as.character(data.capt$animal_ID)
  }
  
  data.capt$ID = as.character(data.capt$ID)
  
  return(list(data.capt = data.capt, dims = list(n.sessions = n.sessions, n.traps = n.traps,
                                                 n.IDs = n.IDs, n.animals = n.animals),
              bucket_info = c(bucket_info, "mrds"[is.mrds])))
  
}


#################################################################################################

trap.fun = function(traps, dims){
  stopifnot(any(is(traps, 'list'), is.data.frame(traps), is.matrix(traps)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  #make traps to a list with length == n.sessions
  traps = df_to_list(traps, n.sessions)
  stopifnot(all(sapply(traps, nrow) == n.traps), all(sapply(traps, ncol) == 2))
  
  
  for(i in 1:n.sessions) {
    traps[[i]] = as.data.frame(traps[[i]])
    colnames(traps[[i]]) = c("trap_x", "trap_y")
    traps[[i]]$session = i
    traps[[i]]$trap = 1:nrow(traps[[i]])
  }
  return(do.call("rbind", traps))
  
}


###################################################################################

mask.fun = function(mask, dims, animal.model, data.traps, data.full, bucket_info, local, sound.speed){
  
  stopifnot(any(is(mask, 'list'), is.data.frame(mask), is.matrix(mask)))
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  n.IDs = dims$n.IDs
  
  mask = df_to_list(mask, n.sessions)
  
  n.masks <- sapply(mask, nrow)
  
  A <- numeric(n.sessions)
  buffer <- numeric(n.sessions)
  #dists records the distances between every mask to every trap
  #thetas here records the angle of the line between every mask and trap
  #in the original package, it is called "bearings" but there is already
  #a function occupies this name, so I changed it
  dists <- vector(mode = "list", length = n.sessions)
  thetas = vector("list", length = n.sessions)
  data.dists.thetas = vector('list', n.sessions)
  data.mask = vector('list', n.sessions)
  
  for (i in 1:n.sessions){
    traps = data.traps[data.traps$session == i, c("trap_x", "trap_y")]
    
    if(ncol(mask[[i]]) != 2) {
      stop("each element of 'mask' must be a matrix or data frame with 2 columns")
    }
    
    colnames(mask[[i]]) = c('x', 'y')
    
    
    data.mask[[i]] = as.data.frame(mask[[i]], stringsAsFactors = FALSE)
    data.mask[[i]]$session = i
    data.mask[[i]]$mask = seq(n.masks[i])
    
    
    mask[[i]] <- as.matrix(mask[[i]])
    dists[[i]] <- distances(as.matrix(traps), mask[[i]])
    #based on the main procedure, at this stage, we don't know whether it is "ss.dir",
    #since it should be determined at later stage when deal with "ss"
    #so we generate bearings anyway.
    thetas[[i]] = bearings(as.matrix(traps), mask[[i]])
    
    data.dists.thetas[[i]] = data.frame(session = rep(i, n.masks[i] * n.traps[i]),
                                        mask = rep(1:n.masks[i], each = n.traps[i]),
                                        trap = rep(1:n.traps[i], n.masks[i]),
                                        dx = as.vector(dists[[i]]),
                                        theta = as.vector(thetas[[i]]))
    
    
    ## Filling in area and buffer, if missing.
    if (is.null(attr(mask[[i]], "area"))){
      mask.dists <- distances(mask[[i]], mask[[i]])
      A[i] <- min(mask.dists[mask.dists > 0])^2/10000
    } else {
      A[i] <- attr(mask[[i]], "area")
    }
    
    if (is.null(attr(mask[[i]], "buffer"))){
      buffer[i] <- max(apply(dists[[i]], 2, min))
    } else {
      buffer[i] <- attr(mask[[i]], "buffer")
    }
    attr(mask[[i]], "area") = A[i]
    attr(mask[[i]], "buffer") = buffer[i]
  }
  #data.mask and data.dists.thetas are naturally in right order here
  data.dists.thetas = do.call('rbind', data.dists.thetas)
  data.mask = do.call('rbind', data.mask)
  
  
  ################################################################################
  #deal with local argument
  
  all.which.local <- vector(mode = "list", length = n.sessions)
  all.n.local <- vector(mode = "list", length = n.sessions)
  
  if ("mrds" %in% bucket_info){
    local <- TRUE
    for (i in 1:n.sessions){
      
      if(n.IDs[i] > 0){
        mrds = data.full[data.full$session == i, c("mrds_x", "mrds_y")]
        mrds = mrds[!duplicated(mrds),]
        all.which.local[[i]] <- find.nearest.mask(as.matrix(mrds), mask[[i]])
        all.n.local[[i]] <- plyr::laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      }
      
    }
  } else if (local){
    
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        if(animal.model){
          #for each animal_ID, find the detectors have at least one record of any call of this animal

          tem = data.full[data.full$session == i, c("bincapt", "animal_ID", "trap")]
          tem = agg_sort(tem, 'bincapt', c('animal_ID', 'trap'), function(x) as.numeric(any(x == 1)))
          bincapt = tem[, c('animal_ID', 'x')]
          colnames(bincapt) = c('animal_ID', 'bincapt')
          bincapt = agg_sort(bincapt, 'bincapt', 'animal_ID', function(x) t(x))
          bincapt = bincapt[, -1]
        } else {
          bincapt = data.full[data.full$session == i, c("bincapt","ID")]
          bincapt = agg_sort(bincapt, 'bincapt', "ID", function(x) t(x))
          bincapt = bincapt[, -1]
        }
        
        all.which.local[[i]] <- find_local(bincapt, dists[[i]], buffer[i])
        all.n.local[[i]] <- plyr::laply(all.which.local[[i]], length)
        all.which.local[[i]] <- c(all.which.local[[i]], recursive = TRUE)
      }
    }
  } else {
    #if not 'local' these two variables are meaningless, just assign them something
    for (i in 1:n.sessions){
      if(n.IDs[i] > 0){
        all.n.local[[i]] <- rep(1, n.IDs[i])
        all.which.local[[i]] <- rep(0, n.IDs[i])
      }
    }
  }

  #prepare a data set for ID*mask
  
  data.ID_mask = vector('list', n.sessions)
  if(local) {
    bucket_info = c(bucket_info, "local")
    tem.df = vector('list', n.sessions)
  }
  
  #to avoid confusing, u.id is animal_ID-ID if animal.model, and the "ID" column is
  #temporarily to be animal_ID-ID as well until in the later step where we restore it
  for(i in 1:n.sessions){
    if(n.IDs[i] > 0){
      tem = subset(data.full, session == i)

      if(animal.model){
        u.id = unique(paste(tem[['animal_ID']], tem[['ID']], sep = '---'))
        u.a = sort(unique(tem$animal_ID))
      } else {
        u.id = unique(tem[["ID"]])
      }
      
      if(local){
        if(animal.model){
          tem.df[[i]] = data.frame(session = rep(i, sum(all.n.local[[i]])),
                                   animal_ID = rep(u.a, all.n.local[[i]]),
                                   mask = all.which.local[[i]], local = 1,
                                   stringsAsFactors = F)
        } else {
          tem.df[[i]] = data.frame(session = rep(i, sum(all.n.local[[i]])),
                                   ID = rep(u.id, all.n.local[[i]]),
                                   mask = all.which.local[[i]], local = 1,
                                   stringsAsFactors = F)
        }
        
        #here column 'ID' is "animal_ID---ID" if animal.model, will be separated later
        data.ID_mask[[i]] = data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id, each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       stringsAsFactors = F)
      } else {
        #here column 'ID' is "animal_ID---ID" if animal.model, will be separated later
        data.ID_mask[[i]] = data.frame(session = rep(i, n.IDs[i] * n.masks[i]),
                                       ID = rep(u.id, each = n.masks[i]),
                                       mask = rep(1:n.masks[i], n.IDs[i]),
                                       local = 1, stringsAsFactors = F)
      }
      
      #assign 'toa_ssq' in this ID_mask data set
      #data.full must be sorted by 'session', 'ID', 'trap'
      #far before this step, it has already been sorted
      #but in case of any disorder, we sort the aggregated result again
      #we take 'order' so careful because 'aggregate' returns
      #results based on sorted 'list(ID = xx)' argument
      if('toa' %in% bucket_info){
        toa = agg_sort(tem, 'toa', c('animal_ID'[animal.model], 'ID'), function(x) t(x))
        toa.ssq = make_toa_ssq(toa$x, dists[[i]], sound.speed)
        #based on the designed structure of 'data.ID_mask' and the structure
        #of toa.ssq, it need to be transposed
        data.ID_mask[[i]][['toa_ssq']] = as.vector(t(toa.ssq))
      } else {
        data.ID_mask[[i]][['toa_ssq']] = 0
      }
    }
  }
  
  data.ID_mask = do.call('rbind', data.ID_mask)
  
  if(animal.model){
    data.ID_mask[['animal_ID']] = gsub("---(.+)", "", data.ID_mask[['ID']])
    data.ID_mask[['ID']] = gsub("^(.+)---", "", data.ID_mask[["ID"]])
  }
  
  
  
  if(local){
    tem.df = do.call('rbind', tem.df)
    if(animal.model){
      data.ID_mask = merge(data.ID_mask, tem.df, by = c('session', 'animal_ID', 'mask'), all = T)
    } else {
      data.ID_mask = merge(data.ID_mask, tem.df, by = c('session', 'ID', 'mask'), all = T)
    }
    
    data.ID_mask[['local']] = ifelse(is.na(data.ID_mask[['local']]), 0, data.ID_mask[['local']])
  }
  
  
  return(list(mask = mask, dists = dists, thetas = thetas, n.masks = n.masks, A = A,
              buffer = buffer, all.which.local = all.which.local, all.n.local= all.n.local,
              data.full = data.full, data.dists.thetas = data.dists.thetas,
              data.mask = data.mask, data.ID_mask = data.ID_mask, bucket_info = bucket_info))
}



###############################################################################

par.extend.fun = function(par.extend, data.full, data.mask, animal.model, dims, fulllist.par, extra_args){
  
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  n.masks = dims$n.masks
  #it is possible that throughout all sessions, the traps/masks remain the same
  #and in this case, users can only provide one set of data frame of traps/masks

  identical_traps = all(n.traps == n.traps[1])
  identical_masks = all(n.masks == n.masks[1])
  
  
  ihd_merge = FALSE
  #if the user is using ADMB version input for ihd model, there will be an argument named 'ihd.opts'
  #in the extra_args, merge it into 'par.extend'
  if('ihd.opts' %in% names(extra_args)){
    #some basic check about 'ihd.opts', this part cannot be done in the later check of 'par.extend'
    #because it is possible that the user provide 'par.extend' as well, and inside of it, there is
    #some 'mask' related extension
    ihd_data = extra_args$ihd.opts$covariates
    ihd_scale = extra_args$ihd.opts$scale
    ihd_model = extra_args$ihd.opts$model
    stopifnot(!is.null(ihd_data) & !is.null(ihd_model))
    if(is.null(ihd_scale)) ihd_scale = TRUE
    
    
    
    if(is.null(par.extend)){
      par.extend = list(data = list(mask = ihd_data), model = list(D=ihd_model), scale = ihd_scale)
    } else {
      if(!is.null(par.extend$model$D)){
        stop('"ihd.opts" and Density extension in "par.extend" are both provided, please choose one of them.')
      } else {
        par.extend$model$D = ihd_model
        par.extend$scale = ihd_scale
        if(is.null(par.extend$data$mask)){
          par.extend$data$mask = ihd_data
        } else {
          ihd_data = covariates_mask_check(ihd_data, n.sessions, n.masks, identical_masks)
          ihd_merge = TRUE
        }
      }
    }
  }
  
  #'fgam' is an object for compatibility with ADMB model
  fgam = NULL
  
  #scale default is set as TRUE
  if(is.null(par.extend$scale)){
    is.scale = TRUE
  } else {
    is.scale = par.extend$scale
  }
  
  namelist.par.extend = fulllist.par[-which(fulllist.par == 'sigma.b0.ss')]
  
  #a list records the number of beta in both data.full and data.mask needed for each parameter
  param_distribution = vector('list', length(fulllist.par))
  names(param_distribution) = fulllist.par
  
  

  if(!is.null(par.extend)){
    ########################################################################################################
    #key component begin
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input.")
    }
    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided.")
    }
    
    name.extend.par = names(par.extend$model)
    if(!all(name.extend.par %in% namelist.par.extend)){
      stop('one or more parameters indicated in "model" are not supported to be extended currently.')
    }
    
    #check 'data' component first
    
    #we could include 'session', 'trap', 'mask' levels covariates,
    #although some parameters such as 'mu' is only extendable under 'animal_ID'
    #level covariates in reality, but we do not make any restriction here
    
    input_data = par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap', 
                                     #'animal_ID', 
                                     'mask'))){
      stop("only 'session', 'trap', or 'mask' level data could be used as input.")
    }
    
    data.par.non.mask = data.frame(session = 1:n.sessions)
    
    #check 'session' level data
    #a data.frame with columns of 'session' and session level covariates
    df.s = input_data$session
    if(!is.null(df.s)){
      stopifnot(is(df.s, 'data.frame'))
      name.s = colnames(df.s)
      stopifnot('session' %in% name.s)
      stopifnot(length(name.s) > 1)
      if(any(duplicated(df.s$session))) stop('duplicated "session" input.')
      stopifnot(df.s$session %in% unique(data.full$session))
      var.s = name.s[-which(name.s == 'session')]
      data.par.non.mask = merge(data.par.non.mask, df.s, by = 'session', all = TRUE)
    } else {
      var.s = NULL
    }
    
    #check 'trap' level data
    df.t = input_data$trap
    if(!is.null(df.t)){
      df.t = extend_dat_check(df.t, 'trap', data.full, n.sessions, n.traps, identical_traps)
      name.t = colnames(df.t)

      var.t = name.t[-which(name.t == 'session'| name.t == 'trap')]
      data.par.non.mask = merge(data.par.non.mask, df.t, by = 'session', all = TRUE)
    } else {
      var.t = NULL
    }
    
    
    if(any(c('x', 'y') %in% c(var.s, var.t))){
      stop("'x' and 'y' are reserved column names, please choose another name.")
    }
    
    #check 'mask' level data
    
    #firstly check whether there is any formula in "model" contains "x" or "y"
    is_x = FALSE
    is_y = FALSE
    for(i in name.extend.par){
      foo = stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
      foo_var = all.vars(foo[[3]])
      if('x' %in% foo_var) is_x = TRUE
      if('y' %in% foo_var) is_y = TRUE
    }
    
    df.m = input_data$mask
    if(!is.null(df.m)){
      df.m = extend_dat_check(df.m, 'mask', data.mask, n.sessions, n.masks, identical_masks)
      
      #when both old argument 'ihd.opts' and par.extend$data$mask are provided,
      #we need to merge them together after confirming they have no identical covariates
      if(ihd_merge){
        name1 = colnames(df.m)
        name1 = name1[-which(name1 == 'session'| name1 == 'mask')]
        name2 = colnames(ihd_data)
        name2 = name2[-which(name2 == 'session'| name2 == 'mask')]
        
        if(any(duplicated(c(name1, name2)))){
          stop('One or more covariates are provided in both "ihd.opts" and "par.extend".')
        }
        
        df.m = merge(df.m, ihd_data, by = c('session', 'mask'))
      }
      
      name.m = colnames(df.m)
      var.m = name.m[-which(name.m == 'session'| name.m == 'mask')]
      data.par.mask = df.m
      #if 'x' or 'y' is in formula but not provided in par.extend$data$mask, then
      #take values from data.mask
      if(is_x & !('x' %in% var.m)){
        tem = data.mask[, c('session', 'mask', 'x')]
        data.par.mask = merge(data.par.mask, tem, by = c('session', 'mask'))
        var.m = c(var.m, 'x')
      }
      
      if(is_y & !('y' %in% var.m)){
        tem = data.mask[, c('session', 'mask', 'y')]
        data.par.mask = merge(data.par.mask, tem, by = c('session', 'mask'))
        var.m = c(var.m, 'y')
      }
      
    } else {
      if(is_x | is_y){
        data.par.mask = data.mask[,c('session', 'mask', 'x'[is_x], 'y'[is_y])]
        var.m = c('x'[is_x], 'y'[is_y])
      } else {
        var.m = NULL
        #data.par.mask is necessary anyway, create an empty data frame
        data.par.mask = data.frame(session = rep(1:n.sessions, n.masks),
                                   mask = do.call('c', lapply(as.list(n.masks), seq)))
      }
    }
    
    #check whether there is duplicated column names across these  data sets
    
    var.ex = c(var.s, var.t, var.m)
    
    if(any(duplicated(var.ex))) {
      stop('duplicated column name across different level data sets.')
    }
    
    #scale.covs needs to deal with each var.ex separately
    #record their mean and std
    var.ex.info = vector('list', length(var.ex))
    names(var.ex.info) = var.ex
    for(i in var.ex){
      if(!(i %in% var.m)){
        tem = data.par.non.mask[, i]
      } else {
        tem = data.par.mask[, i]
      }
      if(is.scale & is.numeric(tem)) var.ex.info[[i]] = c(mean(tem), sd(tem))
    }
    
    scale.covs = scale.closure(var.ex.info)
    
    gam_output = vector('list', length = length(name.extend.par))
    names(gam_output) = name.extend.par
    #check 'model'
    
    #deal the extend parameters one by one
    for(i in name.extend.par){
      foo = stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
      foo_var = all.vars(foo[[3]])
      if(length(foo_var) == 0) stop('please input a valid formula except "~ 1"')
      if(!all(foo_var %in% var.ex)) stop(paste0('one or more specifed covariates in "model"',
                                                ' is not provided in "data".'))
      
      #sigma_toa is not trap extendable
      if(i == 'sigma_toa' & any(foo_var %in% var.t)) stop("'sigma_toa' is not extendable in 'trap' level.")
      
      #in animal_ID embedded model, the parameter 'mu' is not extendable in 'trap' level either.
      if(i == 'mu' & any(foo_var %in% var.t)) stop("'mu' is not extendable in 'trap' level.")
      #in fact there should be much more constraint such as "mu" should not be linked to any 'trap' level covariates.
      #however, from the perspective of code, there is no necessary to make these constraint since all parameters
      #are treated as the same since the data structure does not change. However, 'sigma_toa' is different, because 'toa_ssq'
      #is stored under the data frame of "data.ID_matrix" which is not involved with any parameter extension stuff, so
      #there might be a concern about the work-flow later if we enable parameter 'sigma_toa' trap level extension.
      
      
      tem = formula_separate(foo, var.m)
      
      #"full" here does not mean full information, it does not contain mask information
      
      foo_full = tem$foo_full
      foo_mask = tem$foo_mask
      
      foo_full_var = all.vars(foo_full[[3]])
      foo_mask_var = all.vars(foo_mask[[3]])
      
      gam_output[[i]] = vector('list', 2)
      names(gam_output[[i]]) = c('gam_non_mask', 'gam_mask')
      
      if(length(foo_full_var) > 0){
        #avoid any NA
        tem = as.data.frame(lapply(data.par.non.mask[, foo_full_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index = apply(tem, 1, any)
        if(sum(index) > 0){
          stop(paste0("NA appears in covariates for parameter ", i, ", a review on the data is recommended."))
        }
        
        
        tem.data = data.par.non.mask
        tem.data[['gam.resp']] = 1
        
        for(j in foo_full_var){
          if(is.scale){
            if(is(tem.data[[j]], 'numeric')){
              mu = mean(tem.data[[j]])
              std = stats::sd(tem.data[[j]])
              if(std == 0) stop(paste0("input of ", j, " is constant."))
              tem.data[[j]] = (tem.data[[j]] - mu)/std
            }
          }
        }
        
        tem_model = mgcv::gam(foo_full, data = tem.data, fit = FALSE)
        gam_output[[i]][['gam_non_mask']] = tem_model
        design.matrix = tem_model$X
        n_col_full = ncol(design.matrix)
        design.matrix = as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) = paste(i, "_", colnames(design.matrix), sep = " ")
        design.matrix[['session']] = tem.data[['session']]
        if(!is.null(tem.data[['trap']])) design.matrix[['trap']] = tem.data[['trap']]
        if(!is.null(tem.data[['animal_ID']])) design.matrix[['animal_ID']] = tem.data[['animal_ID']]
        data.par.non.mask = merge(data.par.non.mask, design.matrix,
                                  by = c('session', 'trap'[!is.null(design.matrix[['trap']])],
                                         'animal_ID'[!is.null(design.matrix[['animal_ID']])]), all = TRUE)
        
      } else {
        #no matter whether this parameter is extended under other levels, the intercept always
        #in this data set. On the other hand, no matter it is extended under 'mask' level,
        #there is no intercept column in data.mask
        data.par.non.mask[[paste0(i, ' _ (Intercept)')]] = 1
        n_col_full = 1
      }
      
      
      if(length(foo_mask_var) > 0){
        tem = as.data.frame(lapply(data.par.mask[, foo_mask_var, drop = FALSE], is.na), stringsAsFactors = FALSE)
        index = apply(tem, 1, any)
        if(sum(index) > 0){
          stop(paste0("NA appears in covariates for parameter ", i, ", a review on the data is recommended."))
        }
        
        tem.data = data.par.mask
        tem.data[['gam.resp']] = 1
        
        for(j in foo_mask_var){
          if(is.scale){
            if(is(tem.data[[j]], 'numeric')){
              mu = mean(tem.data[[j]])
              std = stats::sd(tem.data[[j]])
              if(std == 0) stop(paste0("input of ", j, " is constant."))
              tem.data[[j]] = (tem.data[[j]] - mu)/std
            }
          }
        }
        
        #delete the first column, which is intercept
        tem_model = mgcv::gam(foo_mask, data = tem.data, fit = FALSE)
        gam_output[[i]][['gam_mask']] = tem_model
        if(i == 'D') fgam = tem_model
        design.matrix = tem_model$X[, -1, drop = FALSE]
        n_col_mask = ncol(design.matrix)
        design.matrix = as.data.frame(design.matrix, stringsAsFactors = FALSE)
        colnames(design.matrix) = paste(i, "_", colnames(design.matrix), sep = " ")
        design.matrix[['session']] = tem.data[['session']]
        design.matrix[['mask']] = tem.data[['mask']]
        data.par.mask = merge(data.par.mask, design.matrix,
                              by = c('session', 'mask'), all = TRUE)
      } else {
        n_col_mask = 0
      }
      
      param_distribution[[i]] = c(n_col_full, n_col_mask)
    }
    
    
    
    
    #the design matrix generated from gam() has been merged to data, the original
    #columns become redundant
    var.ex.col.index = which(colnames(data.par.non.mask) %in% var.ex)
    if(length(var.ex.col.index) > 0){
      data.par.non.mask = data.par.non.mask[, -var.ex.col.index, drop = FALSE]
    }
    
    var.ex.col.index = which(colnames(data.par.mask) %in% var.ex)
    if(length(var.ex.col.index) > 0){
      data.par.mask = data.par.mask[, -var.ex.col.index, drop = FALSE]
    }
    
    
    data.full = merge(data.full, data.par.non.mask, by = c('session',
                                                           'trap'[!is.null(data.par.non.mask[['trap']])],
                                                           'animal_ID'[!is.null(data.par.non.mask[['animal_ID']])],
                                                           'ID'[!is.null(data.par.non.mask[['ID']])]), all = TRUE)
    data.mask = merge(data.mask, data.par.mask, by = c('session', 'mask'))
    
    for(i in fulllist.par){
      if(!(i %in% name.extend.par)){
        data.full[[paste0(i, ' _ (Intercept)')]] = 1
        param_distribution[[i]] = c(1, 0)
      }
    }
    #key component end
    ########################################################################################################
    
  } else {
    for(i in fulllist.par){
      data.full[[paste0(i, ' _ (Intercept)')]] = 1
      param_distribution[[i]] = c(1, 0)
    }
    name.extend.par = NULL
    gam_output = NULL
    #if nothing extended, scale.covs just return the input
    scale.covs = function(covariates) return(covariates)
  }
  
  #finally, record the link, here we only check whether the name is in the full list of parameters,
  #so technically, user could assign link function to any parameter here, but not recommended
  df.par = data.frame(par = character(0), link = character(0), n_col_full = numeric(0),
                      n_col_mask = numeric(0), stringsAsFactors = FALSE)
  
  link = par.extend$link
  
  if(!is.null(link)){
    if(!all(names(link) %in% fulllist.par)){
      warning("one or more parameters denoted in 'link' are ignored as they are invalid.")
    }
  }
  
  j = 1
  for(i in fulllist.par){
    df.par[j, 'par'] = i
    df.par[j, 'n_col_full'] = param_distribution[[i]][1]
    df.par[j, 'n_col_mask'] = param_distribution[[i]][2]
    if(i %in% names(link)){
      #if 'link' is provided, check it whether is a valid input
      #currently we only support three link functions
      stopifnot(link[[j]] %in% c('log', 'logit', 'identity'))
      
      df.par[j, 'link'] = link[[i]]
    } else {
      #set default 'link' if it is not provided
      df.par[j, 'link'] = default.link(i)
    }
    j = j + 1
  }
  
  #D only accepts log as its link function
  stopifnot(df.par[which(df.par$par == 'D'), 'link'] == 'log')
  
  data.mask = sort.data(data.mask, "data.mask")
  
  #fgam only be valid when 'D' is extended on mask level
  #this component is for compatibility with ADMB model
  if(is.null(fgam)){
    tem = data.frame(gam.resp = rep(1, sum(n.masks)))
    fgam = mgcv::gam(gam.resp~1, data = tem, fit = FALSE)
  }
  
  
  return(list(data.full = data.full, data.mask = data.mask, data.par = df.par, is.scale = is.scale,
              name.extend.par = name.extend.par, fgam = fgam, gam_output= gam_output,
              scale.covs = scale.covs))
}



#############################################################################

ss.fun = function(ss.opts, data.full, data.ID_mask, animal.model, dims, bucket_info, sv, fix){
  cutoff <- ss.opts[["cutoff"]]
  ss.link <- ss.opts[["ss.link"]]
  directional <- ss.opts[["directional"]]
  het.source <- ss.opts[["het.source"]]
  het.source.method <- ss.opts[["het.source.method"]]
  n.dir.quadpoints <- ss.opts[["n.dir.quadpoints"]]
  n.het.source.quadpoints <- ss.opts[["n.het.source.quadpoints"]]
  
  
  if("ss" %in% bucket_info){
    if(is.null(ss.opts)) stop("Argument 'ss.opts' is missing.")
    if(is.null(cutoff)) stop("The 'cutoff' component of 'ss.opts' must be specified.")
    ## Warning for unexpected component names.
    if (!all(names(ss.opts) %in% c("cutoff", "het.source", "het.source.method",
                                   "n.het.source.quadpoints", "directional",
                                   "n.dir.quadpoints", "ss.link"))){
      warning("Components of 'ss.opts' may only consist of \"cutoff\", \"het.source\",
            \"het.source.method\", \"n.het.source.quadpoints\", \"directional\",
            \"n.dir.quadpoints\", \"ss.link\";
            others are being ignored.")
    }
    
    ## Setting default values for het.source and directional, ss.link.
    
    ## By default, directional calling model is only used if b2.ss appears in sv or fix.
    if (is.null(directional)){
      if (is.null(sv$b2.ss) & is.null(fix$b2.ss)){
        directional <- FALSE
      } else {
        directional <- TRUE
      }
    }
    
    #make detailed and final judgment about ss.dir, and settle some basic settings
    if (!directional){
      if (!is.null(sv$b2.ss) | !is.null(fix$b2.ss)){
        warning("Since the 'directional' component of 'ss.opts' is FALSE,
        the values of parameter b2.ss in 'sv' and 'fix' are being ignored")
      }
      n.dir.quadpoints = NULL
    } else {
      if(!is.null(fix$b2.ss)){
        if(fix$b2.ss != 0) {
          bucket_info = c(bucket_info, "ss.dir")
        } else {
          warning("'fix$b2.ss' is zero,
                  all corresponding parameters of 'directional' are ignored")
          directional = FALSE
          n.dir.quadpoints = NULL
        }
      } else {
        bucket_info = c(bucket_info, "ss.dir")
      }
    }
    
    #after settling down all basic info of ss.dir, set default n.dir.quadpoints
    if (directional){
      if (is.null(n.dir.quadpoints)){
        n.dir.quadpoints <- 8
      }
    } else {
      n.dir.quadpoints <- 1
    }
    
    ## By default, heterogeneity source strength model is only
    ## used if sigma.b0.ss appears in sv or fix, or het.source.method is specified
    if (is.null(het.source)){
      if (all(is.null(sv$sigma.b0.ss),
              is.null(fix$sigma.b0.ss),
              is.null(het.source.method))){
        het.source <- FALSE
      } else {
        het.source <- TRUE
      }
    }
    
    #make detailed and final judgment about ss.het, and settle some basic settings
    if (het.source){
      if (is.null(het.source.method)){
        het.source.method <- "GH"
      } else {
        if(het.source.method %in% c("GH", "rect")) stop(paste0("Only 'GH' and 'rect' are ",
                                                               "supported as het.source.method."))
      }
      if(!is.null(fix$sigma.b0.ss)){
        if(fix$sigma.b0.ss != 0) {
          bucket_info = c(bucket_info, "ss.het")
        } else {
          warning("'fix$sigma.b0.ss' is zero,
                  all corresponding parameters of 'het.source' are ignored")
          het.source = FALSE
          het.source.method = NULL
          n.het.source.quadpoints = NULL
        }
      } else {
        bucket_info = c(bucket_info, "ss.het")
      }
    } else {
      ## Fixing sigma.b0.ss to 0 if a heterogeneous source
      ## strength model is not being used.
      if (!is.null(sv$sigma.b0.ss) | !is.null(fix$sigma.b0.ss)){
        warning("As the 'het.source' component of 'ss.opts' is FALSE,
        the values of the parameter sigma.b0.ss in 'sv' and 'fix' are being ignored")
      }
      het.source.method = NULL
      n.het.source.quadpoints = NULL
    }
    
    #after settling down basic settings, assign default value to n.het...
    if(het.source){
      if(is.null(n.het.source.quadpoints)){
        n.het.source.quadpoints <- 15
      }
    } else {
      n.het.source.quadpoints <- 1
    }
    
    
    het.source.nodes <- 0
    het.source.weights <- 0
    if(het.source){
      if(het.source.method == "GH"){
        GHd <- fastGHQuad::gaussHermiteData(n.het.source.quadpoints)
        het.source.nodes <- GHd$x
        het.source.weights <- GHd$w
      }
    }
    
    
    ##################################################################
    #simple check and simple argument adjustment
    if(all(c("ss.dir", "ss.het") %in% bucket_info)){
      stop("Fitting of models with both heterogeneity in source signal strength
       and directional calling is not yet implemented.")
    }
    
    #set default for ss.link
    if (is.null(ss.link)){
      ss.link <- "identity"
    } else if (!(ss.link %in% c("identity", "log", "spherical"))){
      stop("Component 'ss.link' in 'ss.opts' must be \"identity\", \"log\", or \"spherical\".")
    }
    
    ## Error thrown for model with heterogeneity in source strength and a log-link function.
    if ("ss.het" %in% bucket_info){
      if (ss.link != "identity"){
        stop("Fitting of signal strength models with heterogeneity in source
         signal strength is only implemented with an identity link function.")
      }
    }
    
    if('ss.dir' %in% bucket_info){
      if(!('toa' %in% bucket_info)){
        stop('"toa" is required for directional signal strength model.')
      }
    }
    
    
    
    data.ss = subset(data.full, !is.na(ID))
    data.no.det.session = subset(data.full, is.na(ID))
    
    
    #adjust all capture history to zero if "ss" < cutoff, and then remove the observations if
    #there is no capture history for that "ID"
    rem = data.ss$ss < cutoff
    for(i in c('bincapt', 'ss', 'bearing', 'toa', 'dist')){
      data.ss[[i]] = ifelse(rem, 0, data.ss[[i]])
    }
    
    keep = agg_sort(data.ss, 'bincapt', c('session', 'animal_ID'[animal.model], 'ID'), sum)
    
    keep = subset(keep, x > 0)
    if(nrow(keep) == 0) stop("None of signal received is greater than the cut off threshold.")
    colnames(keep)[which(colnames(keep) == 'x')] = "keep"
    data.ss = merge(data.ss, keep, by = c("session", "animal_ID"[animal.model], "ID"), all = TRUE)
    data.ss = subset(data.ss, !is.na(data.ss$keep))
    
    data.ss = data.ss[,-which(colnames(data.ss) == "keep")]
    data.ss = sort.data(data.ss, "data.full")
    data.full = rbind(data.ss, data.no.det.session)
    data.full = sort.data(data.full, "data.full")
    
    #refresh "dims"
    if(animal.model){
      new.n.animals = agg_sort(data.ss, 'animal_ID', 'session', function(x) length(unique(x)))
      new.n.IDs = agg_sort(data.ss, c('animal_ID', 'ID'), 'session', function(x) length(unique(x)))
      new.n.animal.call = agg_sort(data.ss, 'ID', c('session', 'animal_ID'), function(x) length(unique(x)))

      dims$n.animals[new.n.animals$session] = new.n.animals$x
      dims$n.animal.call = new.n.animal.call$x
    } else {
      new.n.IDs = agg_sort(data.ss, 'ID', 'session', function(x) length(unique(x)))
    }
    
    n.removed = sum(dims$n.IDs) - sum(new.n.IDs$x)
    if(n.removed > 0){
      message(n.removed, " capture history entries have no received signal",
              " strengths above the cutoff and have therefore been removed.\n",
              sep = "")
      #a logical value to record whether there is any ID been removed
      ID.removed = TRUE
    } else ID.removed = FALSE
    
    dims$n.IDs[new.n.IDs$session] = new.n.IDs$x
    
    if(ID.removed){
      #when any ID is deleted due to the cut off of ss
      if("animal_ID" %in% colnames(data.full)){
        data.ID_mask = subset(data.ID_mask, paste(session, animal_ID, ID, sep = "-") %in%
                                unique(paste(data.full$session, data.full$animal_ID,data.full$ID, sep = "-")))
      } else {
        data.ID_mask = subset(data.ID_mask, paste(session, ID, sep = "-") %in%
                                unique(paste(data.full$session, data.full$ID, sep = "-")))
      }
      data.ID_mask = sort.data(data.ID_mask, "data.ID_mask")
    }
    
    #modified ss.opts, contains much condensed and regular information
    ss.opts = list(het.source.method = het.source.method,
                   n.dir.quadpoints = n.dir.quadpoints,
                   n.het.source.quadpoints = n.het.source.quadpoints,
                   het.source.nodes = het.source.nodes,
                   het.source.weights = het.source.weights,
                   ss.link = ss.link,
                   cutoff = cutoff)
    
    
  } else {
    if (!is.null(ss.opts)){
      warning("Argument 'ss.opts' is being ignored as signal strength is not provided.")
    }
    
    ss.opts <- NULL
    
    if(any(!is.null(sv$b0.ss),
           !is.null(sv$b1.ss),
           !is.null(sv$b2.ss),
           !is.null(sv$sigma.ss),
           !is.null(sv$sigma.b0.ss),
           !is.null(fix$b0.ss),
           !is.null(fix$b1.ss),
           !is.null(fix$b2.ss),
           !is.null(fix$sigma.ss),
           !is.null(fix$sigma.b0.ss))){
      warning(paste0("Not 'ss' model, coresponding parameters specified in",
                     "'sv' or 'fix' will be ignored"))
    }
  }
  
  
  
  return(list(data.full = data.full, dims = dims, data.ID_mask = data.ID_mask,
              ss.opts = ss.opts, bucket_info = bucket_info))
}

############################################################################
CR_SL = function(cue.rates, survey.length, bucket_info, dims){
  n.sessions = dims$n.sessions
  
  if (is.null(survey.length)){
    survey.length <- rep(1, n.sessions)
    if (!is.null(cue.rates)){
      stop("The use of `cue.rates' without `survey.length' is no longer supported.
           Please provide `survey.length', and ensure `cue.rates' is measured in the same time units.")
    }
  } else {
    if (length(survey.length) != n.sessions){
      stop("The argument `survey.length' must have a value for each session.")
    }
  }
  
  if (!is.null(cue.rates)){
    bucket_info = c(bucket_info, "freqs")
    mean.cue.rates <- mean(cue.rates)
  } else {
    mean.cue.rates = 1
  }
  
  return(list(bucket_info = bucket_info, mean.cue.rates = mean.cue.rates, survey.length = survey.length))
  
}


############################################################################

param.detfn.fun = function(animal.model, sv, fix, bounds, name.extend.par, detfn, data.full, data.mask,
                           data.par, ss.opts, bucket_info, fulllist.par, A, buffer,
                           survey.length, dims){
  #check bounds, sv, fix
  if(!is.null(bounds)){
    if(!is(bounds, 'list')) stop ('"bounds" must be a list.')
    if(!all(names(bounds) %in% fulllist.par)) stop('one or more element in "bounds" is not valid.')
    if(!all(sapply(bounds, function(x) nrow(as.matrix(x))) == 2)){
      stop('each element of "bounds" should be a vector with two values, or a matrix with two rows.')
    }
  } else {
    bounds = vector('list', 0)
  }
  
  if(!is.null(sv)){
    if(!is(sv, 'list')) stop('"sv" must be a list.')
    if(!all(names(sv) %in% fulllist.par)) stop('one or more element in "sv" is not valid.')
  } else {
    sv = vector('list', 0)
  }
  
  if(!is.null(fix)){
    if(!is(fix, 'list')) stop('"fix" must be a list.')
    if(!all(names(fix) %in% fulllist.par)) stop('one or more element in "fix" is not valid.')
  } else {
    fix = vector('list', 0)
  }
  
  
  
  #detfn is the basis of choosing parameters, settle this first
  
  if(is.null(detfn)){
    #the detfn in 'ss' is a little different from original, may be needed
    #to be restored as original at very end of 'fit' function
    #ss_het only support ss.link == 'identity', but this has been checked in ss.fun
    if("ss" %in% bucket_info){
      if("ss.het" %in% bucket_info){
        detfn = 'ss_het'
      } else if("ss.dir" %in% bucket_info){
        detfn = 'ss_dir'
      } else {
        detfn = 'ss'
      }
    } else {
      detfn = "hn"
    }
  } else {
    if("ss" %in% bucket_info){
      if(detfn != 'ss'){
        warning('signal strenth is provided, therefore the specified "detfn" is ignored')
      }
      
      if("ss.het" %in% bucket_info){
        detfn = 'ss_het'
      } else if("ss.dir" %in% bucket_info){
        detfn = 'ss_dir'
      } else {
        detfn = 'ss'
      }
      
    } else {
      if(detfn == 'ss'){
        warning('signal strenth is not provided, therefore the specified "detfn" is replaced to "hn" (half normal),')
        detfn = 'hn'
      } else {
        if(!detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth')){
          warning('invalid "detfn" is specified, the default will be used instead')
          detfn = 'hn'
        }
      }
    }
  }
  
  #then, decide the original names of the parameters will be input into TMB
  param.og = detfn.params(detfn)

  
  if('bearing' %in% bucket_info){
    param.og = c(param.og, 'kappa')
  }
  
  if('dist' %in% bucket_info){
    param.og = c(param.og, 'alpha')
  }
  
  if('toa' %in% bucket_info){
    param.og = c(param.og, 'sigma.toa')
  }
  
  param.og = c(param.og, "D")
  
  if(animal.model) param.og = c(param.og, "mu")
  
  #generate all DATA_MATRIX for every each parameter
  #and in the loop, we also do some basic conflict checking
  data.full = sort.data(data.full, "data.full")
  
  
  design.matrices.full = vector('list', length(fulllist.par))
  names(design.matrices.full) = fulllist.par
  
  tem = strsplit(colnames(data.full), " _ ")
  tem = lapply(tem, function(x) x[[1]])
  clean.names.full = do.call('c', tem)
  
  design.matrices.mask = vector('list', length(fulllist.par))
  names(design.matrices.mask) = fulllist.par
  
  tem = strsplit(colnames(data.mask), " _ ")
  tem = lapply(tem, function(x) x[[1]])
  clean.names.mask = do.call('c', tem)
  
  for(i in fulllist.par){
    index_par = which(data.par$par == i)
    
    if(data.par[index_par, "n_col_mask"] > 0){
      design.matrices.mask[[i]] = as.matrix(data.mask[, which(clean.names.mask == i), drop = FALSE])
    } else {
      design.matrices.mask[[i]] = matrix(0, ncol = 2, nrow = 2)
    }
    
    #because intercept is in "data.full", so this design matrix always has someting
    design.matrices.full[[i]] = as.matrix(data.full[, which(clean.names.full == i), drop = FALSE])
    
    if(i %in% param.og){
      #deal with "fix" first
      if(i %in% name.extend.par){
        #a valid extended parameter should not be fixed
        if(!is.null(fix[[i]])){
          warning(paste0("parameter ", i, " is going to be extended, therefore, it cannot be fixed."))
          fix[[i]] = NULL
        }
      } else {
        #if this valid parameter is not expended, check whether its fixed value is in the range
        if(!is.null(fix[[i]])){
          tem = param.range.validate(i, fix[[i]])
          if(!tem) stop(paste0("fixed value for parameter ", i, " is out of valid range."))
        }
      }
      #then deal with "sv" and "bounds"
      #if a parameter is extended, the start value will be its intercept beta's start value
      #the bounds of an extended parameter will be assigned to its intercept as well
      if(!is.null(bounds[[i]])){
        #the lower and upper bounds all should be in the valid range of this parameter
        tem1 = param.range.validate(i, ifelse(bounds[[i]][1] == 0, 1e-15, bounds[[i]][1]))
        tem2 = param.range.validate(i, ifelse(bounds[[i]][2] == 0, 1e-15, bounds[[i]][2]))
        if(!all(tem1, tem2)) stop(paste0("bounds for parameter ", i, " is out of valid range."))
        if(bounds[[i]][2] <= bounds[[i]][1]) stop(paste0("invalid range for parameter ", i))
        
        if(!is.null(sv[[i]])){
          tem = param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        } else {
          #if sv is not assigned, assign a default value
          sv[[i]] = default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts,
                                              detfn = detfn, dims = dims, animal.model = animal.model,
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
          
          #since the bound is provided, we need the default start value is located within the bound
          #if the default start value is out of the specified bounds, re-assign the start value
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            sv[[i]] = 0.5 * (bounds[[i]][1] + bounds[[i]][2])
          }
        }
        
        #whether the fixed value is inside the assigned bound
        #since we have check the validation of fixed value, no need to check it here like sv
        #if bound is NULL, no need to check fixed value against the defaulted assigned bound
        if(!is.null(fix[[i]])){
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            stop(paste0("the start value of ", i, " is out of the specified corresponding bounds."))
          }
        }
        
      } else {
        bounds[[i]] = default.bounds(i)
        
        if(!is.null(sv[[i]])){
          tem = param.range.validate(i, sv[[i]])
          if(!tem) stop(paste0("start value for parameter ", i, " is out of valid range."))
          #I'm not sure if the assigned start value be out of the default bound, what we should do
          #discuss with Ben later.
          if(!all(sv[[i]] >= bounds[[i]][1], sv[[i]] <= bounds[[i]][2])){
            warning(paste0("specified start value of '", i, "' is extreme, a carefull checking ",
                           "is recommended."))
            bounds[[i]][1] = -2 * abs(sv[[i]])
            bounds[[i]][2] = 2 * abs(sv[[i]])
          }
          
        } else {
          #if sv is not assigned, assign a default value
          
          sv[[i]] = default.sv(i, info = list(data.full = data.full, data.mask = data.mask,
                                              buffer = buffer, A = A, ss.opts = ss.opts,
                                              detfn = detfn, dims = dims, animal.model = animal.model,
                                              survey.length = survey.length, param.og = param.og,
                                              sv = sv))
          
        }
      }
    } else {
      #a start value is needed for every parameter regardless it is used or not
      if(!is.null(sv[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'sv' will be ignored because ",
                       "the detect function is incompatible"))
      }
      sv[[i]] = 0.5
      
      #if this parameter is not used in the model and in the "fix", remove it from "fix"
      #by assigning NULL to it
      if(!is.null(fix[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'fix' will be ignored because ",
                       "the detect function is incompatible"))
        fix[[i]] = NULL
      }
      
      if(!is.null(bounds[[i]])){
        warning(paste0("parameter ", i, " assigned in argument 'bounds' will be ignored because ",
                       "the detect function is incompatible"))
      }
      
      bounds[[i]] = c(0.3, 0.7)
      
      #if this not used parameter specified to be extended, everything refer to extension will be ignored
      if(i %in% name.extend.par){
        warning(paste0("parameter ", i, " specified to be extended, however, this parameter is not used ",
                       "in the model, everything about extension will be ignored"))
        #correct all of the relevent information
        design.matrices.full[[i]] = design.matrices.full[[i]][, 1, drop = FALSE]
        design.matrices.mask[[i]] = matrix(0, ncol = 2, nrow = 2)
        data.par[index_par, 'n_col_mask'] = 0
        data.par[index_par, 'n_col_full'] = 1
        name.extend.par = name.extend.par[-which(name.extend.par == i)]
      }
      
      
    }
  }
  
  #construct the parameters input, which must be transformed
  #they are not needed to be named, especially the betas for extended parameters
  
  #the sv.input contains vectors as element in the list
  #in the vector, the first 'n_col_full' numbers are for the relevant columns in
  #data.full, and the rest for the relevant columns in data.mask
  #the fit.input contains scalars as element in the list, because extended param
  #has no "fit"
  sv.input = vector('list', length(fulllist.par))
  names(sv.input) = fulllist.par
  fix.input = vector('list', length(fulllist.par))
  names(fix.input) = fulllist.par
  #the bounds.input contains matrices as each element in the list
  bounds.input = vector('list', length(fulllist.par))
  names(bounds.input) = fulllist.par
  
  
  for(i in fulllist.par){
    index_par = which(data.par$par == i)
    length_param = data.par[index_par, 'n_col_full'] + data.par[index_par, 'n_col_mask']
    link = data.par[index_par, 'link']
    
    
    #here is a small problem need attention, if the user specify a parameter to be fixed to zero,
    #and its link function is 'log', than, the real input of fixed value in TMB will be log(1e-20)
    #it is not exactly zero. But this problem could be solved by rounding the output to 7 digits,
    #so in the output, this parameter will be zero again.
    if(!is.null(fix[[i]])) {
      fix.input[[i]] = link.fun(link = link, value = fix[[i]])
      #if a parameter is fixed, we set its start value to its fixed value
      #since only not extended parameter could be fixed, so we don't need to
      #worry about its length, it is always a scaler
      sv.input[[i]] = fix.input[[i]]
    } else {
      sv.input[[i]] = numeric(length_param)
      #since the numeric(n) is defaults set to 0, we only need to assign the
      #first element, no matter it is extended or not
      sv.input[[i]][1] = link.fun(link = link, value = sv[[i]])
    }
    
    
    
    #set the default bounds of extended betas as the starter of our bounds matrix
    #and replace the first one (or the only one) by the assigned bounds' transformation
    if(link == 'identity'){
      bounds.input[[i]] = matrix(rep(c(-1e8, 1e8), length_param),
                                 nrow = 2)
    } else {
      bounds.input[[i]] = matrix(rep(c(-log(1e20), log(1e20)), length_param),
                                 nrow = 2)
    }
    
    bounds.input[[i]][1, 1] = link.fun(link = link, value = bounds[[i]][1])
    bounds.input[[i]][2, 1] = link.fun(link = link, value = bounds[[i]][2])
  }

  return(list(design.matrices.full = design.matrices.full,
              design.matrices.mask = design.matrices.mask,
              name.extend.par = name.extend.par,
              data.par = data.par,
              sv.input = sv.input,
              fix.input = fix.input,
              bounds.input = bounds.input,
              detfn = detfn,
              param.og = param.og))
  
}



outFUN = function(data.par, data.full, data.traps, data.mask, data.dists.thetas, detfn, param.og, param.og.4cpp, o, opt,
                  name.fixed.par, name.extend.par, dims, DX.full, DX.mask, fix.input, bucket_info, cue.rates, mean.cue.rates, A,
                  survey.length, sound.speed, par.extend, arg.input, fgam, gam_output, scale.covs, is.scale, ss.link, cutoff){
  ###################################################################################################################
  #sort out output for the function
  out = vector('list', 36)
  names(out) = c("fn", "coefficients", "coeflist", "se", "loglik", "maxgrad", "cor", "vcov", "npar", "npar_re",
                 "npar_sdrpt", "npar_rep", "npar_total", "hes", "eratio", "D.mask", "mm.ihd", "args", "n.sessions",
                 "fit.types", "infotypes", "detpars", "suppars", "D.betapars", "phases", "par.links", "par.unlinks",
                 "fit.ihd", "re.detfn", "fit.freqs", "first.calls", "scale.covs", "model.formula", "fgam", "all.covariates",
                 "output.tmb")
  #create output for TMB model
  out[['output.tmb']] = vector('list', 20)
  names(out[['output.tmb']]) = c('coef_link', 'se_link', 'DX', 'detfn', 'param.og', 'param.extend', 'param.fix',
                                 'param.info.table', 'data.traps', 'data.full', 'data.mask', 'data.dists.thetas', 'dims',
                                 'avg_cue_rates', 'sound.speed', 'area_unit', 'survey.length', 'ss.link', 'cutoff', 'gam_output')
  
  #give an index to each parameter, to make it easier to find it in "data.par"
  par.id = 1:nrow(data.par)
  names(par.id) = data.par$par
  
  param.with.esa = c(param.og, 'esa')
  param.og.4cpp = c(param.og.4cpp, 'esa')
  o_value = vector('list', length(param.with.esa))
  names(o_value) = param.with.esa
  o_sd = vector('list', length(param.with.esa))
  names(o_sd) = param.with.esa
  for(i in 1:length(param.with.esa)){
    o_value[[param.with.esa[i]]] = o$value[which(names(o$value) == param.og.4cpp[i])]
    o_sd[[param.with.esa[i]]] = o$sd[which(names(o$value) == param.og.4cpp[i])]
  }
  
  
  ############################################################################################################
  #the 1st component: "fn"
  out$fn = 'secr_tmb'
  
  
  #############################################################################################################
  #the 2nd component: "coefficients", there are "esa.[session]" not included yet
  #this one can be generated from "name.fited.pars", "name.fixed.par", "name.extend.par"
  
  #create a vector for parameter names used in various of output components
  #adjust it to be the same with the output of current ADMB model
  
  #here "name.fited.pars" is used instead of "param.with.esa" is to make sure the order
  #of parameters are the same as TMB model, and "esa" is not in "param.og"
  #if nothing goes wrong, 'name.fited.pars" shoule be exactly the same with "param.with.esa"
  name.fited.pars = unique(names(o$value))
  name.fited.pars = gsub("_", "\\.", name.fited.pars)
  name_output = vector('list', length(name.fited.pars))
  names(name_output) = name.fited.pars
  
  #names for tmb output, "esa" would be totally separated from other parameters
  name_output_tmb = name_output
  name_output_tmb$esa = NULL
  
  name_split_full = strsplit(colnames(data.full), " _ ")
  name_split_mask = strsplit(colnames(data.mask), " _ ")
  
  #out_coef is a temporary object, which will be assigned to out$coefficients later
  out_coef = vector('list', length(name.fited.pars))
  names(out_coef) = name.fited.pars
  #take 'esa' out from tmb output
  coef_link_tmb = out_coef
  coef_link_tmb$esa = NULL
  
  for(i in name.fited.pars){
    #esa is not in the original parameter list, deal with it separately
    if(i != 'esa'){
      #firstly, sort out names for output coefficients
      index_full = which(sapply(name_split_full, function(x) x[1] == i))
      index_mask = which(sapply(name_split_mask, function(x) x[1] == i))
      name_output[[i]] = c(colnames(data.full)[index_full], colnames(data.mask)[index_mask])
      name_output_tmb[[i]] = name_output[[i]]
      name_output[[i]] = gsub(" _ ", ".", name_output[[i]])
      #if i is "D", then keep the .(Intercept) no matter what model it is
      if(i != "D" & !i %in% name.extend.par){
        name_output[[i]] = gsub("\\.\\(Intercept\\)", "", name_output[[i]])
      }
      
      if(!i %in% name.fixed.par){
        #then, input the estimated values
        out_coef[[i]] = o_value[[i]]
        names(out_coef[[i]]) = paste(name_output[[i]], 'link', sep = "_")
        coef_link_tmb[[i]] = o_value[[i]]
        names(coef_link_tmb[[i]]) = name_output_tmb[[i]]
        
        if(i == "D"){
          #for "D", take it as extended no matter whether it actually is
          tem = out_coef[[i]]
        } else {
          if(!i %in% name.extend.par){
            tem = unlink.fun(data.par[par.id[i], 'link'], out_coef[[i]])
          } else {
            tem = out_coef[[i]]
          }
        }
        names(tem) = name_output[[i]]
        out_coef[[i]] = c(out_coef[[i]], tem)
        
        
        #for homo "D" models, there is a element in the vector named "D"
        if(i == 'D' & !i %in% name.extend.par){
          tem = unlink.fun(data.par[par.id[i], 'link'], tem)
          names(tem) = 'D'
          out_coef[[i]] = c(out_coef[[i]], tem)
        }
      } else {
        coef_link_tmb[[i]] = fix.input[[i]]
        names(coef_link_tmb[[i]]) = name_output_tmb[[i]]
      }
    } else {
      name_output[[i]] = paste('esa', 1:dims$n.sessions, sep = ".")
      out_coef[[i]] = o_value[[i]]
      names(out_coef[[i]]) = name_output[[i]]
    }
  }
  
  out$coefficients = do.call('c', out_coef)
  names(out$coefficients) = do.call('c', lapply(out_coef, function(x) names(x)))
  out[['output.tmb']][['coef_link']] = coef_link_tmb
  #######################################################################################################
  #the 3rd component: "coeflist"
  
  out$coeflist = vector('list', 0)
  
  D_pars = name_output$D
  det_pars = name.fited.pars[which(!name.fited.pars %in% c('D', 'kappa', 'alpha', 'sigma.toa', 'esa'))]
  sup_pars = name.fited.pars[which(name.fited.pars %in% c('kappa', 'alpha', 'sigma.toa'))]
  
  tem = val(o_value[['D']])
  for(i in 1:length(name_output$D)){
    tem_name = paste0('D_betapars_link[', i, ']')
    out$coeflist[[tem_name]] = tem[i]
  }
  
  for(i in 1:length(det_pars)){
    if(!det_pars[i] %in% name.fixed.par){
      tem = val(o_value[[det_pars[i]]])
      tem_name = paste0('detpars_link[', i, ']')
      out$coeflist[[tem_name]] = tem
    }
  }
  
  if(length(sup_pars) > 0){
    for(i in 1:length(sup_pars)){
      if(!sup_pars[i] %in% name.fixed.par){
        tem = val(o_value[[sup_pars[i]]])
        tem_name = paste0('suppars_link[', i, ']')
        out$coeflist[[tem_name]] = tem
      }
    }
  }
  
  
  for(i in 1:dims$n.sessions){
    tem = val(o_value[['esa']][i])
    tem_name = paste0('esa[', i, ']')
    out$coeflist[[tem_name]] = tem
  }
  
  #deal with "D" and "D.mask[s]"
  if(!'D' %in% name.extend.par){
    out$coeflist[['D']] = unlink.fun(data.par[par.id['D'], 'link'], val(o_value[['D']]))
    for(i in 1:dims$n.sessions){
      tem_name = paste0('D_mask[', i, ']')
      out$coeflist[[tem_name]] = rep(out$coeflist[['D']], dims$n.masks[i])
    }
  } else {
    tem_n_col = data.par[par.id['D'], c('n_col_full', 'n_col_mask')]
    tem_full = DX.full[['D']] %*% val(o_value[['D']])[1:tem_n_col[1, 1]]
    if(tem_n_col[1, 2] > 0){
      tem_mask = DX.mask[['D']] %*% val(o_value[['D']])[(tem_n_col[1, 1] + 1):sum(tem_n_col)]
    } else {
      tem_mask = rep(0, sum(dims$n.masks))
    }
    
    for(i in 1:dims$n.sessions){
      tem_name = paste0('D_mask[', i, ']')
      #'D' is only session not trap extendable, we could take the first observation in each session
      #only from 'data.full', which is tem_full below
      out$coeflist[[tem_name]] = tem_full[which(data.full$session ==i)[1]] +
        tem_mask[which(data.mask$session == i)]
      out$coeflist[[tem_name]] = unlink.fun(data.par[par.id['D'], 'link'],
                                            out$coeflist[[tem_name]])
    }
  }
  
  
  
  ####################################################################################################
  #the 4th component: "se"
  
  out_se = vector('list', length(name.fited.pars))
  names(out_se) = name.fited.pars
  
  #in tmb output, not include 'esa'
  se_tmb = out_se
  se_tmb$esa = NULL
  
  for(i in name.fited.pars){
    if(i != 'esa'){
      if(!i %in% name.fixed.par){
        #then, input the estimated values
        out_se[[i]] = o_sd[[i]]
        names(out_se[[i]]) = paste(name_output[[i]], 'link', sep = "_")
        
        se_tmb[[i]] = o_sd[[i]]
        names(se_tmb[[i]]) = name_output_tmb[[i]]
        if(i == "D"){
          #for "D", take it as extended no matter whether it actually is
          tem = out_se[[i]]
        } else {
          if(!i %in% name.extend.par){
            tem = delta.fun(data.par[par.id[i], 'link'], out_se[[i]], o_value[[i]])
          } else {
            tem = out_se[[i]]
          }
        }
        names(tem) = name_output[[i]]
        out_se[[i]] = c(out_se[[i]], tem)
        
        
        #for homo "D" models, there is a element in the vector named "D"
        if(i == 'D' & !i %in% name.extend.par){
          tem = delta.fun(data.par[par.id[i], 'link'], tem, o_value[['D']])
          names(tem) = 'D'
          out_se[[i]] = c(out_se[[i]], tem)
        }
      } else {
        se_tmb[[i]] = NA
        names(se_tmb[[i]]) = name_output_tmb[[i]]
      }
    } else {
      out_se[[i]] = o_sd[[i]]
      names(out_se[[i]]) = name_output[[i]]
    }
    
  }
  
  out$se = do.call('c', out_se)
  names(out$se) = do.call('c', lapply(out_se, function(x) names(x)))
  out[['output.tmb']][['se_link']] = se_tmb
  #####################################################################################################
  #insert the necessary TMB output
  tem = vector('list', 2)
  names(tem) = c('DX_non_mask', 'DX_mask')
  tem[['DX_non_mask']] = DX.full
  tem[['DX_mask']] = DX.mask
  
  out[['output.tmb']][['DX']] = tem
  
  tem = vector('list', 2)
  names(tem) = c('value', 'std')
  tem[['value']] = val(o_value$esa)
  tem[['std']] = o_sd$esa
  
  out[['esa']] = tem
  
  out[['output.tmb']][['param.og']] = param.og
  out[['output.tmb']][['param.extend']] = name.extend.par
  out[['output.tmb']][['param.fix']] = name.fixed.par
  out[['output.tmb']][['param.info.table']] = data.par
  out[['output.tmb']][['data.traps']] = data.traps
  out[['output.tmb']][['dims']] = dims
  out[['output.tmb']][['area_unit']] = A
  out[['output.tmb']][['survey.length']] = survey.length
  out[['output.tmb']][['avg_cue_rates']] = mean.cue.rates
  out[['output.tmb']][['data.full']] = data.full
  out[['output.tmb']][['data.mask']] = data.mask
  out[['output.tmb']][['data.dists.thetas']] = data.dists.thetas
  out[['output.tmb']][['detfn']] = detfn
  out[['output.tmb']][['ss.link']] = ss.link
  out[['output.tmb']][['cutoff']] = cutoff
  out[['output.tmb']][['sound.speed']] = sound.speed
  out[['output.tmb']][['gam_output']] = gam_output
  ######################################################################################################
  #the 5th component: "loglik"
  out$loglik = -1 * opt$objective
  
  #####################################################################################################
  #the 6th component: "maxgrad"
  out$maxgrad = max(abs(o$gradient.fixed))
  
  #####################################################################################################
  #the 7th and 8th component: "cor" and "vcov", after communication it is fine to
  #use the output from TMB directly.
  tem = name_output_tmb
  
  for(i in names(tem)){
    if(i %in% name.fixed.par) tem[i] = list(NULL)
  }
  tem = do.call('c', tem)
  tem = c(tem, paste('esa', seq(dims$n.sessions), sep = "."))
  #remove fixed parameters from cov matrix
  if(length(name.fixed.par) > 0){
    which_fixed = which(names(o$value) %in% name.fixed.par)
    out$vcov = o$cov[-which_fixed, -which_fixed]
  } else {
    out$vcov = o$cov
  }
  
  dimnames(out$vcov) = list(tem, tem)
  
  #convert the covariance matrix to correlation matrix
  tem = sqrt(diag(out$vcov))
  out$cor = out$vcov / outer(tem, tem)
  
  
  ######################################################################################################
  #The 9th component: "npar"
  
  out$npar = sum(sapply(name_output, length))
  
  ######################################################################################################
  #From the 10th component to the 14th component, not sure what they are, skip them temporarily
  
  
  
  
  ######################################################################################################
  #The 15th component: 'eratio'
  #obtain all estimated values excluding 'esa'
  tem = o$value[which(!names(o$value) %in% c('esa', name.fixed.par))]
  out$eratio = max(abs(o$gradient.fixed) / abs(tem))
  
  
  ######################################################################################################
  #The 16th component: "D.mask"
  out$D.mask = vector('list', dims$n.sessions)
  
  for(s in 1:dims$n.sessions){
    tem = paste0("D_mask[", s, "]")
    out$D.mask[[s]] = out$coeflist[[tem]]
    tem = stringr::str_pad(as.character(1:dims$n.masks[s]), nchar(as.character(dims$n.masks[s])),
                           pad = '0', side = 'left')
    names(out$D.mask[[s]]) = paste(paste0('D_mask[', s, ']'), tem, sep = ".")
  }
  
  
  
  ########################################################################################################
  #The 17th component: "mm.ihd"
  
  out$mm.ihd = vector('list', dims$n.sessions)
  
  for(s in 1:dims$n.sessions){
    #firstly got the first index of this session in data.full
    tem = min(which(data.full$session == s))
    #then, got the data from DX.full, which denotes "design matrix for data.full"
    tem = DX.full$D[tem,]
    #finally, generate a matrix with the same length of n.masks and contains the information
    #of "D" under "session - level"
    tem = matrix(rep(tem, dims$n.masks[s]), nrow = dims$n.masks[s], byrow = T)
    
    if(data.par[par.id['D'], 'n_col_mask'] > 0){
      out$mm.ihd[[s]] = cbind(tem, DX.mask$D[data.mask$session == s,])
    } else {
      out$mm.ihd[[s]] = tem
    }
    colnames(out$mm.ihd[[s]]) = gsub("D\\.", "", name_output$D)
    rownames(out$mm.ihd[[s]]) = NULL
  }
  
  
  ######################################################################################################
  #The 18th component: "args"
  
  out$args = arg.input
  
  ######################################################################################################
  #the 19th component: 'n.sessions'
  out$n.sessions = dims$n.sessions
  #####################################################################################################
  #the 20th component: 'fit.types'
  tem = logical(5)
  names(tem) = c('bearing', 'dist', 'ss', 'toa', 'mrds')
  
  tem['bearing'] <- ('bearing' %in% bucket_info)
  tem['dist'] <- ('dist' %in% bucket_info)
  tem['ss'] <- ('ss' %in% bucket_info)
  tem['toa'] <- ('toa' %in% bucket_info)
  tem['mrds'] <- ('mrds' %in% bucket_info)
  
  out$fit.types = tem
  
  #####################################################################################################
  #the 21st component: 'infotypes'
  
  out$infotypes = bucket_info[which(bucket_info %in% c('bearing', 'dist', 'ss', 'toa', 'mrds'))]
  
  
  #####################################################################################################
  #The 22nd component: 'detpars'
  out$detpars = param.og[which(!param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D'))]
  
  
  ######################################################################################################
  #the 23rd component: 'suppars'
  out$suppars = param.og[which(param.og %in% c('kappa', 'alpha', 'sigma.toa'))]
  if(length(out$suppars) == 0) out$suppars = NULL
  
  ######################################################################################################
  #the 24th component: 'D.betapars'
  out$D.betapars = name_output$D
  
  ######################################################################################################
  #the 25th component: 'phases'
  all.pars = c(out$detpars, out$suppars, out$D.betapars)
  
  out$phases = ifelse(all.pars %in% names(out$args$fix), -1, 1)
  names(out$phases) = all.pars
  out$phases = as.list(out$phases)
  
  #####################################################################################################
  #the 26th component: 'par.links'
  #the 27th component: 'par.unlinks'
  
  out$par.links = vector('list', length(all.pars))
  names(out$par.links) = all.pars
  
  out$par.unlinks = vector('list', length(all.pars))
  names(out$par.unlinks) = all.pars
  
  for(i in c(out$detpars, out$suppars)){
    if(!i %in% name.extend.par){
      if(data.par[par.id[i], 'link'] == 'log'){
        out$par.links[[i]] = cus_log
        out$par.unlinks[[i]] = cus_log_unlink
      } else if(data.par[par.id[i], 'link'] == 'logit'){
        out$par.links[[i]] = cus_logit
        out$par.unlinks[[i]] = cus_logit_unlink
      } else {
        out$par.links[[i]] = cus_identity
        out$par.unlinks[[i]] = cus_identity
      }
    } else {
      out$par.link[[i]] = cus_identity
      out$par.unlinks[[i]] = cus_identity
    }
  }
  
  for(i in out$D.betapars){
    out$par.links[[i]] = cus_identity
    out$par.unlinks[[i]] = cus_identity
  }
  
  
  ######################################################################################################
  #the 28th component: 'fit.ihd'
  #only if the model is 'D' extended and the extension is "marks" extension only,
  #in other words, the model which is equivalent to the incumbent "inhomo" model,
  #then this 'ihd' is set to TRUE
  
  if('D' %in% name.extend.par & data.par[par.id['D'], 'n_col_full'] == 1){
    out$fit.ihd = TRUE
  } else {
    out$fit.ihd = FALSE
  }
  
  
  #####################################################################################################
  #the 29th component: 're.detfn', not sure what this is
  #in all of the tested sample data, this is FALSE, so set it as FALSE temporarily
  
  out$re.detfn = FALSE
  
  
  #####################################################################################################
  #the 30th component: 'fit.freqs'
  #this part is not included in the TMB template yet
  
  out$fit.freqs = !is.null(cue.rates)
  
  ######################################################################################################
  #the 31st component: 'first.calls'
  #based on the instruction, set it to FALSE
  
  out$first.calls = FALSE
  
  
  #######################################################################################################
  #the 32nd component: 'scale.covs'
  #from the code, it seems this is a function receive data.frame as input and output a data.frame as well
  #not 100% sure
  
  out$scale.covs = scale.covs
  
  #######################################################################################################
  #the 33rd component: 'model.formula'
  if('D' %in% name.extend.par){
    out$model.formula = stats::as.formula(paste(c('gam.resp', as.character(par.extend$model$D)), collapse = ''))
  } else {
    out$model.formula = stats::as.formula('gam.resp~1')
  }
  
  ########################################################################################################
  #the 34th component: 'fgam'
  out$fgam = fgam
  
  ########################################################################################################
  #the 35th component: 'all.covariates'
  if(is(out$args$mask, 'list')){
    out$all.covariates = do.call('rbind', out$args$mask)
  } else {
    out$all.covariates = out$args$mask
  }
  
  if(is.scale){
    for(i in 1:ncol(out$all.covariates)){
      tem = out$all.covariates[, i]
      out$all.covariates[, i] = (tem - mean(tem)) / sd(tem)
    }
  }
  
  if(ncol(fgam$X) > 1){
    out$all.covariates = cbind(out$all.covariates, fgam$X[,-1])
  }
  out$all.covariates = as.data.frame(out$all.covariates)
  
  
  ############################################################################################################
  return(out)
}
