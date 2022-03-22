#' @export
sim.capt = function(fit, detfn, param, par_extend_model = NULL, traps, control_create_mask = list(), session_cov = NULL, trap_cov = NULL, loc_cov = NULL,
                    control_convert_loc2mask = list(), survey.length = NULL, ss.opts = NULL, cue.rates = NULL, n.sessions = NULL, n.rand = 1,
                    random.location = FALSE, sound.speed = 331){
  if(!missing(fit)){
    #if 'fit' is provided, get all information from the fitted object
    
    detfn = get_detfn(fit)
    
    #get the data frame which contains link function for each parameter
    dat_par = get_data_param(fit)
    
    #then get the parameter's values before back transforming
    param = get_coef(fit)
    
    #the "param" in simulation requires back transformed value
    for(i in names(param)){
      #if the length of one parameter is greater than 1, means it is extended
      #in this case, keep the values unchanged. we only do back transformation
      #when it is not extended
      if(length(param[[i]]) == 1){
        par_link = dat_par[which(dat_par$par == i), 'link']
        param[[i]] = link.fun(par_link, param[[i]])
      }
    }
    
    n.sessions = get_dims_tmb(fit)$n.sessions
    par.extend = get_par_extend(fit)
    traps = get_trap(fit)
    mask = get_mask(fit)
    traps = df_to_list(traps, n.sessions)
    mask = df_to_list(mask, n.sessions)
    survey.length = get_survey_length(fit)
    cue.rates = get_cue_rates(fit)
    sound.speed = get_sound_speed(fit)
    ss.opts = get_ss.opts(fit)
    
  } else {
    stopifnot(all(!missing(detfn), !missing(param), !missing(traps), !is.null(control_create_mask$buffer)))
    #extract n.sessions first and convert traps into list if it is not
    if(is(traps, 'list')){
      if(!is.null(n.sessions)) warning("'n.sessions' will be ignored as 'traps' is a list.")
      n.sessions = length(traps)

    } else {
      stopifnot(any(is(traps, 'matrix'), is(traps, 'data.frame')))
      if(is.null(n.sessions)){
        n.sessions = 1
      }
    }
    
    traps = df_to_list(traps, n.sessions)
    
    
    #create mask, traps is guaranteed to be a list, so mask must be a list
    control_create_mask$traps = traps
    mask = do.call('create.mask', control_create_mask)
    
    #if par_extend_model is assigned, we need to construct "par.extend" for simulation
    par.extend = par_extend_create(par_extend_model, loc_cov = loc_cov, mask = mask,
                                   control_convert_loc2mask = control_convert_loc2mask,
                                   session_cov = session_cov, trap_cov = trap_cov)
  }

  
  tem = sim.data.prepare(detfn, param, par.extend, traps, mask,
                         survey.length, random.location, n.sessions)
  dat_pars = tem$dat_pars
  dat.density = tem$dat.density
  dims = tem$dims
  info.bucket = tem$info.bucket
  
  if(n.rand == 1){
    output = sim.from.param(detfn, dat_pars, dat.density, random.location,
                            dims, info.bucket, ss.opts, cue.rates, sound.speed)
  } else {
    output = vector('list', n.rand)
    pb = utils::txtProgressBar(1, n.rand, style = 3)
    for(i in 1:n.rand){
      output[[i]] = sim.from.param(detfn, dat_pars, dat.density, random.location,
                                   dims, info.bucket, ss.opts, cue.rates, sound.speed)
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  return(output)
  
}

####################################################################################################################################
#simulation from assigned parameters

sim.data.prepare = function(detfn, param, par.extend, traps, mask, survey.length, random.location, n.sessions){
  
  if(is.null(survey.length)){
    survey.length = rep(1, n.sessions)
  } else {
    if(length(survey.length) == 1){
      survey.length = rep(survey.length, n.sessions)
    } else {
      stopifnot(length(survey.length) == n.sessions)
    }
  }
  
  stopifnot(detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth', 'ss'))
  

  #fill with basic information
  n.traps = sapply(traps, nrow)
  n.masks = sapply(mask, nrow)
  buffer = numeric(n.sessions)
  A = numeric(n.sessions)
  for(s in 1:n.sessions){
    buffer[s] = attr(mask[[s]], 'buffer')
    A[s] = attr(mask[[s]], 'area')
  }
  
  identical_traps = all(n.traps == n.traps[1])
  identical_masks = all(n.masks == n.masks[1])

  
  param.og = detfn.params(detfn)
  
  
  animal.model = 'mu' %in% names(param)
  is.toa = 'sigma.toa' %in% names(param)
  is.bearing = 'kappa' %in% names(param)
  is.dist = 'alpha' %in% names(param)
  is.ss = detfn == 'ss'
  
  param.og = c(param.og, 'D')
  if(animal.model) param.og = c(param.og, 'mu')
  if(is.toa) param.og = c(param.og, 'sigma.toa')
  if(is.bearing) param.og = c(param.og, 'kappa')
  if(is.dist) param.og = c(param.og, 'alpha')
  
  #create a session-trap-mask data frame for everything
  dat_list = vector('list', n.sessions)
  for(s in 1:n.sessions){
    
    traps[[s]] = as.matrix(traps[[s]])
    mask[[s]] = as.matrix(mask[[s]])
    
    #if not choose a location randomly inside a mask, then we can
    #just generate 'dx' and 'theta' here, otherwise, it is not necessary
    #to generate them here.
    if(!random.location){
      tem_dx = distances(traps[[s]], mask[[s]])
      tem_theta = bearings(traps[[s]], mask[[s]])
      dat_list[[s]] = data.frame(session = rep(s, n.masks[s] * n.traps[s]),
                                 mask = rep(1:n.masks[s], each = n.traps[s]),
                                 trap = rep(1:n.traps[s], n.masks[s]),
                                 dx = as.vector(tem_dx),
                                 theta = as.vector(tem_theta),
                                 mask_size = A[s])
    } else {
      dat_list[[s]] = data.frame(session = rep(s, n.masks[s] * n.traps[s]),
                                 mask = rep(1:n.masks[s], each = n.traps[s]),
                                 trap = rep(1:n.traps[s], n.masks[s]),
                                 trap_x = rep(traps[[s]][, 1], n.masks[s]),
                                 trap_y = rep(traps[[s]][, 2], n.masks[s]),
                                 x = rep(mask[[s]][,1], each = n.traps[s]),
                                 y = rep(mask[[s]][,2], each = n.traps[s]),
                                 mask_size = A[s])
    }
  }
  dat = do.call('rbind', dat_list)
  dat = sort.data(dat, "data.dists.thetas")
  
  if(is.null(par.extend)){
    name.extend.par = NULL
    #if par.extend is not provided, the param must be a list with scalars of "original values" of all parameters
    for(i in param.og){
      par.value = param[[i]]
      if(is.null(par.value)) stop(paste0('please provide parameter: ', i))
      dat[[i]] = par.value
    }
  } else {
    #if par.extend is provided, the param must be a list. if a parameter is extended, then the corresponding element
    #in "param" should be a vector with all coefficients of its linear predictor; otherwise, the element should be
    #simply scalars of the original value.
    if(!all(names(par.extend) %in% c('data', 'model', 'link', 'scale'))) {
      stop("'par.extend' only accepts 'data', 'model', 'link' and 'scale' as input.")
    }
    
    #since at the beginning of simulation, we do not have observation at all, we cannot make any scaling work.
    #we recommend users to standardize the numeric covariate first 
    if(!is.null(par.extend$scale)){
      warning(paste0("'scale' will be ignored and it will be regarded as FALSE."))
    }
    
    if(!all(c('data', 'model') %in% names(par.extend))){
      stop("'data' and 'model' must be provided.")
    }
    
    name.extend.par = names(par.extend$model)
    if('sigma.b0.ss' %in% name.extend.par) stop('sigma.b0.ss is not supported for parameter extension.')
    if(!all(name.extend.par %in% param.og)){
      stop('one or more parameters assigned in "par.extend" is not valid according to assigned detection function.')
    }

    
    tem_dat = dat[,c('session', 'trap', 'mask')]
    input_data = par.extend$data
    if(!all(names(input_data) %in% c('session', 'trap', 
                                     #'animal_ID', 
                                     'mask'))){
      stop("only 'session', 'trap', or 'mask' level data could be used as input.")
    }
    if(!is.null(input_data$session)){
      stopifnot('session' %in% colnames(input_data$session))
      tem_dat = merge(tem_dat, input_data$session, by = 'session', all.x = T)
    } 
    if(!is.null(input_data$trap)){
      input_data$trap = extend_dat_check(input_data$trap, 'trap', tem_dat, n.sessions, n.traps, identical_traps)
      tem_dat = merge(tem_dat, input_data$trap, by = c('session', 'trap'), all.x = T)
    }
    
    if(!is.null(input_data$mask)){
      input_data$mask = extend_dat_check(input_data$mask, 'mask', tem_dat, n.sessions, n.masks, identical_masks)
      tem_dat = merge(tem_dat, input_data$mask, by = c('session', 'mask'), all.x = T)
    }
    
    #include 'x' and 'y' in case they are used
    data.mask = as.data.frame(do.call('rbind', mask))
    colnames(data.mask) = c('x', 'y')
    data.mask$session = rep(1:n.sessions, n.masks)
    
    tem_mask_index = seq(n.masks[1])
    if(n.sessions > 1){
      for(s in 2:n.sessions) tem_mask_index = c(tem_mask_index, seq(n.masks[s]))
    }
    data.mask$mask = tem_mask_index
    
    tem_dat = merge(tem_dat, data.mask, by = c('session', 'mask'), all.x = T)
    
    #sort it to make sure it has the same order with 'dat'
    tem_dat = sort.data(tem_dat, "data.dists.thetas")
    tem_dat$gam.resp = 1
    for(i in param.og){
      par.value = param[[i]]
      if(is.null(par.value)) stop(paste0('please provide parameter: ', i))
      
      if(i %in% name.extend.par){
        foo = stats::as.formula(paste(c('gam.resp', as.character(par.extend$model[[i]])), collapse = ""))
        tem_model = mgcv::gam(foo, data = tem_dat, fit = FALSE)
        DX = tem_model$X
        cat(paste0("For extended parameter ", i, ", the corresponding covariates are: ", paste(colnames(DX), collapse = ", "), ".\n", 
                   "Please make sure the order of coefficients values in the 'param' for this parameter is correct.\n"))
        if(is.null(par.extend$link[[i]])){
          link = default.link(i)
        } else {
          link = par.extend$link[[i]]
        }
        
        dat[[i]] = unlink.fun(link, DX %*% par.value)
        
      } else {
        dat[[i]] = par.value
      }
    }
    
  }
  
  #D is animal density for no matter what kind of model is
  dat$D = dat$D * dat$mask_size
  dims = list(n.sessions = n.sessions, n.traps = n.traps, n.masks = n.masks, A = A)
  #create a data frame with session-mask and 'D' only, we will use this to generate
  #simulated number of animals
  dat.density = dat[, c('session', 'mask', 'x'[random.location], 'y'[random.location], 'D', 'mu'[animal.model])]
  dat.density = dat.density[!duplicated(dat.density[,c('session', 'mask')]),]
  dat.density = sort.data(dat.density, 'data.mask')
  dat.density$survey.length = rep(survey.length, n.masks)
  #remove useless columns from dat (because they are transferred to dat.density, or just useless)
  dat = dat[, -which(colnames(dat) %in% c('x', 'y', 'D', 'mu', 'mask_size'))]
  
  return(list(dat_pars = dat, dat.density = dat.density, dims = dims, 
              info.bucket = list(animal.model = animal.model, is.ss = is.ss,
                                is.dist = is.dist, is.bearing = is.bearing,
                                is.toa = is.toa, name.extend.par = name.extend.par)))

}


sim.from.param = function(detfn, dat_pars, dat.density, random.location, dims, info.bucket,
                          ss.opts, cue.rates, sound.speed){
  #firstly generate the number of animals for each mask
  dat.density$n_animals = rpois(nrow(dat.density), dat.density$D)
  dat.density = subset(dat.density, n_animals != 0)
  
  
  if(nrow(dat.density) == 0){
    #if no simulated animal, return an empty data frame
    return(data.frame(session = numeric(0), trap = numeric(0), ID = numeric(0)))
  }
  
  #split multiple animals if they are both in the same mask, so we have one row for each individual
  dat.density = split_item(dat.density, 'n_animals')
  
  if(random.location){
      tem_list = vector('list', dims$n.sessions)
    for(s in 1:dims$n.sessions){
      tem_density = subset(dat.density, session == s)
      
      #calculate the side length of each mask, since the mask points are evenly spread,
      #and each mask is a square, we only need to calculate average interval of x or y
      #which should be exactly equal to the side length
      tem = unique(tem_density$x)
      side_len = (max(tem) - min(tem))/(length(tem) - 1)
      
      #randomly select an coordination from the mask for each ID
      tem_density$x = tem_density$x + runif(nrow(tem_density), -0.5 * side_len, 0.5 * side_len)
      tem_density$y = tem_density$y + runif(nrow(tem_density), -0.5 * side_len, 0.5 * side_len)
      
      tem_list[[s]] = tem_density
    }
    dat.density = do.call('rbind', tem_list)
  }
  
  #generate random n_calls
  if(info.bucket$animal.model){
    dat.density$animal_ID = 1:nrow(dat.density)
    #since "n_animals" has been split, each row only denotes one animal, so here we could
    #literally use mu*survey.length as lambda for the Poisson distribution
    dat.density$n_calls = rpois(nrow(dat.density), dat.density$mu * dat.density$survey.length)
  } else {
    #when cue.rates is not assigned, we regard this as a "cue density model", which means n_calls is equivalent to n_animals * survey.length
    if(is.null(cue.rates)){
      #since we have split n_animals, now each row represents one animal, so n_animals === 1
      dat.density$n_calls = dat.density$survey.length
    } else {
      #if cue.rates is assigned, we are doing a "animal density model", so we generate random n_calls based on cue.rates and survey.length
      cue.rates = mean(cue.rates)
      dat.density$n_calls = rpois(nrow(dat.density), cue.rates * dat.density$survey.length)
    }
  }
  
  #split n_calls and assign ID
  dat.density = split_item(dat.density, 'n_calls')
  dat.density$ID = 1:nrow(dat.density)
  
  #drop useless columns
  dat.density = dat.density[,-which(colnames(dat.density) %in% c('n_animals', 'n_calls', 'survey.length'))]
  
  #merge the dat.density (simulated animals and calls) and the dat_pars (values of each parameter) to generate our simulated capture history
  dat_capt = merge(dat.density, dat_pars, by = c('session', 'mask'), all.x = T)
  
  if(random.location){
    #re-calculate distance and theta
    dat_capt$dx = sqrt((dat_capt$trap_x - dat_capt$x) ^ 2 +
                            (dat_capt$trap_y - dat_capt$y) ^ 2)
    
    dat_capt$theta = bearings_by_vec(dat_capt$trap_x, dat_capt$trap_y,
                                     dat_capt$x, dat_capt$y)
    #drop useless columns
    dat_capt = dat_capt[,-which(colnames(dat_capt) %in% c('x', 'y', 'trap_x', 'trap_y'))]
  }
  
  #calculate the probability of each trap detects each call
  name.param.det = detfn.params(detfn)
  param.det = vector('list', length(name.param.det))
  names(param.det) = name.param.det
  for(i in name.param.det) param.det[[i]] = dat_capt[[i]]
  
  if(info.bucket$is.ss){
    ss.link = ss.opts[['ss.link']]
    if(is.null(ss.link)) ss.link = 'identity'
    cut_off = ss.opts[['cutoff']]
    stopifnot(!is.null(cut_off))
    #ss is a little different, if the simulated ss > cut off, the detection probability is 1
    #and 0 otherwise
    essx = det_prob(detfn, param.det, dat_capt$dx, ss.link)
    dat_capt[['ss']] = rnorm(nrow(dat_capt), essx, param.det[["sigma.ss"]])
    dat_capt[['det_prob']] = ifelse(dat_capt[['ss']] > cut_off, 1, 0)
  } else {
    dat_capt[['det_prob']] = det_prob(detfn, param.det, dat_capt$dx)
  }
  
  dat_capt$bincapt = stats::rbinom(nrow(dat_capt), 1, dat_capt[['det_prob']])
  dat_capt = subset(dat_capt, bincapt == 1)
  
  if(nrow(dat_capt) == 0){
    #if no simulated detection, return an empty data frame
    return(data.frame(session = numeric(0), trap = numeric(0), ID = numeric(0)))
  }
  
  
  output = data.frame(session = dat_capt$session, ID = dat_capt$ID, occasion = 1, trap = dat_capt$trap)
  if(info.bucket$animal.model) output$animal_ID = dat_capt$animal_ID
  if(info.bucket$is.ss) output$ss = dat_capt$ss
  #simulate extra info
  if(info.bucket$is.bearing){
    output[['bearing']] = 0

    if('kappa' %in% info.bucket$name.extend.par){
      #if kappa is extended, we need to calculate bearing one by one since rvm() does not
      #support vectorize computing with different parameters
      for(i in 1:nrow(dat_capt)){
        output[i, 'bearing'] = CircStats::rvm(1, dat_capt[i, 'theta'], dat_capt[i, 'kappa'])
      }
    } else {
      #if kappa is not extended, we could just take the first value since they are all the same
      #and use vectorize computing
      k = dat_capt[['kappa']][1]
      output[['bearing']] = (dat_capt[['theta']] + CircStats::rvm(nrow(dat_capt), 0, k)) %% (2*pi)
    }
  }
  
  if(info.bucket$is.dist){
    shape = dat_capt[['alpha']]
    scale = dat_capt[['dx']] / shape
    output[['dist']] = rgamma(nrow(dat_capt), shape = shape, scale = scale)
  }
  
  if(info.bucket$is.toa){
    output[['toa']] = rnorm(nrow(dat_capt), dat_capt[['dx']]/sound.speed, dat_capt[['sigma.toa']])
  }
  
  if(info.bucket$animal.model){
    output = output[order(output$session, output$animal_ID, output$ID, output$trap), ]
  } else {
    output = output[order(output$session, output$ID, output$trap), ]
  }
  
  return(output)
  
}


