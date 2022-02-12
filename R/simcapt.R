#' @export
sim.capt = function(fit, detfn, param, par.extend = NULL, traps, mask, survey.length = NULL, ss.opts = NULL, cue.rates = NULL,
                    n.sessions = NULL, n.rand = 1, random.location = FALSE, sound.speed = 331){
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
      #in this case, keep the values the same. we only do back transformation
      #when it is not extended
      if(length(param[[i]]) == 1){
        par_link = dat_par[which(dat_par$par == i), 'link']
        param[[i]] = link.fun(par_link, param[[i]])
      }
    }
    
    par.extend = get_par_extend(fit)
    traps = get_trap(fit)
    mask = get_mask(fit)
    survey.length = get_survey_length(fit)
    cue.rates = get_cue_rates(fit)
    sound.speed = get_sound_speed(fit)
    ss.opts = get_ss.opts(fit)
  } else {
    stopifnot(all(!missing(detfn), !missing(param), !missing(traps), !missing(mask)))
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
  stopifnot(detfn %in% c('hn', 'hhn', 'hr', 'th', 'lth', 'ss'))
  
  #confirm n.sessions firstly
  if(is(traps, 'list')){
    if(!is.null(n.sessions)) warning("'n.sessions' will be ignored as 'traps' is a list.")
    n.sessions = length(traps)
    if(is(mask, 'list')){
      stopifnot(length(mask) == n.sessions)
    } else {
      stopifnot(any(is(mask, 'matrix'), is(mask, 'data.frame')))
      tem = mask
      mask = vector('list', n.sessions)
      for(s in 1:n.sessions) mask[[s]] = tem
    }
  } else {
    if(is.null(n.sessions)){
      if(is(mask, 'list')){
        n.sessions = length(mask)
      } else {
        stopifnot(any(is(mask, 'matrix'), is(mask, 'data.frame')))
        mask = list(mask)
        n.sessions = 1
      }
    } else {
      if(is(mask, 'list')){
        stopifnot(length(mask) == n.sessions)
      } else {
        stopifnot(any(is(mask, 'matrix'), is(mask, 'data.frame')))
        tem_mask = mask
        mask = vector('list', n.sessions)
        for(s in 1:n.sessions){
          mask[[s]] = tem_mask
        }
      }
    }
    
    if(is.null(survey.length)){
      survey.length = rep(1, n.sessions)
    } else {
      if(length(survey.length) == 1){
        survey.length = rep(survey.length, n.sessions)
      } else {
        stopifnot(length(survey.length) == n.sessions)
      }
    }
    
    stopifnot(any(is(traps, 'matrix'), is(traps, 'data.frame')))
    tem_traps = traps
    traps = vector('list', n.sessions)
    for(s in 1:n.sessions){
      traps[[s]] = tem_traps
    }
    
  }

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
    #if par.extend is not provided, the param must be a list with "original values" of all parameters
    for(i in param.og){
      par.value = param[[i]]
      if(is.null(par.value)) stop(paste0('please provide parameter: ', i))
      dat[[i]] = par.value
    }
  } else {
    #if par.extend is provided, the param must be a list. if a parameter is extended, then the corresponding element
    #in "param" should be a vector with all coefficients of its linear predictor; otherwise, the element should be
    #simply the original value.
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
      input_data$trap = extend_dat_check(input_data$trap, 'trap', tem_dat, n.sessions, identical_traps)
      tem_dat = merge(tem_dat, input_data$trap, by = c('session', 'trap'), all.x = T)
    } 
    if(!is.null(input_data$mask)){
      input_data$mask = extend_dat_check(input_data$mask, 'mask', tem_dat, n.sessions, identical_masks)
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

######################################################################################################################################
#simulation from model output, animal_ID embedded model is not included yet.

sim.from.fit = function(fit){

  #get pop density for each mask
  data_density = get_D_tmb(fit)
  #area for each mask
  area_unit = get_area_unit_tmb(fit)
  
  
  dims = get_dims_tmb(fit)
  
  #calculate lambda of each mask for poisson distribution
  data_density$lambda_unit = data_density$D_est * rep(area_unit, dims$n.masks)
  
  
  data_density$n_animals = 0
  data_density$n_calls = 0
  
  for(s in 1:dims$n.sessions){
    index_data_density_session = which(data_density$session == s)
    #get the lambda for this session, which sums up the lambda for each mask in this session
    lambda_unit_session = data_density$lambda_unit[index_data_density_session]
    lambda_session = sum(lambda_unit_session)
    #sample a random number from Poisson distribution with lambda_session
    pop_session = rpois(1, lambda_session)
    
    #allocate this number into each mask by a multinomial distribution
    #with (lambdas of each mask / lambda_session) as their probabilities
    
    #actually the operation above with Poisson(lambda_session) first and
    #multinomial later is equivalent to sample a number for each mask with
    #Poisson(lambda_each_mask)
    
    tem = as.vector(rmultinom(1, pop_session, prob = lambda_unit_session))
    data_density[index_data_density_session, 'n_animals'] = tem
    
  }
  #simulate capture history based on the simulated density and detection function
  data_capt = sim_det_history(fit, data_density)
  
  fit_types = get_fit_type(fit)
  
  if(sum(data_capt[['n_det']]) == 0){
    #if there is no simulated detection, return a empty data frame
    col_names = c('session', 'ID', 'occasion', 'trap',
                  'toa'[fit_types['toa']],
                  'ss'[fit_types['ss']],
                  'dist'[fit_types['dist']],
                  'bearing'[fit_types['bearing']])
    output = vector('list', length(col_names))
    names(output) = col_names
    for(i in col_names) output[[i]] = numeric(0)
    output = as.data.frame(output)
    
  } else {
    #simulate extra info (bearing, dist, ss, toa)
    data_capt = sim_extra_info(fit, data_capt)
    
    #take out the observations with no detection
    o = subset(data_capt, data_capt$n_det > 0)
    
    #pack all result into a data frame
    output = data.frame(session = o$session, ID = o$ID, occasion = 1, trap = o$trap)
    
    if(fit_types['toa']){
      output[['toa']] = o$toa
    }
    
    if(fit_types['ss']){
      output[['ss']] = o$ss
    }
    
    if(fit_types['dist']){
      output[['dist']] = o$dist
    }
    
    if(fit_types['bearing']){
      output[['bearing']] = o$bearing
    }
    
    output = output[order(output$session, output$ID, output$trap),]
  }
  
  return(output)
}


###########################################################################################################
#simulate detection history, the key component of simulation
sim_det_history = function(fit, data_density){
  #source('support_functions.r', local = TRUE)
  #source('detfn_tmb.r', local = TRUE)
  data_full = get_data_full(fit)
  data_mask = get_data_mask(fit)
  data_trap = get_data_trap(fit)
  #get design matrices for non-mask level and mask level
  #indicated as "full" and "mask" respectively
  DX_full = get_DX_full(fit)
  DX_mask = get_DX_mask(fit)
  #get the small data frame which records the link function
  #for each parameter and the numbers of covariates for each
  #parameter in non-mask level (in "DX_full") and mask level
  #(in "DX_mask")
  param_table = get_data_param(fit)
  #get all estimated values for each parameter's covariates
  param = get_coef(fit)
  det_fn = get_detfn(fit)
  #get original parameters' names for detection function
  #and extra information (will be used in another function)
  #since "D" is not used in detection function nor extra info
  #remove it
  param_name = get_param_og(fit)
  param_name = param_name[-which(param_name == 'D')]
  
  dims = get_dims_tmb(fit)
  
  #survey length for each session
  survey_length = get_survey_length(fit)
  #average cue rate
  avg_cue_rates = get_cue_rates(fit)
  
  #because in the data.full, if one session has no detection
  #it will still record n.trap's rows, a special n.id is necessary
  n_ID_data_full = ifelse(dims$n.IDs == 0, 1, dims$n.IDs)
  
  #create a new data.full but this time the columns are
  #back-transformed values for each parameters
  data_full_param_value = data_full[, c('session', 'ID', 'trap')]
  data_mask_param_value = data_mask[, c('session', 'mask')]
  
  #get the data frame records distance and angles between each mask and trap
  data_dist_theta = get_dist_theta(fit)
  
  for(i in param_name){
    tem = subset(param_table, param_table$par == i)
    link = tem$link
    #n_col essentially means number of covariates
    n_col_full = tem$n_col_full
    n_col_mask = tem$n_col_mask
    
    #get the design matrices
    tem_DX_full = DX_full[[i]]
    tem_DX_mask = DX_mask[[i]]
    
    #get the values for covariates and make it into a matrix with dim of k*1
    tem_coef_full = matrix(param[[i]][1:n_col_full], ncol = 1)
    
    #design matrix multiple parameter vector for non-mask level to get the
    #non-mask level part of parameters values
    data_full_param_value[[i]] = tem_DX_full %*% tem_coef_full
    
    #since 'ID' is not extendable, shrink this data set
    tem = data_full_param_value[!duplicated(data_full_param_value[, c('session', 'trap')]),
                                c('session', 'trap', i)]
    #merge the data with distance and parameters values into one data set
    data_dist_theta = merge(data_dist_theta, tem, by = c('session', 'trap'))
    
    #if there is any mask level covariates, do the same operation with non-mask level
    #covariates to obtain the parameters values in mask level part
    if(n_col_mask > 0){
      tem_coef_mask = matrix(param[[i]][(n_col_full + 1) : (n_col_full + n_col_mask)],
                             ncol = 1)
      data_mask_param_value[[paste0(i, '_mask')]] = tem_DX_mask %*% tem_coef_mask
      
      tem = data_mask_param_value[, c('session', 'mask', paste0(i, '_mask'))]
      tem = merge(data_dist_theta, data_mask_param_value, by = c('session', 'mask'))
      #combine the two parts to get the final parameters values
      data_dist_theta[[i]] = tem[[i]] + tem[[paste0(i, '_mask')]]
    }
    
    #based on link function, back transform the values of parameters
    data_dist_theta[[i]] = unlink.fun(link = link, value = data_dist_theta[[i]])
  }
  
  #merge the distances, parameters values and simulated number of animals together
  data_dist_theta = merge(data_dist_theta, data_density[,c('session', 'mask', 'n_animals')],
                          by = c('session', 'mask'))
  
  #shrink the data set by dropping the masks with no simulated animal
  data_capt = subset(data_dist_theta, data_dist_theta$n_animals != 0)
  
  #merge the coordination of each trap and each mask to this data frame
  data_capt = merge(data_capt, data_trap, by = c('session', 'trap'), all.x = TRUE)
  
  data_capt = merge(data_capt, data_mask, by = c('session', 'mask'), all.x = TRUE)
  
  #sort it firstly, in order to assign ID to each call after the "split" operation much easier
  data_capt = data_capt[order(data_capt$session, data_capt$mask, data_capt$trap),]
  
  
  #separate n_animals, e.g: if in one mask, the n_animals == 2, then duplicate this row 2 times
  data_capt = split_item(data_capt, "n_animals")
  
  tem_list = vector('list', dims$n.sessions)
  
  #randomly allocate each animal inside the area of that mask point, and generate n_calls
  for(s in 1:dims$n.sessions){
    tem_capt = subset(data_capt, session == s)
    tem_mask = subset(data_mask, session == s)
    
    #calculate the width of each mask
    h_x = unique(diff(sort(unique(tem_mask$x))))
    if(length(h_x) != 1) stop('masks are not evenly splited.')
    #calculate the height of each mask
    h_y = unique(diff(sort(unique(tem_mask$y))))
    if(length(h_y) != 1) stop('masks are not evenly splited.')
    
    #randomly select an coordination from the mask for each ID
    tem_capt$x = tem_capt$x + runif(nrow(tem_capt), -0.5 * h_x, 0.5 * h_x)
    tem_capt$y = tem_capt$y + runif(nrow(tem_capt), -0.5 * h_y, 0.5 * h_y)
    
    #re-calculate distance and theta
    tem_capt$dx = sqrt((tem_capt$trap_x - tem_capt$x) ^ 2 +
                         (tem_capt$trap_y - tem_capt$y) ^ 2)
    
    tem_capt$theta = bearings_by_vec(tem_capt$trap_x, tem_capt$trap_y,
                                     tem_capt$x, tem_capt$y)
    
    #about the theta part, there is one point needs attention. For e.g.,
    #when the original theta is 1.469693, then the new theta after taking
    #a random location within a mask is 4.612040, although they seems
    #to be different, but if we take tan(theta), we could see they are
    #quite similar, one is 9.857147 and the other one is 9.931751
    
    ##############################################################################
    #then simulate n_calls
    
    #temporarily, we simplify the simulation of calls of one animal by using just
    #a sample from a poisson distribution with only "average cue rate", this part
    #may be more complicated in the future
    
    n_animals_vec = tem_capt[['n_animals']]
    
    if(fit$fit.freqs){
      cue_rates = rpois(length(n_animals_vec), avg_cue_rates)
    } else {
      #if the original cue.rate input in the model fitting is "NULL" then it means
      #we don't have the concept of "animal" and "calls" in the first place
      #so we do not proceed with sampling number of calls, just fix it to be 1
      cue_rates = 1
    }
    
    #get number of calls for each mask
    tem_capt[['n_calls']] = n_animals_vec * cue_rates * survey_length[s]
    
    tem_list[[s]] = tem_capt
  }
  
  data_capt = do.call('rbind', tem_list)
  
  #then split n_calls just the same as splitting n_animals
  data_capt = split_item(data_capt, "n_calls")
  
  
  
  #since data is well sorted, we could easily assign ID to each call
  data_capt[['ID']] = 0
  for(s in 1:dims$n.sessions){
    index_s = which(data_capt[['session']] == s)
    data_capt[['ID']][index_s] = rep(seq(length(index_s) / dims$n.traps[s]), each = dims$n.traps[s])
  }
  
  dx = data_capt[['dx']]
  
  #create variables for each parameter will be used in det_fn part
  det_par = vector('list', length(param_name))
  names(det_par) = param_name
  for(i in param_name) det_par[[i]] = data_capt[[i]]
  #calculate detection probability
  if(det_fn != 'ss'){
    data_capt[['det_prob']] = det_prob(det_fn, det_par, dx)
  } else {
    ss.link = get_ss_link(fit)
    cut_off = get_cutoff(fit)
    #ss is a little different, if the simulated ss > cut off, the detection probability is 1
    #and 0 otherwise
    essx = det_prob(det_fn, det_par, dx, ss.link)
    data_capt[['ss']] = rnorm(nrow(data_capt), essx, det_par[["sigma.ss"]])
    data_capt[['det_prob']] = ifelse(data_capt[['ss']] > cut_off, 1, 0)
  }
  
  #remove used variables to release RAM
  rm(det_par)
  
  data_capt[['n_det']] = 0
  data_capt[['n_det']] = rbinom(nrow(data_capt), data_capt[['n_calls']], data_capt[['det_prob']])
  return(data_capt)
}





sim_extra_info = function(fit, data_capt){
  fit_types = get_fit_type(fit)
  #check which rows have detection, only simulate for those rows
  index_det = which(data_capt[['n_det']] > 0)
  n_rows = length(index_det)
  if(fit_types['toa']){
    sound_speed = get_sound_speed(fit)
    data_capt[['toa']] = 0
    data_capt[['toa']][index_det] = rnorm(n_rows, data_capt[['dx']][index_det] / sound_speed,
                                          data_capt[['sigma.toa']][index_det])
  }
  
  
  if(fit_types['dist']){
    data_capt[['dist']] = 0
    shape = data_capt[['alpha']][index_det]
    scale = data_capt[['dx']][index_det] / shape
    data_capt[['dist']][index_det] = rgamma(n_rows, shape = shape, scale = scale)
  }
  
  if(fit_types['bearing']){
    
    data_capt[['bearing']] = 0
    #function 'rvm' does not support vector computation
    if(length(get_coef(fit)$kappa) > 1){
      #if kappa is extended, do a loop to each row with a detection
      for(i in index_det){
        data_capt[['bearing']][i] = CircStats::rvm(1, data_capt[['theta']][i], data_capt[['kappa']][i])
      }
    } else {
      #if kappa is not extended, just take the first one as its value
      k = data_capt[['kappa']][1]
      data_capt[['bearing']][index_det] = (data_capt[['theta']][index_det] +
                                             CircStats::rvm(n_rows, 0, k)) %% (2*pi)
    }
    
  }
  
  return(data_capt)
}
