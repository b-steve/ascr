#' @export
sim.capt = function(fit){

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
  
  #simulate extra info (bearing, dist, ss, toa)
  data_capt = sim_extra_info(fit, data_capt)
  
  #take out the observations with no detection
  o = subset(data_capt, data_capt$n_det > 0)
  
  #pack all result into a data frame
  output = data.frame(session = o$session, ID = o$ID, occasion = 1, trap = o$trap)
  fit_types = get_fit_type(fit)
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
    mu = det_prob(det_fn, det_par, dx, ss.link)
    data_capt[['ss']] = rnorm(nrow(data_capt), mu, det_par[["sigma.ss"]])
    data_capt[['det_prob']] = ifelse(data_capt[['ss']] > cut_off, 1, 0)
  }
  
  #remove used variables to release RAM
  rm(det_par)
  
  data_capt[['n_det']] = 0
  
  #at least we need some detection history, otherwise we cannot proceed with any further procedure
  while(sum(data_capt[['n_det']]) == 0){
    data_capt[['n_det']] = rbinom(nrow(data_capt), data_capt[['n_calls']], data_capt[['det_prob']])
  }
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



###################################################################################################
#bellows are only used for development, delete them when formally building the package

sim_result_agg = function(sim_result){
  tem = sim_result[[1]]
  n.sim = length(sim_result)
  
  output = vector('list', length(tem))
  names(output) = names(tem)
  
  for(i in names(output)){
    output[[i]] = vector('list', length(tem[[i]]))
    names(output[[i]]) = names(tem[[i]])
    for(j in names(output[[i]])) {
      output[[i]][[j]] = numeric(n.sim)
      for(n in 1:n.sim) output[[i]][[j]][n] = sim_result[[n]][[i]][j]
    }
  }
  
  return(output)
}

#input is the aggregated result
sim_result_plot = function(sim_agg, coef_real, pre_dir){
  for(i in names(coef_real)){
    tem_sim = sim_agg[[i]]
    tem_real = coef_real[[i]]
    for(j in names(tem_real)){
      sim_vec = tem_sim[[j]]
      real_val = tem_real[j]
      dir = paste0(pre_dir, "/", j, ".jpeg")
      jpeg(filename = dir, width = 1920, height = 1080)
      hist(sim_vec, main = paste0('simulation of ', j), xlab = j)
      abline(v = real_val, col = 2)
      dev.off()
    }
  }
  return(0)
}
