boot.ascr = function(fit, N, n.cores = 1, M = 10000, infotypes = NULL, seed = NULL){
  arguments = fit$args
  
  #deal with infotypes. Here we directly change the fit object because this object will
  #be used for large scale simulation later, if one info such as 'bearing' is assigned not
  #to be boot, then it should not be simulated during simulation to save resources
  
  if(!is.null(infotypes)){
    all_infotypes = names(fit$fit.types)
    if(any(!infotypes %in% all_infotypes)){
      stop(paste0('argument infotypes only accept following characters: ', 
                  paste(all_infotypes, collapse = ', ')))
    }
    
    for(i in all_infotypes){
      fit$fit.types[i] = (i %in% infotypes)
    }
    
    fit$infotypes = infotypes
  } else {
    #if NULL, then default setting is all info types used in the model
    infotypes = fit$infotypes
  }
  
  
  #since in parametric bootstrap we know the true values of each parameter
  #there are benefits to use the true values as start values
  sv.og = arguments$sv
  
  arguments$sv = get_sv_for_boot(fit)
  
  #extract some important components from original fit, maybe they need be adjusted
  #or they need be retained
  dims= get_dims_tmb(fit)
  traps = arguments$traps
  masks = arguments$mask
  
  ##mrds, I'm not sure how to deal with 'mrds' in bootstrap yet, because the number
  ##of detection from simulation is random, however 'mrds' has provided the number
  ##of detection with their location. how to solve this conflict?
  is.mrds = fit$fit.types['mrds']
  
  ###############################################################################
  #temporarily force it to be FALSE since not sure how to deal with it yet.
  is.mrds = FALSE
  if('mrds' %in% infotypes){
    infotypes = infotypes[-which(infotypes == 'mrds')]
    fit$fit.types['mrds'] = FALSE
    fit$infotypes = infotypes
    if(fit$n.sessions == 1){
      arguments$capt[['mrds']] = NULL
    } else {
      for(s in 1:fit$n.sessions) arguments$capt[[s]][['mrds']] = NULL
    }
    
  }

  ##############################################################################
  
  if(is.mrds){
    if(fit$n.sessions == 1){
      mrds.loc = arguments$capt[['mrds']]
    } else {
      mrds.loc = lapply(arguments$capt, function(x) x[['mrds']])
    }
  }

  
  
  ##cue.rate
  cue.rates = arguments$cue.rates
  
  #coefficients names
  tem = coef(fit, 'linked')
  n_pars = length(tem)
  par_names = names(tem)
  
  #set seed
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    set.seed(sample(1:1e8, size = 1))
  }
  
  seed_boot = sample(1:1e8, size = N)
  seed_mce = sample(1:1e8, size = 1)
  
  
  if(n.cores == 1){
    res = matrix(NA, nrow = N, ncol = n_pars + 1)
    colnames(res) = c(par_names, 'maxgrad')
    
    for(n in 1:N){
      set.seed(seed_boot[n])
      capture_sim = sim.capt(fit)
      if(nrow(capture_sim) > 0){
        arguments$capt = get_capt_for_boot(capture_sim, dims, infotypes)
      } else {
        
      }
    }
  }

  
  
}