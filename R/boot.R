#' Title
#'
#' @param fit 
#' @param N 
#' @param n.cores 
#' @param M 
#' @param infotypes 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
boot.ascr = function(fit, N, n.cores = 1, M = 10000, infotypes = NULL, seed = NULL){
 
  fit.og = fit
  arguments = fit$args
  
  #deal with infotypes
  
  if(!is.null(infotypes)){
    #all possible infotypes
    all_infotypes = fit$infotypes
    
    if(is.null(all_infotypes)){
      warning('argument "infotypes" will be ignored as there is no extra infomation in the original model.')
      infotypes = NULL
    } else {
      #check it is a list or character vector
      stopifnot(any(is(infotypes, 'character'), is(infotypes, 'list')))
      if(is(infotypes, 'list')) stopifnot(all(sapply(infotypes, function(x) is(x, 'character'))))
      
      #if it is a character vector, convert it to a list to make future operation easier
      if(is(infotypes, 'character')) infotypes = list(infotypes)
      
      if(any(!do.call('c', infotypes) %in% all_infotypes)){
        stop(paste0('argument infotypes only accept following characters: ', 
                    paste(all_infotypes, collapse = ', ')))
      }
    }
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
  
  ##cue.rate, when cue.rate is NULL, length(cue.rates) return 0.
  cue.rates <- arguments$cue.rates
  is_sim_cue <- (length(cue.rates) > 1)
  
  
  ###############################################################################
  ##mrds, I'm not sure how to deal with 'mrds' in bootstrap yet, because the number
  ##of detection from simulation is random, however 'mrds' has provided the number
  ##of detection with their location. how to solve this conflict?
  #temporarily force it to be FALSE since not sure how to deal with it yet.
  if(fit$fit.types['mrds']){
    fit$infotypes = fit$infotypes[which(fit$infotypes!='mrds')]
    if(fit$n.sessions == 1){
      arguments$capt[['mrds']] = NULL
    } else {
      for(s in 1:fit$n.sessions) arguments$capt[[s]][['mrds']] = NULL
    }
  }
  
  if(!is.null(infotypes)){
    for(i in 1:length(infotypes)){
      if('mrds' %in% infotypes[[i]]) infotypes[[i]] = infotypes[[i]][which(infotypes[[i]] != 'mrds')]
    }
  }
  
  
  fit$fit.types['mrds'] = FALSE
  ##############################################################################
  is.mrds = fit$fit.types['mrds']
  if(is.mrds){
    if(fit$n.sessions == 1){
      mrds.loc = arguments$capt[['mrds']]
    } else {
      mrds.loc = lapply(arguments$capt, function(x) x[['mrds']])
    }
  }
  

  
  #coefficients names
  coefs = coef(fit, 'linked')
  n_pars = length(coefs)
  par_names = names(coefs)
  
  #set seed
  if(!is.null(seed)){
    set.seed(seed)
  } else {
    set.seed(sample(1:1e8, size = 1))
  }
  
  seed_boot = sample(1:1e8, size = N)
  seed_mce = sample(1:1e8, size = 1)
  
  
  if(n.cores == 1){
    ncol_res = n_pars + 1
    colnames_res = c(par_names, 'maxgrad')
    res = matrix(NA, nrow = N, ncol = ncol_res)
    colnames(res) = colnames_res
    

    ########################################################
    res = boot_step(seed = seed_boot,
                    N = N,
                    fit = fit,
                    arguments = arguments,
                    dims = dims,
                    infotypes = fit$infotypes,
                    is_sim_cue = is_sim_cue,
                    cue.rates = cue.rates,
                    len_output = ncol_res,
                    name_output = colnames_res)
    ########################################################
    
    
    #Additional bootstraps
    if(!is.null(infotypes)){
      extra.res <- vector(mode = "list", length = length(infotypes))
      names(extra.res) <- names(infotypes)
      
      for (i in seq(from = 1, by = 1, along.with = infotypes)){
        new_args <- arguments
        new_args$capt <- arguments$capt[c("bincapt", infotypes[[i]])]
        new_fit <- suppressWarnings(do.call("fit_og", new_args))
        tem = coef(new_fit, 'linked')
        
        new_ncol_res = length(tem) + 1
        new_colnames_res = c(names(tem), 'maxgrad')
        extra.res[[i]] <- matrix(NA, nrow = N, ncol = new_ncol_res)
        colnames(extra.res[[i]]) <- new_colnames_res
        
 
        ####################################################################################  
        extra.res[[i]] <- suppressWarnings(boot_step(seed = seed_boot,
                                                     N = N,
                                                     fit = new_fit,
                                                     arguments = new_args,
                                                     dims = dims,
                                                     infotypes = infotypes[[i]],
                                                     is_sim_cue = is_sim_cue,
                                                     cue.rates = cue.rates,
                                                     len_output = new_ncol_res,
                                                     name_output = new_colnames_res))
        ####################################################################################
        
        
      }
      
      #end of if(!is.null(infotypes))
    } else {
      extra.res = NULL
    }
    
    #end of n.core == 1
  }
  
  
  ## Calculating bootstrapped standard errors, correlations and
  ## covariances.
  maxgrads <- res[, ncol(res)]
  ## Removing maximum gradient component.
  res <- res[, -ncol(res), drop = FALSE]
  se <- apply(res, 2, sd, na.rm = TRUE)
  names(se) <- par_names
  colnames(res) <- par_names
  corr <- diag(n_pars)
  dimnames(corr) <- list(par_names, par_names)
  vcov <- diag(se^2)
  dimnames(vcov) <- list(par_names, par_names)
  for (i in 1:(n_pars - 1)){
    for (j in (i + 1):n_pars){
      corr[i, j] <- corr[j, i] <- cor(res[, i], res[, j], use = "na.or.complete")
      vcov[i, j] <- vcov[j, i] <- corr[i, j]*se[i]*se[j]
    }
  }
  bias <- apply(res, 2, mean, na.rm = TRUE) - coefs
  ## Bootstrap to calculate MCE for bias and standard errors.
  bias.mce <- se.mce <- numeric(n_pars)
  names(bias.mce) <- names(se.mce) <- par_names
  if (M > 0){
    set.seed(seed_mce)
    converged <- which(!is.na(res[, 1]))
    n.converged <- length(converged)
    mce.boot <- matrix(sample(converged, size = n.converged*M,
                              replace = TRUE), nrow = M,
                       ncol = n.converged)
    #browser()
    for (i in par_names){
      par.boot <- matrix(res[mce.boot, i], nrow = M, ncol = n.converged)
      bias.mce[i] <- sd(apply(par.boot, 1, mean) - coefs[i] - bias[i])
      se.mce[i] <- sd(apply(par.boot, 1, sd))
    }
  } else {
    bias.mce <- NA
    se.mce <- NA
  }
  out <- fit.og
  boot <- list(boots = res, se = se, se.mce = se.mce, cor = corr, vcov = vcov,
               bias = bias, bias.mce = bias.mce, maxgrads = maxgrads,
               extra.boots = extra.res)
  out$boot <- boot
  class(out) <- c("ascr.boot", class(fit))
  return(out)
}

#continue from here========================================================================
boot_step = function(seed, N, fit, arguments, dims, infotypes, is_sim_cue, 
                     cue.rates, len_output, name_output){
  set.seed(seed)
  output = matrix(NA, nrow = N, ncol = len_output)
  #main bootstrap
  capture_sim = sim.capt(fit = fit, n.rand = N)$capt
  for(n in 1:N){
    if(nrow(capture_sim[[n]]) > 0){
      arguments$capt = get_capt_for_boot(capture_sim[[n]], dims, infotypes)
      if(is_sim_cue){
        arguments$cue.rates = sample(cue.rates, replace = TRUE)
      }
      
      fit_boot = suppressWarnings(try(do.call('fit_og', arguments), silent = TRUE))
      #If unconverged, refit model with default start values.
      if ("try-error" %in% class(fit_boot) || fit_boot$maxgrad < -0.01){
        arguments$sv <- NULL
        fit_boot <- suppressWarnings(try(do.call("fit_og", arguments), silent = TRUE))
      }
      #If still unconverged, give up and report NA, which is the default value of the output,
      #so we need to do nothing. And if it converged, assign result to the output
      if (!("try-error" %in% class(fit_boot) || fit_boot$maxgrad < -0.01)){
        output[n,] = c(coef(fit_boot, 'linked'), fit_boot$maxgrad)
      } 
      
    } else {
      #default of output is NA, so no need to modify the columns of parameters excluding "D" related parameter
      #change the intercept or its equivalent to -Inf, and other "D" related coefficients to 0
      i_D_int = which(name_output %in% c('D_link', 'D.(Intercept)_link'))
      i_D_other = setdiff(which(grepl("^D", name_output)), i_D_int)
      output[n, i_D_int] = -Inf
      output[n, i_D_other] = 0
    }
  }
  
  
  
  return(output)
}


