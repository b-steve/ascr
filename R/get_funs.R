get_D_tmb = function(fit){
  
  #estimated population density for each mask is recorded in fit$D.mask
  #extract the values from it and make it into a usable format data frame
  #with session and mask indices
  
  dims = fit$output.tmb$dims
  output = data.frame(session = rep(1:dims$n.sessions, dims$n.masks))
  tem = vector('list', dims$n.sessions)
  for(s in 1:dims$n.sessions){
    tem[[s]] = data.frame(mask = 1:dims$n.masks[s], D_est = fit$D.mask[[s]])
  }
  tem = do.call('rbind', tem)
  output = cbind(output, tem)
  row.names(output) = NULL
  data_mask = get_data_mask(fit)
  output = merge(output, data_mask, by = c('session', 'mask'))
  output = output[order(output$session, output$mask),]
  return(output)
}

get_esa = function(fit){
  output = fit$esa
  return(output)
}

get_data_full = function(fit){
  output = fit$output.tmb$data.full
  return(output)
}

get_data_mask = function(fit){
  output = fit$output.tmb$data.mask
  return(output)
}

get_data_trap = function(fit){
  output = fit$output.tmb$data.traps
  return(output)
}

get_par_extend = function(fit){
  output = fit$args$par.extend
  return(output)
}

get_par_extend_name = function(fit){
  output = fit$output.tmb$param.extend
  return(output)
}

get_area_unit_tmb = function(fit){
  output = fit$output.tmb$area_unit
  return(output)
}

get_dims_tmb = function(fit){
  output = fit$output.tmb$dims
  return(output)
}

get_cue_rates = function(fit){
  output = fit$output.tmb$avg_cue_rates
  return(output)
}

get_survey_length = function(fit){
  output = fit$args$survey.length
  return(output)
}

get_dist_theta = function(fit){
  output = fit$output.tmb$data.dists.thetas
  return(output)
}

get_detfn = function(fit){
  output = fit$output.tmb$detfn
  return(output)
}

get_DX_full = function(fit){
  output = fit$output.tmb$DX$DX_non_mask
  return(output)
}

get_DX_mask = function(fit){
  output = fit$output.tmb$DX$DX_mask
  return(output)
}

get_data_param = function(fit){
  output = fit$output.tmb$param.info.table
  return(output)
}

get_coef = function(fit){
  output = fit$output.tmb$coef_link
  return(output)
}

get_param_og = function(fit){
  output = fit$output.tmb$param.og
  return(output)
}

get_ss_link = function(fit){
  output = fit$output.tmb$ss.link
  return(output)
}

get_cutoff = function(fit){
  output = fit$output.tmb$cutoff
  return(output)
}

get_buffer = function(fit){
  output = numeric(fit$n.sessions)
  for(s in 1:fit$n.sessions){
    buffer = attr(fit$args$mask[[s]], 'buffer')
    if(!is.null(buffer))
      output[s] = buffer
  }
  return(output)
}



get_sound_speed = function(fit){
  output = fit$output.tmb$sound.speed
  return(output)
}

get_fit_type = function(fit){
  output = fit$fit.types
  return(output)
}

get_gam = function(fit, param = NULL, which.level = NULL){
  if(!is.null(param)){
    output = fit$output.tmb$gam_output[[param]]
    if(!is.null(which.level)){
      output = output[[which.level]]
    }
  } else {
    output = fit$output.tmb$gam_output
    if(!is.null(which.level)){
      for(i in names(output)){
        output[[i]] = output[[i]][[which.level]]
      }
    }
  }
  
  return(output)
}

get_DX_new_gam = function(mod, newdata){
  newdata = newdata[, names(mod$var.summary), drop = FALSE]
  
  newdata$gam.resp = 1
  olddata = mod$mf
  dat = rbind(olddata, newdata)
  output = mgcv::gam(mod$formula, data = dat, fit = FALSE)$X
  output = output[(nrow(olddata) + 1):nrow(dat), , drop = FALSE]
  return(output)
}


get_extended_par_value = function(gam, n_col_full, n_col_mask, par_value_linked, new_covariates){
  if(n_col_full > 1){
    gam.model = gam$gam_non_mask
    DX_full_new = get_DX_new_gam(gam.model, new_covariates)
  } else {
    DX_full_new = matrix(1, ncol = 1, nrow = nrow(new_covariates))
  }
  
  if(n_col_mask > 0){
    gam.model = gam$gam_mask
    DX_mask_new = get_DX_new_gam(gam.model, new_covariates)
    #get rid of the first column since the intercept is not here
    DX_mask_new = DX_mask_new[, -1, drop = FALSE]
  } else {
    DX_mask_new = NULL
  }
  
  DX_new = cbind(DX_full_new, DX_mask_new)
  output = as.vector(DX_new %*% as.matrix(par_value_linked, ncol = 1))
  return(output)
}

