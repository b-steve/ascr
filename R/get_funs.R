get_D_tmb = function(fit){
  
  #estimated population density for each mask is recorded in fit$D.mask
  #extract the values from it and make it into a usable format data frame
  #with session and mask indices
  
  dims = get_dims_tmb(fit)
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
  output = sort.data(output, 'data.mask')
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

get_mask = function(fit){
  output = fit$args$mask
  return(output)
}

get_data_trap = function(fit){
  output = fit$output.tmb$data.traps
  return(output)
}

get_trap = function(fit){
  df.traps = get_data_trap(fit)
  colnames(df.traps)[1:2] = c('x', 'y')
  n.sessions = max(df.traps$session)
  output = vector('list', n.sessions)
  
  for(s in 1:n.sessions){
    output[[s]] = as.matrix(subset(df.traps, session == s, select = c('x', 'y')))
  }
  
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

get_ss.opts = function(fit){
  output = fit$args$ss.opts
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

get_scale_cov = function(fit){
  output = fit$scale.covs
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

#for numeric columns in the data frames in par.extend$data, the mean will be
#used as the default value, and for non-numeric columns, the first element
#will be used as the default value.
#Actually if the argument for the fit.ascr(par.extend = xxx) sets scale = TRUE
#then the default values for numeric columns will be just 0, however, if scale = FALSE
#this function could help

get_default_covariates = function(fit){
  dat = fit$args$par.extend$data
  dat_name = names(dat)
  o = vector('list', length(dat_name))
  names(o) = dat_name
  
  for(i in dat_name){
    tem = dat[[i]]
    tem = tem[, which(!colnames(tem) %in% c('session','ID','trap','mask')), drop = FALSE]
    tem = as.data.frame(apply(tem, 2, mean_diy, simplify = FALSE))
    tem = fit$scale.covs(tem)
    o[[i]] = tem
  }
  
  output = o[[dat_name[1]]]
  if(length(dat_name) > 1){
    for(i in 2:length(dat_name)) output = cbind(output, o[[dat_name[i]]])
  }
  return(output)
}

get_sv_for_boot = function(fit){
  #depend on whether there is any parameter be extended in this model,
  #it can be determined whether use default covariates (mean for each covariate) values 
  is.extended = !is.null(get_par_extend_name(fit))
  if(is.extended){
    covaraites_default = get_default_covariates(fit)
    output = as.list(coef(fit, 'fitted', new_covariates = covaraites_default))
  } else {
    output = as.list(coef(fit, 'fitted'))
  }
  
  df_parm = get_data_param(fit)
  fixed_par = names(fit$args$fix)

  #avoid extreme start values if it is not fixed
  for(i in names(output)){
    if(!i %in% fixed_par){
      link = df_parm[which(df_parm$par == i), 'link']
      if(link == 'logit') output[[i]] = max(min(0.95, output[[i]]), 0.05)
      if(link == 'log') output[[i]] = max(0.05, output[[i]])
    }
  }
  
  return(output)
}

#a simple version of "create.capt()", input is the simulated capture history if there is any,
#which means nrow(capture) > 0
get_capt_for_boot = function(captures, dims, infotypes){
  #dims$n.IDs is not valid here, do not use it
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps

  #deal with irregular IDs and remove useless info types
  captures = captures[,c('session', 'ID', 'trap', infotypes)]
  tem_capt = vector('list', n.sessions)
  n.IDs = numeric(n.sessions)
  tem_data_head = vector('list', n.sessions)
  for(s in 1:n.sessions){
    tem_capt[[s]] = subset(captures, captures$session == s)
    if(nrow(tem_capt[[s]] > 0)){
      tem_capt[[s]]$ID = as.numeric(as.factor(tem_capt[[s]]$ID))
      n.IDs[s] = max(tem_capt[[s]]$ID)
      tem_data_head[[s]] = data.frame(session = s, ID = rep(seq(n.IDs[s]), each = n.traps[s]), 
                                      trap = rep(seq(n.traps[s]), n.IDs[s]))
    } else {
      n.IDs[s] = 0
    }

  }
  captures = do.call('rbind', tem_capt)
  captures$bincapt = 1
  data_head = do.call('rbind', tem_data_head)
  
  o = merge(data_head, captures, by = c('session', 'ID', 'trap'), all.x = TRUE)
  o = sort.data(o, 'data.full')
  for(i in c('bincapt', infotypes)) o[,i] = ifelse(is.na(o[,i]), 0, o[,i])
  #prepare the structure of the output
  if(n.sessions == 1){
    output = vector('list', length(infotypes) + 1)
    names(output) = c('bincapt', infotypes)
    for(i in names(output)){
      output[[i]] = matrix(o[,i], ncol = n.traps, byrow = TRUE)
      rownames(output[[i]]) = seq(n.IDs)
    }
  } else {
    output = vector('list', n.sessions)
    for(s in 1:n.sessions){
      output[[s]] = vector('list', length(infotypes) + 1)
      names(output[[s]]) = c('bincapt', infotypes)
      tem = subset(o, o$session == s)
      for(i in names(output[[s]])){
        if(n.IDs[s] > 0){
          output[[s]][[i]] = matrix(tem[,i], ncol = n.traps[s], byrow = TRUE)
          rownames(output[[s]][[i]]) = seq(n.IDs[s])
        } else {
          output[[s]][[i]] = matrix(nrow = 0, ncol = n.traps[s])
        }
      }
    }
  }
  
  return(output)
  
}

