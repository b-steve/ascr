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

  return(output)
}

get_loc_cov = function(fit){
  output = fit$arg_input$loc_cov
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

get_mask_from_data = function(dat){
  output = dat$mask
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

get_trap_from_data = function(dat){
  output = dat$traps
  return(output)
}

get_par_extend = function(fit){
  output = fit$args$par.extend
  return(output)
}

get_par_extend_data = function(fit){
  output = get_par_extend(fit)
  output = output$data
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

get_infotypes = function(fit){
  output = fit$infotypes
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
  var_names = names(mod$var.summary)
  try({newdata = newdata[, var_names, drop = FALSE]})
  #make sure the new data are in the right data type
  for(i in var_names){
    if(is.null(newdata[[i]])) stop(paste0('covariate: ', i, " is not provided."))
    if(is(mod$var.summary[[i]], 'factor')) newdata[[i]] = as.character(newdata[[i]])
    if(is(mod$var.summary[[i]], 'numeric')) newdata[[i]] = as.numeric(newdata[[i]])
  } 
  
  #construct the "old data", make sure it has all levels for all categorical variables
  olddata = vector('list', length(var_names))
  names(olddata) = var_names
  for(i in var_names){
    if(is(mod$var.summary[[i]], 'factor')){
      #var.summary[[i]] is a factor is the explanatory variable is a categorical one
      #and the levels() returns a character vector, just what we want
      #one potential problem, if the gam object was built by the data which contains "factor" class
      #column in the first place, and it contains for example 5 levels, but in the data, this column
      #only has 4 levels. However, I guess this scenario is extremely unlikely to happen.
      olddata[[i]] = levels(mod$var.summary[[i]])
    } else if(is(mod$var.summary[[i]], 'numeric')){
      olddata[[i]] = 1
    } else {
      #actually I never see any other class, just leave an "in-case" step here
      stop(paste0('unknown data type in one of gam model objects, the variable name is: ', i, 
                  ', please check the relevant data.'))
    }
  }
  
  #adjust the length, make sure all element in this list have the same length
  len = max(sapply(olddata, 'length'))
  
  for(i in var_names) olddata[[i]] = rep(olddata[[i]], length = len)
  
  
  olddata = as.data.frame(olddata)
  
  dat = rbind(olddata, newdata)
  dat$gam.resp = 1
  output = mgcv::gam(mod$formula, data = dat, fit = FALSE)$X
  output = output[(nrow(olddata) + 1):nrow(dat), , drop = FALSE]
  return(output)
}

get_extended_par_value = function(gam, n_col_full, n_col_mask, par_value_linked, new_covariates,
                                  DX_output = FALSE, matrix_par_value = FALSE){
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
  if(!matrix_par_value){
    output = as.vector(DX_new %*% as.matrix(par_value_linked, ncol = 1))
  } else {
    output = par_value_linked %*% t(DX_new)
  }
  
  if(DX_output){
    return(list(output = output, DX = DX_new))
  } else {
    return(output)
  }
  
}

#for numeric columns in the data frames in par.extend$data, the mean will be
#used as the default value, and for non-numeric columns, the first element
#will be used as the default value.
#Actually if the argument for the fit.ascr(par.extend = xxx) sets scale = TRUE
#then the default values for numeric columns will be just 0, however, if scale = FALSE
#this function could help

get_default_covariates = function(fit){
  dat = get_par_extend_data(fit)
  dat_name = names(dat)
  o = vector('list', length(dat_name))
  names(o) = dat_name
  
  for(i in dat_name){
    tem = dat[[i]]
    tem = tem[, which(!colnames(tem) %in% c('session','ID','trap','mask')), drop = FALSE]
    tem = as.data.frame(apply(tem, 2, mean_diy, simplify = FALSE))
    o[[i]] = tem
  }
  
  output = o[[dat_name[1]]]
  if(length(dat_name) > 1){
    for(i in 2:length(dat_name)) output = cbind(output, o[[dat_name[i]]])
  }
  return(output)
}

get_fixed_par_name = function(fit){
  output = fit$output.tmb$param.fix
  if(length(output) == 0) output = NULL
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
get_capt_for_boot = function(captures, dims, infotypes, animal.model){
  #since captures is the simulated data instead of original capture history used in the original model fitting
  #dims$n.IDs and dims$animal_ID are not valid here, do not use it
  n.sessions = dims$n.sessions
  n.traps = dims$n.traps
  #deal with irregular IDs and remove useless info types
  captures = captures[,c('session', 'animal_ID'[animal.model], 'ID', 'trap', infotypes)]
  #convert ID and animal_ID to natural number to make the later process to build head data easier
  captures = convert_natural_number(captures, animal.model)
  captures$bincapt = 1
  
  #prepare a head data which contains all combinations of session-animal_ID[animal.model]-ID-trap

  
  tem_data_head = vector('list', n.sessions)
  for(s in 1:n.sessions){
    capt_session = subset(captures, session == s)
    if(nrow(capt_session) > 0){
      if(animal.model){
        n.animal = max(capt_session$animal_ID)
        tem = vector('list', n.animal)
        for(a in 1:n.animal){
          n.ID = max(subset(capt_session, animal_ID == a)$ID)
          tem[[a]] = data.frame(session = s, ID = rep(seq(n.ID), each = n.traps[s]), 
                                trap = rep(seq(n.traps[s]), n.ID), animal_ID = a)
        }
        
        tem_data_head[[s]] = do.call('rbind', tem)
        
      } else {
        n.ID = max(capt_session$ID)
        tem_data_head[[s]] = data.frame(session = s, ID = rep(seq(n.ID), each = n.traps[s]), 
                                        trap = rep(seq(n.traps[s]), n.ID))
        
      }
      
    } else {
      tem_data_head[[s]] = data.frame(session = s, ID = NA, trap = 1:n.traps[s])
      if(animal.model) tem_data_head[[s]]$animal_ID = NA
    }

  }

  data_head = do.call('rbind', tem_data_head)
  
  if(animal.model){
    o = merge(data_head, captures, by = c('session', 'ID', 'trap', 'animal_ID'), all.x = TRUE)
  } else {
    o = merge(data_head, captures, by = c('session', 'ID', 'trap'), all.x = TRUE)
  }
  
  o = sort.data(o, 'data.full')
  for(i in c('bincapt', infotypes)) o[,i] = ifelse(is.na(o[,i]) & !is.na(o[,'ID']), 0, o[,i])
  
  if(!animal.model){
    n.IDs = numeric(n.sessions)
    for(s in 1:n.sessions){
      tem = subset(captures, session == s)
      if(nrow(tem) > 0) n.IDs[s] = max(tem$ID)
    }
    
    
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
  } else {
    output = o
  }
  
  
  return(output)
  
}

#obtain capture history from the ascr_data object
get_capt_for_plot = function(dat){
  all.types <- c("bearing", "dist", "ss", "toa")
  
  capt = dat$capt
  
  #according to current data conversion function, if capt is a data frame, it means
  #the data is individual identification embedded data, or "animal_ID" is included.
  if(is(capt, 'data.frame')){
    #firstly convert animal_ID and ID to natural successive numbers
    data.capt = convert_natural_number(capt, TRUE, "both")
    
  } else {
    if("bincapt" %in% names(capt)){
      n.sessions = 1
      n.traps = ncol(capt$bincapt)
      n.IDs = nrow(capt$bincapt)
      extra_info <- all.types[all.types %in% names(capt)]
      is.mrds = "mrds" %in% names(capt)
    } else {
      n.sessions = length(capt)
      n.traps = sapply(capt, function(x) ncol(x$bincapt))
      n.IDs = sapply(capt, function(x) nrow(x$bincapt))
      extra_info <- all.types[all.types %in% names(capt[[1]])]
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
      for(j in c("trap", "bincapt", extra_info)) tem.df[[j]] = numeric(number.row)
      if(is.mrds) tem.df[,c('mrds_x', 'mrds_y')] = 0
      
      if(number.row > 0){
        tem.df$ID = rep(1:nrow(tem$bincapt), n.traps[i])
        tem.df$trap = rep(1:n.traps[i], each = n.IDs[i])
        for(k in c('bincapt', extra_info)) tem.df[[k]] = as.vector(tem[[k]])
        
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
    data.capt = convert_natural_number(data.capt, FALSE, "ID")
  }
  
  data.capt = sort.data(data.capt, "data.full")
  data.capt = subset(data.capt, bincapt == 1)
  data.capt = data.capt[, -which(colnames(data.capt) == 'bincapt')]
  
  return(data.capt)
}


get_boot_res = function(fit, pars){
  output = fit$boot$boots
  if(!is.null(pars)){
    name_og = ori_name(colnames(output))
    output = output[, which(name_og %in% pars), drop = FALSE]
  }
  
  return(output)
}

get_boot_res_esa = function(fit){
  output = fit$boot$res_esa
  colnames(output) = gsub("_", "\\.", colnames(output))
  return(output)
}


#get bias from a single bootstrapped matrix
get_bias <- function(res, coefs){
  f = function(x) mean(x, na.rm = T)
  bias = apply(res, 2, f) - coefs

  return(bias)
}

###################################################################################
#obtain a subset of coefs, which is generated by coef.ascr_tmb(), based on "type"
#currently is not used
get_coef_base_type = function(coefs, type){
  coef_name = names(coefs)
  
  if(type == 'linked'){
    output = coefs[grepl('_link$', coef_name)]
  } else if(type == 'fitted'){
    is.link = grepl('_link$', coef_name)
    is.derived = grepl('^esa', coef_name)
    output = coefs[!(is.link | is.derived)]

  } else {
    output = coefs[grepl('^esa', coef_name)]
  }
  
  return(output)
  
}


###################################################################################################