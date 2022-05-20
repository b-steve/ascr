#' Extract coefficients from the output of ascr_tmb model
#'
#' @param object a fitted model from fit.ascr.tmb
#' @param types a character vector, accept any subset from <'all', 'fitted', 'linked', 'derived'>, default is 'linked'
#' @param pars a character vector containing any parameter names
#' @param new_covariates a data frame containing the values of covariates of any extended parameter
#'
#' @return a named numeric vector
#' @export
#'
coef.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, ...){
  #source('get_funs.r', local = TRUE)
  #source('support_functions.r', local = TRUE)
  
  #deal with default setting for 'types'
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))){
    stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  } 
  
  if('esa' %in% pars){
    pars = pars[-which(pars == 'esa')]
    #if pars = 'esa' only, then regard it as NULL and add 'derived' to 'types'
    if(length(pars) == 0){
      pars = NULL
    } else {
      #if also other parameters be assigned, we need to add the default setting for 'types' here
      #because after this step 'types' will not be NULL anyway, the default setting later will not work
      types = c(types, 'linked')
    }
    types = c(types, 'derived')
  }
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  if(is.null(types)) types = 'linked'
  
  is.fitted = 'fitted' %in% types
  is.linked = 'linked' %in% types
  is.derived = 'derived' %in% types
  
  
  #extract some key information for the object
  param_values_og = get_coef(object)
  name_og = get_param_og(object)
  name_extend = get_par_extend_name(object)
  df_param = get_data_param(object)
  esa = get_esa(object)
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  output = NULL
  
  for(i in pars){
    values_link = param_values_og[[i]]
    if(!i %in% name_extend){
      names(values_link) = paste(i, "link", sep = "_")
    } else {
      names(values_link) = gsub(" _ ", "\\.", names(values_link))
      names(values_link) = paste(names(values_link), "link", sep = "_")
    }
    
    if(is.linked){
      output = c(output, values_link)
    }
    
    if(is.fitted){
      par_info = subset(df_param, par == i)
      link = par_info$link
      if(i %in% name_extend){
        if(is.null(new_covariates)){
          #if there is no new covariates assigned and the parameter is extended, then regards all of the
          #relevant beta as identity linked parameters
          values_fitted = values_link
          names(values_fitted) = gsub("_link", "", names(values_fitted))
          link = 'identity'
        } else {
          if(nrow(new_covariates) != 1){
            stop('Argument "new_covariates" can only accept 1 row.')
          } else {
            gam = get_gam(object, i)
            tem = try({values_fitted = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask, values_link, new_covariates)})
            if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided.
                                          Defaulty all parameters are assigned as "par".
                                          And please make sure all categorical variables do not contain any new category.')
            names(values_fitted) = i
          }
        }
      } else {
        values_fitted = values_link
        names(values_fitted) = i
      }
      values_fitted = unlink.fun(link = link, value = values_fitted)
      output = c(output, values_fitted)
    }
  }
  
  if(is.derived){
    value_esa = esa$value
    names(value_esa) = paste('esa', seq(length(value_esa)), sep = '.')
    output = c(output, value_esa)
  }
  
  class(output) = 'coef_ascr_tmb'
  
  return(output)
}

##################################################################################################################

#' S3 method for output of coef.ascr_tmb
#'
#' @export
#'
print.coef_ascr_tmb = function(x){
  output = as.matrix(x)
  colnames(output) = "Est"
  print(output)
}

##################################################################################################################

#' Extract variance covariance matrix of the estimated parameters from ascr.tmb models
#'
#' @param object a fitted model from fit.ascr.tmb
#' @param types a character vector, accept any subset from <'all', 'fitted', 'linked', 'derived'>, default is 'linked'
#' @param pars a character vector containing any parameter names
#' @param new_covariates a data frame containing the values of covariates of any extended parameter
#'
#' @return a list with matrices as its elements if multiple 'types'. a matrix if only one 'types'.
#' @export
vcov.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, ...){
  #source('get_funs.r', local = TRUE)
  #source('support_functions.r', local = TRUE)
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))) stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  if('esa' %in% pars){
    pars = pars[-which(pars == 'esa')]
    #if pars = 'esa' only, then regard it as NULL and add 'derived' to 'types'
    if(length(pars) == 0){
      pars = NULL
    } else {
      #if also other parameters be assigned, we need to add the default setting for 'types' here
      #because after this step 'types' will not be NULL anyway, the default setting later will not work
      types = c(types, 'linked')
    }
    types = c(types, 'derived')
  }
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  
  if(is.null(types)) types = 'linked'
  
  # 'og' below means original output from the object
  param_values_og = get_coef(object)
  cov_og = object$vcov
  name_dim_og = colnames(cov_og)
  name_dim = gsub(" _ ", "\\.",name_dim_og)
  which.derived = which(substr(name_dim, 1, 3) == 'esa')
  name_dim = name_dim[-which.derived]
  name_dim_og = name_dim_og[-which.derived]
  cov_derived = cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked = cov_og[-which.derived, -which.derived, drop = FALSE]
  dimnames(cov_linked) = list(name_dim, name_dim)
  name_extend = get_par_extend_name(object)
  
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  #the original name is needed, in "name_og", "D.(Intercept)", "D.x1", "D.xxx" will be regarded as 'D' only
  #and the same to the other extendable parameters
  
  #the 'og' here means the original parameters name
  name_og = ori_name(name_dim)
  
  #extract the link functions for all element in the linked covariance matrix
  #if new_covariates is provided, this "link_funs" will be useless as cov matrix of fitted and
  #cov matrix of linked will have different dimentions
  df_param = get_data_param(object)
  link_funs = character(length(name_og))
  #extracted linked estimated values, the reason why not doing do.call('c', param_values_og)
  #is that, the fixed parameters could screw this up
  param_values = numeric(length(name_og))
  
  for(i in 1:length(name_dim)){
    param_values[i] = param_values_og[[name_og[i]]][name_dim_og[i]]
    if(!name_og[i] %in% name_extend){
      name_dim[i] = name_og[i]
      link_funs[i] = df_param[which(df_param$par == name_og[i]), 'link']
    } else {
      link_funs[i] = 'identity'
    }
  }
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = unique(name_og)
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters, which is not fixed, in this model.")
  }
  
  #obtain the indices for the assigned "pars"
  index_par = which(name_og %in% pars)
  name_og = name_og[index_par]
  name_dim = name_dim[index_par]
  cov_linked = cov_linked[index_par, index_par, drop = FALSE]
  link_funs = link_funs[index_par]
  param_values = param_values[index_par]
  
  types = unique(types)
  output = vector('list', length(types))
  names(output) = types
  
  for(type in types){
    if(type == "linked"){
      output[[type]] = cov_linked
      tem = paste(name_dim, 'link', sep = '_')
      dimnames(output[[type]]) = list(tem, tem)
    } else if(type == 'derived'){
      output[[type]] = cov_derived
    } else if(type == 'fitted'){
      #if 'new_covariates' is not provided, just keep the extended covariates as parameters with identity link function
      if(is.null(new_covariates)){
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, link_funs = link_funs)
        dimnames(output[[type]]) = list(name_dim, name_dim)
      } else {
        gam.output = get_gam(object)
        tem = try({output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates, pars = pars,
                                                          name_og = name_og, name_extend = name_extend, df_param = df_param,
                                                          gam.output = gam.output)})
        if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided. Defaulty all parameters are assigned as "par".')
        dimnames(output[[type]]) = list(pars, pars)
      }
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}

##################################################################################################################

#' Extract standard errors of the estimated parameters from ascr.tmb models
#'
#' @param object a fitted model from fit.ascr.tmb
#' @param types a character vector, accept any subset from <'all', 'fitted', 'linked', 'derived'>, default is 'linked'
#' @param pars a character vector containing any parameter names
#' @param new_covariates a data frame containing the values of covariates of any extended parameter
#'
#' @return a named numeric vector
#' @export
stdEr.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, ...){
  output_vcov = vcov(object = object, types = types, pars = pars, new_covariates = new_covariates)
  if(is(output_vcov, 'list')){
    output = NULL
    for(i in names(output_vcov)){
      output = c(output, sqrt(diag(output_vcov[[i]])))
    }
  } else {
    output = sqrt(diag(output_vcov))
  }
  class(output) = 'std_ascr_tmb'
  return(output)
}

##################################################################################################################

#' S3 method for the output of stdEr.ascr_tmb
#'
#' @export
#'
print.std_ascr_tmb = function(x){
  output = as.matrix(x)
  colnames(output) = "Std"
  print(output)
}

##################################################################################################################

#' Extract confidence interval for ascr.tmb models
#'
#' @param object a fitted model from fit.ascr.tmb
#' @param types a character vector, accept any subset from <'all', 'fitted', 'linked', 'derived'>, default is 'linked'
#' @param level confident level, default is 0.95
#' @param method 'default' only
#' @param linked logical value. if set to be TRUE, inverse of linked function will be applied to parameters confident 
#'               intervals as the 'fitted' confidence interval. if FALSE, delta method will be applied to obtain the 
#'               'fitted' confidence interval
#' @param pars a character vector containing any parameter names
#' @param new_covariates a data frame containing the values of covariates of any extended parameter
#' @param qqplot Not available yet
#' @param ask Not available yet
#'
#' @return a matrix
#' @export
confint.ascr_tmb = function(object, types = NULL, level = 0.95, method = 'default',
                            linked = FALSE, pars = NULL, new_covariates = NULL,
                            qqplot = FALSE, ask = FALSE, ...){

  #modify here to make it compatible with boot
  
  if(method != 'default') stop('Please apply bootstrap to the output of the model before using other methods.')
  
  #############################################################################################
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))) stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  
  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }
  
  if('esa' %in% pars){
    pars = pars[-which(pars == 'esa')]
    #if pars = 'esa' only, then regard it as NULL and add 'derived' to 'types'
    if(length(pars) == 0){
      pars = NULL
    } else {
      #if also other parameters be assigned, we need to add the default setting for 'types' here
      #because after this step 'types' will not be NULL anyway, the default setting later will not work
      types = c(types, 'linked')
    }
    types = c(types, 'derived')
  }
  
  if(!is.null(new_covariates) & linked){
    warning('Argument "linked" is set to TRUE, "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  if(linked & (!'fitted' %in% types)) types = c(types, 'fitted')
  if(is.null(types)) types = 'linked'
  
  types = unique(types)
  
  stopifnot(all(level < 1 & level > 0))
  p_upper = 0.5 + 0.5 * level
  p_lower = 0.5 - 0.5 * level
  col_name = paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  q_upper = qnorm(p_upper)
  q_lower = qnorm(p_lower)
  
  #if linked = FALSE, use coef(fit, 'fitted') + 1.96 * stdEr(fit, 'fitted')
  #if linked = TRUE, use unlink.fun(coef(fit, 'linked') + 1.96 * stdEr(fit, 'linked'))
  
  if(linked){
    df_linked = confint_gaussian_cal(object = object, types = 'linked', pars = pars,
                                     new_covariates = NULL, q_lower = q_lower,
                                     q_upper = q_upper)
    output_fitted = as.matrix(df_linked[, c('lower', 'upper')])
    df_param = get_data_param(object)
    name_par_extend = get_par_extend_name(object)
    
    for(i in 1:nrow(output_fitted)){
      par_name = ori_name(df_linked$par[i])
      if(!par_name %in% name_par_extend){
        link = df_param[which(df_param$par == par_name), 'link']
        output_fitted[i,] = unlink.fun(link = link, value = output_fitted[i,])
      }
    }
    rownames(output_fitted) = gsub("_link", "", df_linked$par)
    colnames(output_fitted) = col_name
    
    #linked = TRUE is a special method for types='fitted', so after dealing with it, we
    #can remove 'fitted' from 'types'
    if('fitted' %in% types) types = types[-which(types == 'fitted')]
    
    if(length(types) != 0){
      df_other = confint_gaussian_cal(object = object, types = types, pars = pars,
                                      new_covariates = new_covariates, q_lower = q_lower,
                                      q_upper = q_upper)
      output_other = as.matrix(df_other[, c('lower', 'upper')])
      rownames(output_other) = df_other$par
      colnames(output_other) = col_name
    } else {
      output_other = NULL
    }
    
    output = rbind(output_fitted, output_other)
  } else {
    df_all = confint_gaussian_cal(object = object, types = types, pars = pars,
                                  new_covariates = new_covariates, q_lower = q_lower,
                                  q_upper = q_upper)
    output = as.matrix(df_all[, c('lower', 'upper')])
    rownames(output) = df_all$par
    colnames(output) = col_name
  }
  
  return(output)
}


#' Title
#'
#' @param object 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
AIC.ascr_tmb = function(object, k = 2){
  if(object$fit.freqs){
    message("NOTE: Use of AIC for this model relies on independence between locations of calls from the same animal, which may not be appropriate.")
    return(NA)
  }
  return(k * length(coef(object)) - 2 * object$loglik)
}





#' Title
#' 
#' @param object 
#' @param newdata 
#' @param session 
#' @param type 
#' @param se.fit 
#' @param confidence 
#' @param level 
#' @param linked 
#' @param ... 
#' 
#' @return
#' @export
#'
#' @examples
predict.ascr_tmb = function(object, newdata = NULL, session = NULL, type = 'link', se.fit = FALSE, confidence = FALSE, level = 0.95, linked = FALSE, ...){
  
  if(is.null(newdata) & is.null(session)) session = 1
  if(!is.null(newdata) & !is.null(session)){
    warning("argument 'newdata' and 'session' could not be assigned at the same time, 'session' will be ignored.")
    session = NULL
  }
  
  stopifnot(type %in% c('link', 'response'))
  
  if(linked){
    type = 'response'
    confidence = TRUE
  } 
  
  if(confidence){
    se.fit = TRUE
    stopifnot(all(level < 1 & level > 0))
  }
  
  
  
  par_info = get_data_param(object)
  par_info = subset(par_info, par == 'D')
  
  if('D' %in% get_par_extend_name(object)){
    
    if(is.null(newdata)){
      #use the selected session's masks and corresponding covariates as newdata
      
      mask = get_mask(object)
      if(is(mask, 'list')) mask = mask[[session]]
      colnames(mask) = c('x', 'y')
      newdata = as.data.frame(mask)
      cov_mask = get_par_extend_data(object)$mask
      if(!is.null(cov_mask)){
        if('session' %in% colnames(cov_mask)) cov_mask = cov_mask[cov_mask$session == session, ,drop = FALSE]
        if(any(c('x', 'y') %in% colnames(cov_mask))) cov_mask = cov_mask[, which(!colnames(cov_mask) %in% c('x', 'y')),drop = FALSE]
        stopifnot(nrow(mask) == nrow(cov_mask))
        newdata = cbind(newdata, cov_mask)
      }
      
      cov_session = get_par_extend_data(object)$session
      if(!is.null(cov_session)){
        cov_session = cov_session[cov_session$session == session, ,drop = FALSE]
        cov_session = cov_session[,which(colnames(cov_session) != 'session'), drop = FALSE]
        stopifnot(nrow(cov_session) == 1)
        newdata = cbind(newdata, cov_session)
      }
    }
    
    
    gam.model = get_gam(object, 'D')
    values_link = as.vector(coef(object, types = 'linked', pars = 'D'))

    tem = get_extended_par_value(gam.model, par_info$n_col_full, par_info$n_col_mask, values_link, newdata, DX_output = TRUE)
    if(type == 'link'){
      est = tem$output
    } else if(type == 'response'){
      est = unlink.fun(link = par_info$link, value = tem$output)
    }
    
    if(linked) est_link = tem$output
      
    if(se.fit){
      DX = tem$DX
      vcov_matrix = vcov(object, types = 'linked', pars = 'D')
      if(type == 'link' | linked){
        log_scale = TRUE
      } else if(type == 'response' & !linked){
        log_scale = FALSE
      }
      std = sqrt(delta_for_pred(DX, values_link, vcov_matrix, log_scale))
    }
    
    
  } else {
    #when D is not extended, it is just a constant
    if(type == 'link'){
      tem = 'linked'
    } else if(type == 'response'){
      tem = 'fitted'
    }
    
    if(linked) est_link = coef(object, types = 'linked', pars = 'D')
    est <- coef(object, types = tem, pars = 'D')
    if(se.fit){
      if(!linked){
        std = stdEr(object, types = tem, pars = 'D')
      } else {
        std = stdEr(object, types = 'linked', pars = 'D')
      }
      
    }
    
  }
  output = est
  col_name = c('est')
  if(se.fit & !linked){
    col_name = c(col_name, 'std')
    output = cbind(est, std)
    colnames(output) = col_name
  }
  
  
  if(confidence){
    p_upper = 0.5 + 0.5 * level
    p_lower = 0.5 - 0.5 * level
    
    q_upper = qnorm(p_upper)
    q_lower = qnorm(p_lower)
    if(!linked){
      upper_lim = est + q_upper * std
      lower_lim = est + q_lower * std
    } else {
      upper_lim = unlink.fun(link = par_info$link, value = est_link + q_upper * std)
      lower_lim = unlink.fun(link = par_info$link, value = est_link + q_lower * std)
    }
    
    
    col_name = c(col_name, paste(c(p_lower, p_upper) * 100, "%", sep = " "))
    output = cbind(output, lower_lim)
    output = cbind(output, upper_lim)
    colnames(output) = col_name

  }
  
  return(output)
}





#predict density with information of locations coordinates, only used in the plot currently
#this function is a little bit ambiguous that I'm not sure this function should be here or the "support_functions.R"
#leave it here temporarily.
predict_with_location = function(fit, session_select = 1, new_data = NULL, D_cov = NULL, xlim = NULL, ylim = NULL,
                                  x_pixels = 50, y_pixels = 50, se_fit = FALSE, log_scale = FALSE, set_zero = NULL, 
                                  control_convert_loc2mask = NULL, ...){
  
  
  if(!is.null(set_zero) & se_fit){
    warning('when se_fit is TRUE, set_zero cannot be assigned.')
    set_zero = NULL
  }
  
  original_mask = FALSE
  
  if(is.null(new_data)){
    
    #if new_data is not provided, but xlim and ylim provided, we use xlim and ylim to build mask instead of 
    #extracting mask from fit
    if(any(!is.null(xlim), !is.null(ylim))){
      if(!all(!is.null(xlim), !is.null(ylim))) stop('please provide both xlim and ylim.')
      x = seq(from = xlim[1], to = xlim[2], length.out = x_pixels)
      y = seq(from = ylim[1], to = ylim[2], length.out = y_pixels)
      mask = data.frame(x = rep(x, each = y_pixels), y = rep(y, x_pixels))
    } else {
      mask = get_mask(fit)[[session_select]]
      original_mask = TRUE
    }
    
  } else {
    #if new_data is provided, then xlim and ylim will just do what they are supposed to do,
    #to trim the plot instead of building "mask"
    stopifnot(any(is(new_data, 'data.frame'), is(new_data, 'matrix')))
    stopifnot(all(c('x', 'y') %in% colnames(new_data)))
    mask = new_data[, c('x', 'y')]
  }
  
  
  if('D' %in% get_par_extend_name(fit)){
    
    #build the old_covariates
    ##we interpolate it again no matter there is new mask grid or not because in theory, user
    ##could use the same mask grid but different control_convert_loc2mask
    
    old_covariates = as.data.frame(mask)
    
    old_loc_cov = get_loc_cov(fit)
    
    #when we cannot find loc_cov from the model fitting object, it is possible the D is
    #not mask-level extended, it is also possible that the model fitting object comes
    #from simulation_study/demo, in the demo, we directly assign mask-level covariates
    #on each mask point
    mask_level_dat_extract = FALSE
    if(is.null(old_loc_cov)){
      old_loc_cov = get_par_extend_data(fit)$mask
      if(!is.null(old_loc_cov)){
        mask_level_dat_extract = TRUE
        if('session' %in% colnames(old_loc_cov)){
          old_loc_cov = subset(old_loc_cov, session == session_select)
        }
        
        #when locations/mask in prediction is assigned by xlim/ylim or new_data, we use
        #the original mask_level data as loc_cov, otherwise, we do not need to do
        #conversion in the first place
        if(!original_mask){
          
          if('session' %in% colnames(old_loc_cov)){
            old_loc_cov = old_loc_cov[,which(colnames(old_loc_cov)!='session'),drop = FALSE]
          }
          if('mask' %in% colnames(old_loc_cov)){
            old_loc_cov = old_loc_cov[,which(colnames(old_loc_cov)!='mask'),drop = FALSE]
          }
          
          tem_mask = get_mask(fit)
          if(is(tem_mask, 'list')) tem_mask = tem_mask[[session_select]]
          stopifnot(nrow(old_loc_cov) == nrow(tem_mask))
          old_loc_cov = cbind(tem_mask, old_loc_cov)
          
        }
      }

    }
    
    
    if(!is.null(old_loc_cov)){
      #like mentioned right above, if we are using original mask data, and mask-level covariates data
      #we do not need to do conversion at all
      if(mask_level_dat_extract & original_mask){
        cov_mask = old_loc_cov
      } else {
        if(is.null(control_convert_loc2mask)){
          control_convert_loc2mask = vector('list', 2)
          names(control_convert_loc2mask) = c('mask', 'loc_cov')
        }
        control_convert_loc2mask$mask = list(mask)
        control_convert_loc2mask$loc_cov = old_loc_cov
        
        cov_mask = do.call('location_cov_to_mask', control_convert_loc2mask)
      }
      
      if(any(colnames(cov_mask) %in% c('session', 'mask'))){
        old_covariates = cbind(old_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      } else {
        old_covariates = cbind(old_covariates, cov_mask)
      }
      
    }
    
    
    
    session_data_in_model = get_par_extend_data(fit)$session
    ##for session related covariates, it is simple, just take them from the input argument of fit
    
    if(!is.null(session_data_in_model)){
      tem = session_data_in_model[which(session_data_in_model$session == session_select), , drop = FALSE]
      tem = tem[, -which(colnames(tem) == 'session'), drop = FALSE]
      for(i in colnames(tem)) old_covariates[[i]] = tem[1,i]
    }
    
    
    if(!is.null(D_cov) | !is.null(new_data)){
      #firstly deal with all scenarios that there may be any new covariate provided 
      if(!is.null(D_cov)){
        stopifnot(is(D_cov, 'list'))
        stopifnot(any(c('session', 'location') %in% names(D_cov)))
      }
      
      #considering the possibility that user may only want to change part of covariates and remains other as the same
      #as the model, for example, 3 covariates related to D, 'weather' (session related), 'noise' and 'forest_type'(loc related)
      #and user only want to change weather, or to change noise, or to change noise and forest_type, and etc.
      
      #so we create 2 data frame, one for all new covariates, and one for old covariates. We must make sure these two data frame
      #contains the same mask points, and then replace the covariates in the "old" data frame with the same column in the "new",
      #if the same covariate appears in the "new" data frame.
      if(!is.null(new_data)){
        #in new_data, user could include any location related covariates directly, and it also contains x and y,
        #so we could directly use it instead of mask
        new_covariates = as.data.frame(new_data)
      } else {
        new_covariates = as.data.frame(mask)
      }
      
      
      
      #build the new_covariates based on all information we could have
      if(!is.null(D_cov$location)){
        if(is.null(control_convert_loc2mask)){
          control_convert_loc2mask = vector('list', 2)
          names(control_convert_loc2mask) = c('mask', 'loc_cov')
        }
        control_convert_loc2mask$mask = list(mask)
        control_convert_loc2mask$loc_cov = D_cov$location
        
        
        cov_mask = do.call('location_cov_to_mask', control_convert_loc2mask)
        new_covariates = cbind(new_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      }
      
      
      #since we only plot one session, the number of row for D_cov$session should be only 1
      if(!is.null(D_cov$session)){
        stopifnot(nrow(D_cov$session) == 1)
        for(i in colnames(D_cov$session)) new_covariates[[i]] = D_cov$session[1,i]
        
      }
      
      #update the old_covariates by the new_covariates
      
      for(i in colnames(old_covariates)){
        if(all(i != 'x', i != 'y', i %in% colnames(new_covariates))){
          old_covariates[[i]] = new_covariates[[i]]
        }
      }
    }
    
    par_info = get_data_param(fit)
    par_info = subset(par_info, par == 'D')
    gam.model = get_gam(fit, 'D')
    values_link = as.vector(coef(fit, types = 'linked', pars = 'D'))
    if(!is.null(set_zero)) values_link[set_zero] = 0

    tem = get_extended_par_value(gam.model, par_info$n_col_full, par_info$n_col_mask, values_link, old_covariates, DX_output = TRUE)
    D.mask = unlink.fun(link = par_info$link, value = tem$output)
    
    if(se_fit){
      DX = tem$DX
      vcov_matrix = vcov(fit, types = 'linked', pars = 'D')
      D.se = sqrt(delta_for_pred(DX, values_link, vcov_matrix, log_scale))
    }
    
    
  } else {
    #when D is not extended, it is literally a constant, just take the 1st estimated D from session 1
    #if D varies between sessions, it is extended, so there is no problem to just take session1
    D.mask <- rep(fit$D.mask[[1]][1], nrow(mask))
    if(se_fit){
      if(log_scale){
        type = "linked"
      } else {
        type = "fitted"
      }
      D.se = as.vector(rep(stdEr(fit, types = type, pars = 'D'), nrow(mask)))
    }
    
  }
  
  output = as.data.frame(mask)
  output$est = D.mask
  if(se_fit) output$std = D.se
  
  
  return(output)
}



#' Title
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
summary.ascr_tmb = function(object, ...){
  coefs = coef(object, types = 'fitted')
  derived = coef(object, types = 'derived')
  coefs_se = stdEr(object, types = 'fitted')
  derived_se = stdEr(object, types = 'derived')
  infotypes = get_infotypes(object)
  detfn = get_detfn(object)
  n.sessions = object$n.sessions
  
  pars_ext = get_par_extend_name(object)
  
  if(!is.null(pars_ext)){
    df_parm = get_data_param(object)
    df_parm = df_parm[df_parm$par %in% pars_ext, , drop = FALSE]
    pars_ext_links = df_parm$link
    names(pars_ext_links) = df_parm$par
  } else {
    pars_ext_links = NULL
  }
  
  if('ss' %in% infotypes){
    ss_link = get_ss_link(object)
    ss_cutoff = get_ss.opts(object)$cutoff
    ss_opts = list(cutoff = ss_cutoff, link = ss_link)
  } else {
    ss_opts = NULL
  }
  
  output = list(coefs = coefs, derived = derived,
                coefs_se = coefs_se, derived_se = derived_se,
                infotypes = infotypes, detfn = detfn,
                n.sessions = n.sessions, pars_ext_links = pars_ext_links,
                ss_opts = ss_opts)
  class(output) = c("summary_ascr_tmb", class(output))
  return(output)
}


#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.summary_ascr_tmb = function(x, ...){
  n.coefs <- length(x$coefs)
  n.derived <- length(x$derived)
  mat <- matrix(0, nrow = n.coefs + n.derived + 1, ncol = 2)
  mat[1:n.coefs, 1] <- as.vector(x$coefs)
  mat[1:n.coefs, 2] <- as.vector(x$coefs_se)
  mat[n.coefs + 1, ] <- NA
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), ] <- c(x$derived, x$derived_se)
  rownames(mat) <- c(names(x$coefs), "---", names(x$derived))
  colnames(mat) <- c("Estimate", "Std. Error")
  detfn <- c(hn = "Halfnormal", hhn = "Hazard halfnormal", hr = "Hazard rate", th = "Threshold",
             lth = "Log-link threshold", ss = "Signal strength")[x$detfn]
  infotypes <- c(bearing = "Bearings", dist = "Distances", ss = "Signal strengths",
                 toa = "Times of arrival", mrds = "Exact locations")[x$infotypes]
  
  pars_ext_links = x$pars_ext_links
  ss_opts = x$ss_opts
  cutoff = ss_opts$cutoff
  ss_link = ss_opts$link
  
  cat("Detection function:", detfn, "\n")
  cat("Number of sessions:", x$n.sessions, "\n")
  cat("Information types: ")
  cat(infotypes, sep = ", ")

  cat("\n", "\n", "Parameters:", "\n")
  stats::printCoefmat(mat, na.print = "")
  
  if(!is.null(pars_ext_links)){
    cat("\n", "\n", "Extended parameters link functions:", "\n")
    for(i in names(pars_ext_links)){
      cat(i, ":", pars_ext_links[i], "\n")
    }
  }
  
  if(!is.null(ss_opts)){
    cat("\n", "Signal strength related information:", "\n")
    cat("cutoff: ", cutoff, "\n")
    cat("link function:", ss_link, "\n")
  }
}
