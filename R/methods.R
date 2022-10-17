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
  
  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  
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
  
  #since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  #to be consistent with the default order of all parameters
  pars = fulllist.par.generator()[fulllist.par.generator() %in% pars]
  
  
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
          #if(nrow(new_covariates) != 1) stop('Argument "new_covariates" can only accept 1 row.')
          #browser()
          gam = get_gam(object, i)
          tem = try({values_fitted = get_extended_par_value(gam, par_info$n_col_full,
                                                            par_info$n_col_mask, values_link, new_covariates)})
          if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided.
                                        Defaulty all parameters are assigned as "par".
                                        And please make sure all categorical variables do not contain any new category.')
          names(values_fitted) = rep(i, length(values_fitted))
          
        }
      } else {
        values_fitted = values_link
        names(values_fitted) = i
      }
      
      if(back_trans){
        values_fitted = unlink.fun(link = link, value = values_fitted)
      }
      
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


#' Title
#'
#' @param object 
#' @param types 
#' @param pars 
#' @param new_covariates 
#' @param correct_bias 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
coef.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, correct_bias = FALSE, ...){

  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  name_og = get_param_og(object)
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  #since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  #to be consistent with the default order of all parameters
  pars = fulllist.par.generator()[fulllist.par.generator() %in% pars]
  
  
  if(correct_bias){
    res = get_boot_res(object, pars)
    coefs = coef.ascr_tmb(object, types = 'linked', pars)
    bias = get_bias(res, coefs)
    coefs_linked_bias_corrected = coefs - bias
    output = vector('list', length(types))
    names(output) = types
    
    for(i in types){
      
      if(i == 'linked'){
        output[[i]] = coefs_linked_bias_corrected
      } else if (i == 'derived'){
        coefs_esa = coef.ascr_tmb(object, types = 'derived')
        res_esa = get_boot_res_esa(object)
        bias_esa = get_bias(res_esa, coefs_esa)
        output[[i]] = coefs_esa - bias_esa
        
      } else {
        #when types is 'fitted', use the 'coefs_linked_bias_corrected' and new_covariates to
        #calculate the estimations.
        
        #get foundation information needed for this method
        name_extend = get_par_extend_name(object)
        df_param = get_data_param(object)
        
        linked_name = names(coefs)
        original_name = ori_name(linked_name)
        
        
        for(j in pars){
          values_link = coefs_linked_bias_corrected[original_name == j]
          
          par_info = subset(df_param, par == j)
          link = par_info$link
          
          if(j %in% name_extend){
            if(is.null(new_covariates)){
              #if there is no new covariates assigned and the parameter is extended, then regards all of the
              #relevant beta as identity linked parameters
              values_fitted = values_link
              names(values_fitted) = gsub("_link", "", names(values_fitted))
              link = 'identity'
            } else {
              if(nrow(new_covariates) != 1) stop('Argument "new_covariates" can only accept 1 row.')
              
              gam = get_gam(object, j)
              tem = try({values_fitted = get_extended_par_value(gam, par_info$n_col_full,
                                                                par_info$n_col_mask, values_link, new_covariates)})
              if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided.
                                        Defaulty all parameters are assigned as "par".
                                        And please make sure all categorical variables do not contain any new category.')
              names(values_fitted) = j
              
            }
          } else {
            values_fitted = values_link
            names(values_fitted) = j
          }
          
          if(back_trans){
            values_fitted = unlink.fun(link = link, value = values_fitted)
          }
          
          output[[i]] = c(output[[i]], values_fitted)
        }
        #end of if(i == 'fitted')
      }
      
      #end of for i in types
    }

    names(output) = NULL
    output = do.call('c', output)
    
  } else {
    output = coef.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates,
                           back_trans = back_trans)
  }
  
  class(output) = "coef_ascr_tmb"
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
  
  #back_trans is a hidden argument that we may not allow user to use it, it is only used
  #for confident interval calculation. When new_covariates is provided, in order to calculate
  #"fitted" scale CI, we need delta method to calculate sd and build CI with new covariates included
  #but not back transformed. For example, D ~ D_beta1, we need sd(D_int + D_beta1 * x1) and build CI
  #of log(D). We don't need sd(exp(D_int + D_beta1 * x1)) in this scenario.
  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  
  
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
  
  
  #since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  #to be consistent with the default order of all parameters
  pars = fulllist.par.generator()[fulllist.par.generator() %in% pars]
  
  #obtain the indices for the assigned "pars"
  index_par = which(name_og %in% pars)
  name_og = name_og[index_par]
  name_dim = name_dim[index_par]
  cov_linked = cov_linked[index_par, index_par, drop = FALSE]
  link_funs = link_funs[index_par]
  param_values = param_values[index_par]
  
  output = vector('list', length(types))
  names(output) = types
  #browser()
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
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, link_funs = link_funs, back_trans = back_trans)
        dimnames(output[[type]]) = list(name_dim, name_dim)
      } else {
        gam.output = get_gam(object)
        tem = try({output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates, pars = pars,
                                                          name_og = name_og, name_extend = name_extend, df_param = df_param,
                                                          gam.output = gam.output, back_trans = back_trans)})
        if(is(tem, 'try-error')) stop('Please make sure all covariates needed for assigned "par" are provided.
                                      Defaulty all parameters are assigned as "par".')
        dimnames(output[[type]]) = list(pars, pars)
      }
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}


#' Title
#'
#' @param object 
#' @param types 
#' @param pars 
#' @param new_covariates 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
vcov.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, from_boot = TRUE, ...){
  #This hidden argument just left here in order to make this method be consistent with
  #vcov.ascr_tmb. Because we don't use standard error to build confidence interval
  #for boot strap cases in the first place, this is propbably useless.
  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  if(!from_boot){
    output = vcov.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates,
                           back_trans = back_trans)
    return(output)
  }
  
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  all_par_name_og = get_param_og(object)
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = all_par_name_og
  } else {
    if(any(!pars %in% all_par_name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  #since the loop below depends on the order of pars, we need to re-order 'pars' to make its order
  #to be consistent with the default order of all parameters
  pars = fulllist.par.generator()[fulllist.par.generator() %in% pars]

  #browser()
  if(any(c('linked', 'fitted') %in% types)){
    res = get_boot_res(object, pars)
    fixed_par = get_fixed_par_name(object)
    pars = pars[!pars %in% fixed_par]
  }
  
  

  output = vector('list', length(types))
  names(output) = types
  for(i in types){
    if(i == 'linked'){
      output[[i]] = var_from_res(res, fixed_par)
    } else if(i == 'derived'){
      #browser()
      res_esa = get_boot_res_esa(object)
      output[[i]] = var_from_res(res_esa)
    } else {
      res_fitted = res_transform(res, new_covariates, pars, object, back_trans)
      output[[i]] = var_from_res(res_fitted, fixed_par)
    }
  }
  
  #when only one types, directly output it without a list structure
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
  
  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  output_vcov = vcov.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates,
                              back_trans = back_trans)
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


#' Title
#'
#' @param object 
#' @param types 
#' @param pars 
#' @param new_covariates 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
stdEr.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, from_boot = TRUE, ...){
  extra_args = list(...)
  back_trans = extra_args$back_trans
  if(is.null(back_trans)){
    back_trans = TRUE
  }
  
  output_vcov = vcov.ascr_boot(object = object, types = types, pars = pars, new_covariates = new_covariates,
                               from_boot = from_boot, back_trans = back_trans)
  
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

######################################################################################################################

#' Extract confidence interval for ascr.tmb models
#'
#' @param object a fitted model from fit.ascr.tmb
#' @param types a character vector, accept any subset from <'all', 'fitted', 'linked', 'derived'>, default is 'linked'
#' @param level confident level, default is 0.95
#' @param pars a character vector containing any parameter names
#' @param new_covariates a data frame containing the values of covariates of any extended parameter
#' @param ... 
#'
#' @return a matrix
#' @export
confint.ascr_tmb = function(object, types = NULL, level = 0.95, pars = NULL, new_covariates = NULL, ...){

  #browser()
  #############################################################################################
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars

  
  stopifnot(all(level < 1 & level > 0))
  p_upper = 0.5 + 0.5 * level
  p_lower = 0.5 - 0.5 * level
  col_name = paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  q_upper = qnorm(p_upper)
  q_lower = qnorm(p_lower)
  
  #if linked = FALSE, use coef(fit, 'fitted') + 1.96 * stdEr(fit, 'fitted')
  #if linked = TRUE, use unlink.fun(coef(fit, 'linked') + 1.96 * stdEr(fit, 'linked'))
  

  output = confint_gaussian_cal(object = object, types = types, pars = pars,
                                new_covariates = new_covariates, q_lower = q_lower,
                                q_upper = q_upper)

  df_param = get_data_param(object)
  #browser()
  if('fitted' %in% types){
    tem = output[['fitted']]
    name_og = get_param_og(object)
    for(i in 1:nrow(tem)){
      p = tem[i, 'par']
      if(p %in% name_og){
        link = df_param[which(df_param$par == p), 'link']
      } else {
        link = 'identity'
      }
      
      tem[i, c('lower', 'upper')] = unlink.fun(link, tem[i, c('lower', 'upper')])
    }
    output[['fitted']] = tem
  }
  
  for(i in names(output)){
    par_name = output[[i]]$par
    output[[i]] = as.matrix(output[[i]][, c('lower', 'upper'), drop = FALSE])
    rownames(output[[i]]) = par_name
  }
  
  output = do.call('rbind', output)
  
  return(output)
}




#' Title
#'
#' @param object 
#' @param types 
#' @param level 
#' @param pars 
#' @param new_covariates 
#' @param correct_bias 
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
confint.ascr_boot = function(object, types = NULL, level = 0.95, pars = NULL, new_covariates = NULL,
                             correct_bias = FALSE, from_boot = TRUE, ...){
  if(!from_boot){
    if(correct_bias) message(paste0("The argument 'correct_bias' will be ignored as the confidence",
      " interval is not from bootstrap results."))
    
    output = confint.ascr_tmb(object = object, types = types, level = level, pars = pars,
                              new_covariates = new_covariates)
    return(output)
  }
  
  stopifnot(all(level < 1 & level > 0))
  name_og = get_param_og(object)
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  
  
  #check which parameter will be displayed
  if(is.null(pars)){
    pars = name_og
  } else {
    if(any(!pars %in% name_og)) stop("Argument 'pars' only accept parameters' name in this model.")
  }
  
  #make sure pars are in the right order
  pars = fulllist.par.generator()[fulllist.par.generator() %in% pars]

  
  #obtain the bootstrap result
  res = get_boot_res(object, pars)
  coef_link = coef.ascr_tmb(object, types = 'linked', pars)
  
  res = res_mod_for_CI(res, coef_link, correct_bias)

  
  #prepare for output
  output = vector('list', length(types))
  #deal with different output types one by one
  for(i in types){
    if(i == 'linked'){
      output[[i]] = boot_res_to_CI(res, level)
    } else if(i == 'derived'){
      res_esa = get_boot_res_esa(object)
      coef_esa = coef.ascr_tmb(object, types = 'derived')
      res_esa = res_mod_for_CI(res_esa, coef_esa, correct_bias)
      
      output[[i]] = boot_res_to_CI(res_esa, level)
    } else {
      res_fitted = res_transform(res, new_covariates, pars, object)
      output[[i]] = boot_res_to_CI(res_fitted, level)
      
    }
  }
  
  output = do.call('rbind', output)

  return(output)
}



######################################################################################################################


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
    message("NOTE: Use of AIC for this model relies on independence between locations of
            calls from the same animal, which may not be appropriate.")
    return(NA)
  }
  return(k * length(coef(object)) - 2 * object$loglik)
}



predict.ascr_tmb = function(fit, type = 'link', newdata = NULL, se.fit = TRUE, confidence = TRUE,
                            level = 0.95, pars = NULL, ...){
  
}



















#for reference only, delete this function after finishing predict method.
predict_D = function(object, newdata = NULL, session = NULL, type = 'link', se.fit = FALSE, 
                     confidence = FALSE, level = 0.95, ...){

  if(is.null(newdata) & is.null(session)) session = 1
  if(!is.null(newdata) & !is.null(session)){
    warning("arguments of 'newdata' and 'session' could not be assigned at the same time, 'session' will be ignored.")
    session = NULL
  }
  
  stopifnot(type %in% c('link', 'response'))
  
  
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
        if(any(c('x', 'y') %in% colnames(cov_mask))) cov_mask = cov_mask[, which(!colnames(cov_mask) %in% c('x', 'y')), drop = FALSE]
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
  
  #if argument 'linked' is TRUE, the relationship between std and CI is not obvious
  #to avoid confusion, here we do not show std under this scenario.
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
  CI = confint(object, types = 'fitted')
  CI_derived = confint(object, types = 'derived')
  is_boot = is(object, 'ascr_boot')
  
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
                CI = CI, CI_derived = CI_derived,
                is_boot = is_boot, infotypes = infotypes, detfn = detfn,
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
  mat <- matrix(0, nrow = n.coefs + n.derived + 1, ncol = 4)
  mat[1:n.coefs, 1] <- as.vector(x$coefs)
  mat[1:n.coefs, 2] <- as.vector(x$coefs_se)
  mat[1:n.coefs, 3:4] = x$CI
  mat[n.coefs + 1, ] <- NA
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 1:2] <- cbind(x$derived, x$derived_se)
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 3:4] <- x$CI_derived
  rownames(mat) <- c(names(x$coefs), "---", names(x$derived))
  colnames(mat) <- c("Estimate", "Std. Error", colnames(x$CI))
  detfn <- c(hn = "Halfnormal", hhn = "Hazard halfnormal", hr = "Hazard rate", th = "Threshold",
             lth = "Log-link threshold", ss = "Signal strength")[x$detfn]
  infotypes <- c(bearing = "Bearings", dist = "Distances", ss = "Signal strengths",
                 toa = "Times of arrival", mrds = "Exact locations")[x$infotypes]
  
  pars_ext_links = x$pars_ext_links
  ss_opts = x$ss_opts
  cutoff = ss_opts$cutoff
  ss_link = ss_opts$link
  CI_method = ifelse(x$is_boot, "Percentile", "Wald")
  
  cat("Detection function:", detfn, "\n")
  cat("Number of sessions:", x$n.sessions, "\n")
  cat("Information types: ")
  cat(infotypes, sep = ", ")
  cat("\n")
  cat("Confidence interval method:", CI_method, "\n")
  
  
  
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



