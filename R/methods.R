#' Extract coefficients from the output of ascr_tmb model
#'
#' @param object a fitted model from fit.ascr().
#' @param types a character vector, accept any subset from ('all', 'fitted', 'linked', 'derived'), default is 'linked'. In details:
#'              linked - There is a link function attached to each parameter, either 'log', 'logit' or 'identical'. The 'linked' 
#'                       estimation is the estimation before back-transformation according to the link functions. For any extended
#'                       parameters, the estimations of the coefficients of any covariates will be always under link function of
#'                       identity. For example, if we have D ~ x1, the estimations of D_int and D_x1 are 1 and 0.5, if the value
#'                       of 'x1' is not provided in the argument of 'new_covariates', then no matter what the 'types' is, these
#'                       estimations will be shown all the time. If x1 = 3, then when 'linked', the function will return 1 + 0.5 * 3
#'                       as the estimation of log(D); when 'fitted', the function will return exp(1 + 0.5 * 3) as the estimation of D.
#'              fitted - The back transformed estimations from the 'linked' estimations and their link function.
#'              derived - The estimations of 'esa' for each session.
#'                       
#' @param pars a character vector containing any parameter names.
#' @param new_covariates a data frame containing the values of covariates of any extended parameter.
#'
#' @return a named numeric vector
#' @export
#'
coef.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, ...){
  
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
  extend_par_covariates = get_par_extend_covariate(object)
  
  for(i in pars){
    par_info = subset(df_param, par == i)
    values_link = param_values_og[[i]]
    #this is an indicator about whether we calculate \beta \times X, which means, for the extended parameter
    #whether we evaluate the link function
    do_BX = FALSE
    
    if(!i %in% name_extend){
      names(values_link) = paste(i, "link", sep = "_")
    } else {
      names(values_link) = gsub(" _ ", "\\.", names(values_link))
      names(values_link) = paste(names(values_link), "link", sep = "_")
      if(all(extend_par_covariates[[i]] %in% colnames(new_covariates))){
        do_BX = TRUE
      }
    }
    
    if(is.linked){
      if(!do_BX){
        output = c(output, values_link)
      } else {
        gam = get_gam(object, i)
        values_link_evaluated = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                       values_link, new_covariates)
        names(values_link_evaluated) = rep(paste(i, "link", sep = "_"), length(values_link_evaluated))
        output = c(output, values_link_evaluated)
      }

    }
    
    if(is.fitted){
      link = par_info$link
      if(i %in% name_extend){
        if(!do_BX){
          #if there is no new covariates assigned and the parameter is extended, then regards all of the
          #relevant beta as identity linked parameters
          values_fitted = values_link
          names(values_fitted) = gsub("_link", "", names(values_fitted))
          link = 'identity'
        } else {
          gam = get_gam(object, i)
          values_fitted = get_extended_par_value(gam, par_info$n_col_full,
                                                 par_info$n_col_mask, values_link, new_covariates)

          names(values_fitted) = rep(i, length(values_fitted))
          
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


#' Title
#'
#' @param object an object generated by the bootstrap function "boot.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param correct_bias logical. if TRUE, apply a bias correction method to the bootstrap results; otherwise,
#'                     return the original estimation directly. The bias correct method is
#'                     hat(coef) - (mean(boot_results) - hat(coef)).
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
coef.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, correct_bias = FALSE, ...){
  
  if(!correct_bias){
    output = coef.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates)
    return(output)
  }

  
  ########################################################################################################
  #some preparation
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
  
  #get foundation information needed for this method
  name_extend = get_par_extend_name(object)
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  

  #######################################################################################################
  

  df_param = get_data_param(object)
  res = get_boot_res(object, pars)
  coefs = coef.ascr_tmb(object, types = 'linked', pars)
  bias = get_bias(res, coefs)
  coefs_linked_bias_corrected = coefs - bias
  
  linked_name = names(coefs)
  original_name = ori_name(linked_name)
  extend_par_covariates = get_par_extend_covariate(object)
  
  output = vector('list', length(types))
  names(output) = types
  
  for(i in types){
    
    if(i == 'linked'){
      if(is.null(new_covariates)){
        output[[i]] = coefs_linked_bias_corrected
      } else {
        #if new_covariates is provided, then we need to go each parameters one by one
        for(j in pars){
          values_link = coefs_linked_bias_corrected[original_name == j]
          par_info = subset(df_param, par == j)
          
          if(j %in% name_extend){
            if(!all(extend_par_covariates[[j]] %in% colnames(new_covariates))){
              values_link_evaluated = values_link
            } else {
              #if(nrow(new_covariates) != 1) stop('Argument "new_covariates" can only accept 1 row.')
              gam = get_gam(object, j)
              values_link_evaluated = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                     values_link, new_covariates)
              
              names(values_link_evaluated) = rep(paste(j, 'link', sep = '_'), length(values_link_evaluated))
            }
          } else {
            values_link_evaluated = values_link
          }
          
          output[[i]] = c(output[[i]], values_link_evaluated)
        }
        
      }
      
    } else if (i == 'derived'){
      coefs_esa = coef.ascr_tmb(object, types = 'derived')
      res_esa = get_boot_res_esa(object)
      bias_esa = get_bias(res_esa, coefs_esa)
      output[[i]] = coefs_esa - bias_esa
      
    } else {
      #when types is 'fitted', use the 'coefs_linked_bias_corrected' and new_covariates to
      #calculate the estimations.

      for(j in pars){
        values_link = coefs_linked_bias_corrected[original_name == j]
        
        par_info = subset(df_param, par == j)
        link = par_info$link
        
        if(j %in% name_extend){
          if(!all(extend_par_covariates[[j]] %in% colnames(new_covariates))){
            #if there is no new covariates assigned and the parameter is extended, then regards all of the
            #relevant beta as identity linked parameters
            values_fitted = values_link
            names(values_fitted) = gsub("_link", "", names(values_fitted))
            link = 'identity'
          } else {
            #if(nrow(new_covariates) != 1) stop('Argument "new_covariates" can only accept 1 row.')
            
            gam = get_gam(object, j)
            values_fitted = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask,
                                                   values_link, new_covariates)
            
            names(values_fitted) = rep(j, length(values_fitted))
          }
        } else {
          values_fitted = values_link
          names(values_fitted) = j
        }
        
        values_fitted = unlink.fun(link = link, value = values_fitted)

        output[[i]] = c(output[[i]], values_fitted)
      }
      #end of if(i == 'fitted')
    }
    
    #end of for i in types
  }

  names(output) = NULL
  output = do.call('c', output)

  
  class(output) = "coef_ascr_tmb"
  return(output)
  
}




##################################################################################################################

#' S3 method for output of coef.ascr_tmb and coef.ascr_boot
#'
#' @param x an object generated by the function "coef.ascr_tmb()" or "coef.ascr_boot()".
#'
#' @export
#'
print.coef_ascr_tmb = function(x){
  output = as.matrix(x)
  colnames(output) = "Est"
  print(output)
}

##################################################################################################################

#' Extract variance covariance matrix of the estimated parameters from ascr models
#'
#' @param object a fitted model from fit.ascr().
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param show_fixed_par a logical value. To control whether to include the fixed parameters in the covariance matrix.
#'                       It is TRUE by default.
#' @param ... 
#'
#' @return a list with matrices as its elements if multiple 'types'. a matrix if only one 'types'.
#' @export
vcov.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, show_fixed_par = TRUE, ...){
  
  #deal with default setting for 'types'
  tem = types_pars_sol(types, pars, new_covariates)
  types = tem$types
  pars = tem$pars
  
  
  # 'og' below means original output from the object
  param_values_og = get_coef(object)
  cov_og = object$vcov
  name_dim_og = colnames(cov_og)
  name_dim = gsub(" _ ", "\\.", name_dim_og)
  which.derived = which(substr(name_dim, 1, 3) == 'esa')
  name_dim = name_dim[-which.derived]
  name_dim_og = name_dim_og[-which.derived]
  cov_derived = cov_og[which.derived, which.derived, drop = FALSE]
  cov_linked = cov_og[-which.derived, -which.derived, drop = FALSE]
  dimnames(cov_linked) = list(name_dim, name_dim)
  name_extend = get_par_extend_name(object)
  extend_par_covariates = get_par_extend_covariate(object)
  
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  #browser()
  name_extned_covariate_provided = ext_par_in_new_df(name_extend, new_covariates, extend_par_covariates)
  
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
  #name_og is naturally in the right order because it is extracted from the TMB var-cov matrix output
  
  #obtain the indices for the assigned "pars"
  index_par = which(name_og %in% pars)
  name_og = name_og[index_par]
  name_dim = name_dim[index_par]
  cov_linked = cov_linked[index_par, index_par, drop = FALSE]
  link_funs = link_funs[index_par]
  param_values = param_values[index_par]
  fixed_par = get_fixed_par_name(object)
  
  
  output = vector('list', length(types))
  names(output) = types
  #browser()
  for(type in types){
    if(type == "linked"){
      if(is.null(new_covariates)){
        output[[type]] = cov_linked
        tem = paste(name_dim, 'link', sep = '_')
        dimnames(output[[type]]) = list(tem, tem)
      } else {
        gam.output = get_gam(object)
        #browser()
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates,
                                               name_og = name_og, name_extend = name_extend,
                                               par_ext_cov_provided = name_extned_covariate_provided, df_param = df_param,
                                               gam.output = gam.output, back_trans = FALSE)
        
        new_name = character(0)
        for(p in pars){
          if(p %in% name_extend){
            if(p %in% name_extned_covariate_provided){
              tem = rep(paste0(p, "_link"), nrow(new_covariates))
            } else {
              tem = paste(name_dim[which(name_og == p)], "link", sep = "_")
            }
          } else {
            tem = paste0(p, "_link") 
          }
          new_name = c(new_name, tem)
        }
        
        dimnames(output[[type]]) = list(new_name, new_name)
        
      }
      
      
      if(show_fixed_par){
        output[[type]] = vcov_fixed_par_add(output[[type]], fixed_par, type)
      }

    } else if(type == 'derived'){
      output[[type]] = cov_derived
    } else if(type == 'fitted'){
      #if 'new_covariates' is not provided, just keep the extended covariates as parameters with identity link function
      if(is.null(new_covariates)){
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, link_funs = link_funs, back_trans = TRUE)
        dimnames(output[[type]]) = list(name_dim, name_dim)
      } else {
        gam.output = get_gam(object)
        #browser()
        output[[type]] = delta_method_ascr_tmb(cov_linked, param_values, new_covariates = new_covariates,
                                               name_og = name_og, name_extend = name_extend,
                                               par_ext_cov_provided = name_extned_covariate_provided, df_param = df_param,
                                               gam.output = gam.output, back_trans = TRUE)
        
        new_name = character(0)
        for(p in pars){
          if(p %in% name_extend){
            if(p %in% name_extned_covariate_provided){
              tem = rep(p, nrow(new_covariates))
            } else {
              tem = name_dim[which(name_og == p)]
            }
          } else {
            tem = p 
          }
          new_name = c(new_name, tem)
        }
        
        dimnames(output[[type]]) = list(new_name, new_name)
      }
      
      if(show_fixed_par){
        output[[type]] = vcov_fixed_par_add(output[[type]], fixed_par, type)
      }
      
    }
  }
  
  #if only one type be assigned, change the output from a list to a matrix only
  if(length(types) == 1) output = output[[types]]
  
  return(output)
  
}


#' Title
#'
#' @param object an object generated by the bootstrap function "boot.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param from_boot a logical value, TRUE by default, to control whether the covariance matrix is calculated
#'                  from the bootstrap results. If FALSE, the method "vcov.ascr_tmb()" will be called.
#' @param show_fixed_par a logical value, the same as "vcov.ascr_tmb()".
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
vcov.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, from_boot = TRUE, show_fixed_par = TRUE, ...){
  
  if(!from_boot){
    output = vcov.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates, show_fixed_par = show_fixed_par)
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
  
  #browser()


  output = vector('list', length(types))
  names(output) = types
  for(i in types){
    if(i == 'linked'){
      res_linked = res_transform(res, new_covariates, pars, object, back_trans = FALSE)
      output[[i]] = var_from_res(res_linked, fixed_par)
      
      if(show_fixed_par){
        output[[i]] = vcov_fixed_par_add(output[[i]], fixed_par, i)
      }
      
    } else if(i == 'derived'){
      #browser()
      res_esa = get_boot_res_esa(object)
      output[[i]] = var_from_res(res_esa)
    } else {
      res_fitted = res_transform(res, new_covariates, pars, object, back_trans = TRUE)
      output[[i]] = var_from_res(res_fitted, fixed_par)
      
      if(show_fixed_par){
        output[[i]] = vcov_fixed_par_add(output[[i]], fixed_par, i)
      }
      
    }
  }
  
  #when only one types, directly output it without a list structure
  if(length(types) == 1) output = output[[types]]
 
  return(output)

    
}



##################################################################################################################

#' Extract standard errors of the estimated parameters from ascr models
#'
#' @param object a fitted model from "fit.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param show_fixed_par a logical value, the same as "vcov.ascr_tmb()".
#' @param ... 
#'
#' @return a named numeric vector
#' @export
stdEr.ascr_tmb = function(object, types = NULL, pars = NULL, new_covariates = NULL, show_fixed_par = TRUE, ...){
  
  output_vcov = vcov.ascr_tmb(object = object, types = types, pars = pars, new_covariates = new_covariates,
                              show_fixed_par = show_fixed_par)
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
#' @param object an object generated by the bootstrap function "boot.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param from_boot a logical value, the same as "vcov.ascr_boot()".
#' @param show_fixed_par a logical value, the same as "vcov.ascr_tmb()".
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
stdEr.ascr_boot = function(object, types = NULL, pars = NULL, new_covariates = NULL, from_boot = TRUE, show_fixed_par = TRUE, ...){

  output_vcov = vcov.ascr_boot(object = object, types = types, pars = pars, new_covariates = new_covariates,
                               from_boot = from_boot, show_fixed_par = show_fixed_par)
  
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
#' @param object a fitted model from "fit.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
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
#' @param object an object generated by the bootstrap function "boot.ascr()".
#' @param types a character vector, the same as "coef.ascr_tmb()".
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param pars a character vector, the same as "coef.ascr_tmb()".
#' @param new_covariates a data frame, the same as "coef.ascr_tmb()".
#' @param correct_bias a logical value indicates whether to apply a bias correction method to
#'                     the confidence interval. Default is FALSE, and the naive percentiles from
#'                     the bootstrap results will used as confidence interval; if TRUE, the 
#'                     bootstrap results will be corrected by beta_boot = 2 * beta_est - beta_boot,
#'                     where beta_boot is the bootstrap results values and beta_est is the
#'                     estimations from "coef.ascr_tmb()", and the percentiles from these corrected
#'                     bootstrap results are used as the confidence interval.
#' @param from_boot a logical value indicates whether to use the bootstrap results to construct 
#'                  the confidence matrix. Default is TRUE; if FALSE, "confint.ascr_tmb()" will
#'                  be called.
#' @param ...
#' 
#'
#' @return
#' @export
#'
#' @examples
#' 
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
      res_linked = res_transform(res, new_covariates, pars, object, back_trans = FALSE)
      output[[i]] = boot_res_to_CI(res_linked, level)
    } else if(i == 'derived'){
      res_esa = get_boot_res_esa(object)
      coef_esa = coef.ascr_tmb(object, types = 'derived')
      res_esa = res_mod_for_CI(res_esa, coef_esa, correct_bias)
      
      output[[i]] = boot_res_to_CI(res_esa, level)
    } else {
      res_fitted = res_transform(res, new_covariates, pars, object, back_trans = TRUE)
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
  
  n_fixed = length(get_fixed_par_name(object))
  
  return(k * (length(coef(object)) - n_fixed) - 2 * object$loglik)
}



#' Title
#'
#' @param fit a fitted model from "fit.ascr()".
#' @param type a character vector, could be either "response" or "link". The "response" correspondence to
#'             the type of "fitted" in the function "coef.ascr_tmb()" and the "link" correspondence to
#'             the type of "linked" in the function "coef.ascr_tmb()".
#' @param newdata a data frame, the same as "coef.ascr_tmb()". If there is any extended parameters, and
#'                the corresponding covariates is not provided here, the estimation for such parameters
#'                will not be shown.
#' @param se.fit a logical value indicates whether to show standard error.
#' @param confidence a logical value indicates whether to show confidence interval.
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param realnames a character vector containing any parameter names.
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
predict.ascr_tmb = function(fit, type = 'response', newdata = NULL, se.fit = TRUE, confidence = TRUE,
                            level = 0.95, realnames = NULL, ...){
  stopifnot(all(type %in% c('link', 'response'), length(type) == 1))
  
  if(type == 'response'){
    type = 'fitted'
  } else if(type == 'link'){
    type = 'linked'
  } else {
    stop('invalid "type", it could be either "response" or "link".')
  }
  
  name_og = get_param_og(fit)
  
  #check which parameter will be displayed
  if(is.null(realnames)){
    realnames = name_og
  } else {
    if(any(!realnames %in% name_og)) stop("Argument 'realnames' only accept parameters' name in this model.")
  }
  
  #make sure pars are in the right order
  realnames = fulllist.par.generator()[fulllist.par.generator() %in% realnames]
  
  par_fixed_name = get_fixed_par_name(fit)
  par_ext_name = get_par_extend_name(fit)
  par_ext_vars = get_par_extend_covariate(fit)
  newdata_vars = colnames(newdata)
  
  output = vector('list', length(realnames))
  names(output) = realnames
  
  for(i in realnames){
    if(i %in% par_ext_name){
      if(!all(par_ext_vars[[i]] %in% newdata_vars)){
        n_row = 0
        output[[i]] = paste0('Estimates not provided, because ', i, ' depends on ', paste(par_ext_vars[[i]], collapse = ', '),
                             ', which are not all provided in the argument "newdata".')
      } else {
        #if the condition above is FALSE, then newdata must not be NULL
        n_row = nrow(newdata)
      }
      
    } else {
      n_row = 1
    }
    
    if(n_row != 0){
      output[[i]] = data.frame(Estimate = numeric(n_row))

      output[[i]]$Estimate = as.vector(coef(fit, types = type, new_covariates = newdata, pars = i))
      
      if(se.fit){
        if(i %in% par_fixed_name){
          output[[i]]$StdError = 0
        } else {
          output[[i]]$StdError = as.vector(stdEr(fit, types = type, new_covariates = newdata, pars = i, show_fixed_par = FALSE))
        }
        
      } 
      if(confidence){
        if(i %in% par_fixed_name){
          output[[i]]$Lower = output[[i]]$Estimate
          output[[i]]$Upper = output[[i]]$Estimate
        } else {
          output[[i]]$Lower = numeric(n_row)
          output[[i]]$Upper = numeric(n_row)
          output[[i]][, c('Lower', 'Upper')] = confint(fit, types = type, new_covariates = newdata, pars = i)
        }

      }
    }

  }
  
  class(output) = 'predict_ascr_tmb'
  return(output)
  
}



#' Title
#'
#' @param fit an object generated by the bootstrap function "boot.ascr()".
#' @param type a character vector, the same as "predict.ascr_tmb()".
#' @param newdata a data frame, the same as "predict.ascr_tmb()".
#' @param se.fit a logical value indicates whether to show standard error.
#' @param confidence a logical value indicates whether to show confidence interval.
#' @param level a numeric value indicates the confident level, default is 0.95.
#' @param realnames a character vector, the same as "predict.ascr_tmb()".
#' @param correct_bias a logical value indicates whether to apply the bias correction method to the
#'                     point estimations and the confidence intervals. The method is described in the
#'                     functions "coef.ascr_boot()" and "confint.ascr_boot()".
#' @param from_boot a logical value indicates whether to extract standard error and confidence interval
#'                  from the bootstrap results.
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
predict.ascr_boot = function(fit, type = 'response', newdata = NULL, se.fit = TRUE, confidence = TRUE,
                            level = 0.95, realnames = NULL, correct_bias = FALSE, from_boot = TRUE, ...){
  if(!from_boot){
    if(correct_bias) message(paste0("The argument 'correct_bias' will be ignored as the confidence",
                                    " interval is not from bootstrap results."))
    correct_bias = FALSE
  }
  
  
  stopifnot(all(type %in% c('link', 'response'), length(type) == 1))
  
  if(type == 'response'){
    type = 'fitted'
  } else {
    type = 'linked'
  }
  
  name_og = get_param_og(fit)
  
  #check which parameter will be displayed
  if(is.null(realnames)){
    realnames = name_og
  } else {
    if(any(!realnames %in% name_og)) stop("Argument 'realnames' only accept parameters' name in this model.")
  }
  
  #make sure pars are in the right order
  realnames = fulllist.par.generator()[fulllist.par.generator() %in% realnames]
  
  par_fixed_name = get_fixed_par_name(fit)
  par_ext_name = get_par_extend_name(fit)
  par_ext_vars = get_par_extend_covariate(fit)
  newdata_vars = colnames(newdata)
  
  output = vector('list', length(realnames))
  names(output) = realnames
  
  for(i in realnames){
    if(i %in% par_ext_name){
      if(!all(par_ext_vars[[i]] %in% newdata_vars)){
        n_row = 0
        output[[i]] = paste0('Estimates not provided, because ', i, ' depends on ', paste(par_ext_vars[[i]], collapse = ', '),
                             ', which are not all provided in the argument "newdata".')
      } else {
        #if the condition above is FALSE, then newdata must not be NULL
        n_row = nrow(newdata)
      }
      
    } else {
      n_row = 1
    }
    
    if(n_row != 0){
      output[[i]] = data.frame(Estimate = numeric(n_row))

      output[[i]]$Estimate = as.vector(coef(fit, types = type, new_covariates = newdata, pars = i, correct_bias = correct_bias))
      
      if(se.fit){
        
        if(i %in% par_fixed_name){
          output[[i]]$StdError = 0
        } else {
          output[[i]]$StdError = as.vector(stdEr(fit, types = type, new_covariates = newdata, pars = i, from_boot = from_boot, show_fixed_par = FALSE))
        }
        
      }

      if(confidence){
        if(i %in% par_fixed_name){
          output[[i]]$Lower = output[[i]]$Estimate
          output[[i]]$Upper = output[[i]]$Estimate
        } else {
          output[[i]]$Lower = numeric(n_row)
          output[[i]]$Upper = numeric(n_row)
          output[[i]][, c('Lower', 'Upper')] = confint(fit, types = type, new_covariates = newdata, pars = i, correct_bias = correct_bias,
                                                       from_boot = from_boot)
        }
      }
    }
    
  }
  
  class(output) = 'predict_ascr_tmb'
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
print.predict_ascr_tmb = function(x, ...){
  for(i in names(x)){
    cat(paste0(i, ": \n"))
    if(is(x[[i]], 'data.frame')){
      print(x[[i]], row.names = FALSE)
      cat("\n")
    } else {
      cat(paste0(x[[i]], "\n"))
      cat("\n")
    }
  }
}




#' Title
#'
#' @param object an object generated from the model fitting function "fit.ascr_tmb()" or
#'               the bootstrap process "boot.ascr()".
#' @param derived_print a logical value indicates whether to show the estimations of "esa".
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.ascr_tmb = function(object, derived_print = FALSE, ...){
  coefs = coef(object, types = 'fitted')
  derived = coef(object, types = 'derived')
  coefs_se = stdEr(object, types = 'fitted', show_fixed_par = FALSE)
  derived_se = stdEr(object, types = 'derived', show_fixed_par = FALSE)
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
                ss_opts = ss_opts, derived_print = derived_print)
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
  
  derived_print = x$derived_print
  
  n.coefs <- length(x$coefs)
  n.derived <- length(x$derived)
  mat <- matrix(0, nrow = n.coefs + n.derived + 1, ncol = 4)
  mat[1:n.coefs, 1] <- as.vector(x$coefs)
  for(i in 1:n.coefs){
    name_coef = names(x$coefs)[i]
    index_coefs_se = which(names(x$coefs_se) == name_coef)
    if(length(index_coefs_se) != 0){
      mat[i, 2] <- as.vector(x$coefs_se[index_coefs_se])
    } else {
      mat[i, 2] <- 0
    }
    
  }
  
  mat[1:n.coefs, 3:4] = x$CI
  mat[n.coefs + 1, ] <- NA
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 1:2] <- cbind(x$derived, x$derived_se)
  mat[(n.coefs + 2):(n.coefs + n.derived + 1), 3:4] <- x$CI_derived
  rownames(mat) <- c(names(x$coefs), "---", names(x$derived))
  colnames(mat) <- c("Estimate", "Std. Error", colnames(x$CI))
  
  if(!derived_print){
    mat = mat[1:n.coefs,,drop = FALSE]
  }
  
  
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

  if(length(infotypes) != 0){
    cat(infotypes, sep = ", ")
  } else {
    cat("NULL")
  }
  
  
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



