#dx and each component of det_par should be vectors with the same length
#of course, if any one of them is only a scaler, there will be no problem
det_prob = function(det_fn, det_par, dx, ss.link = NULL){
  
  for(i in names(det_par)) assign(i, det_par[[i]])
  
  if(det_fn == 'hn'){
    out = g0 * exp((-1 * dx^2) / (2 * sigma^2))
  } else if(det_fn == 'hhn'){
    out = 1 - exp((-1) * lambda0 * exp(dx^2 * (-0.5) / sigma^2))
  } else if(det_fn == 'hr'){
    out = g0 * (1 - exp((-1) * (dx/sigma)^(-z)))
  } else if(det_fn == 'th'){
    out = 0.5 - 0.5 * erf(dx / scale - shape)
  } else if(det_fn == 'lth'){
    out = 0.5 - 0.5 * erf(shape.1 - exp(shape.2 - scale * dx))
  } else if(det_fn == 'ss'){
    #ss is a little different, if the simulated ss > cut off, the detection probability is 1
    #and 0 otherwise, so we return E(SS) instead of detected probability
    if(ss.link == 'identity'){
      out = b0.ss - b1.ss * dx
    } else if(ss.link == 'log'){
      out = exp(b0.ss - b1.ss * dx)
    } else if(ss.link == 'spherical'){
      if(dx > 1){
        out = b0.ss - 20 * log10(dx) - b1.ss * (dx - 1)
      } else {
        out = b0.ss
      }
    }
  }
  
  return(out)
  
}


#only two arguments needs some explanation
#param_extend_skip: if the user wants to skip any extended parameter, could input it as a char vector here
#new_covariates: contains any covariates will be used for all extended parameters (if not be skipped), a data.frame
#for skipped extended parameters, use its intercept as the value for this parameter

show_detfn_tmb <- function(fit, new_covariates = NULL, param_extend_skip = NULL, xlim = NULL, ylim = NULL,
                           main = NULL, xlab = NULL, ylab = NULL, col = NULL, add = FALSE, ...){
  #library(mgcv)
  #source('get_funs.r', local = TRUE)
  #source('support_functions.r', local = TRUE)
  det_fn = get_detfn(fit)
  
  if(is.null(xlab)) xlab = "Distance (m)"
  if(is.null(ylab)){
    if(det_fn == 'ss') ylab = "E(ss)"
    if(det_fn != 'ss') ylab = "Detection probability"
  }
  
  #get some essential component from output
  param_name = get_param_og(fit)
  param_name = param_name[-which(param_name == 'D')]
  data_param = get_data_param(fit)
  par_extend_name = get_par_extend_name(fit)
  param_values = get_coef(fit)
  ss_link = get_ss_link(fit)
  
  #deal with xlim and the values will be used as x-axis (distances)
  if (is.null(xlim)){
    buffer <- get_buffer(fit)
    #get_buffer return a vector with length of n.sessions
    #and if the buffer for on session is not provided, corresponding element is zero
    if (!any(buffer == 0)){
      x.max <- max(buffer)
    } else {
      x.max <- max(get_dist_theta(fit)$dx)
    }
    xlim <- c(0, x.max)
  }
  dists <- seq(xlim[1], xlim[2], length.out = 1000)
  
  #deal with skipped extended parameters and new_covariates
  
  if(is.null(new_covariates) & !is.null(param_extend_skip)){
    if(!all(param_extend_skip == par_extend_name[-which(par_extend_name == 'D')])){
      warning('No new data of covariates provided, all extended parameters will be skipped.')
    }
  }
  
  #if no new_covariates assigned, skip all extended parameters
  if(is.null(new_covariates) & is.null(param_extend_skip)){
    if(length(par_extend_name) > 0){
      print('No new data of covariates provided, only the intercept will be used for all extended parameters.')
      param_extend_skip = par_extend_name
    }
  }
  
  if(!is.null(new_covariates)){
    #standardize new_covariates first
    new_covariates = fit$scale.covs(new_covariates)
    print(paste0("Covariates for detect function's parameters are provided. Please make sure that, ",
                 "for each extended parameter which are not included in the argument 'param_extend_skip',",
                 "all covariates are provided"))
  }
  
  #set the values of every each detection function parameters
  det_param_input = vector('list', length(param_name))
  names(det_param_input) = param_name
  
  for(i in param_name){
    #get the basic information for this parameter
    tem = data_param[which(data_param$par == i),]
    link = tem$link
    n_col_full = tem$n_col_full
    n_col_mask = tem$n_col_mask
    
    values = param_values[[i]]
    
    if((i %in% par_extend_name) & (!i %in% param_extend_skip)){
      gam = get_gam(fit, i)
      det_param_input[[i]] = get_extended_par_value(gam, n_col_full, n_col_mask, values, new_covariates)
    } else {
      det_param_input[[i]] = values[1]
    }
    
    #use 'link' to back-transform the parameter's value
    
    det_param_input[[i]] = unlink.fun(link = link, value = det_param_input[[i]])
    
    names(det_param_input[[i]]) = NULL
  }
  
  #be careful, each component in det_param_input may have different length
  #but as.data.frame() will automatically solve this
  tem = as.data.frame(det_param_input)
  n_lines = nrow(tem)
  for(i in param_name) det_param_input[[i]] = tem[, i]
  
  probs = vector('list', n_lines)
  tem_det_par = vector('list', length(param_name))
  names(tem_det_par) = param_name
  
  for(i in 1:n_lines){
    #extract i'th value for each parameter from det_param_input
    for(j in param_name) tem_det_par[[j]] = det_param_input[[j]][i]
    probs[[i]] = det_prob(det_fn, tem_det_par, dists, ss_link)
  }
  
  #based on all probs values, determine ylim
  if(is.null(ylim)){
    if(det_fn != 'ss'){
      ylim = c(0, 1)
    } else {
      tem = do.call('c', probs)
      ylim = range(tem)
    }
  }
  
  if(is.null(col)){
    col = 1:n_lines
  } else {
    col = rep(col, length = n_lines)
  }
  
  if (!add){
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, xaxs = 'i', yaxs = "r")
    axis(1)
    axis(2)
    box()
    if(det_fn != 'ss') abline(h = c(0, 1), col = "lightgrey")
    title(main = main, xlab = xlab, ylab = ylab)
  }
  
  for(i in 1:n_lines){
    lines(dists, probs[[i]], col = col[i], ...)
  }
  
}
