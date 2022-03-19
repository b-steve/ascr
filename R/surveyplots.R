
#' Title
#'
#' @param fit 
#' @param new_covariates data.frame; contains any covariates will be used for all extended parameters (if not be skipped)
#' @param param_extend_skip character; skip extended parameter, for skipped extended parameters, use its intercept as the value for this parameter
#' @param xlim 
#' @param ylim 
#' @param main 
#' @param xlab 
#' @param ylab 
#' @param col 
#' @param add 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
show_detfn_tmb <- function(fit, new_covariates = NULL, param_extend_skip = NULL, xlim = NULL, ylim = NULL,
                           main = NULL, xlab = NULL, ylab = NULL, col = NULL, add = FALSE, ...){

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




#' Title
#'
#' @param fit 
#' @param session 
#' @param xlim 
#' @param ylim 
#' @param zlim 
#' @param scale 
#' @param plot.contours 
#' @param add 
#' @param D_cov a list; 2 possible components: "session", a data frame with only one row which contains
#' session level extended covariates, and "mask", a data frame contains location (coordinators by columns x and y)
#' and the location related covariates
#' @param control_convert a list; used for control the convert from D_cov$mask to model fitting mask
#'
#' @return
#' @export
#'
#' @examples
show_Dsurf <- function(fit, session = 1, new_data = NULL, D_cov = NULL, xlim = NULL, ylim = NULL,
                        x_pixels = 50, y_pixels = 50, zlim = NULL, scale = 1, plot_contours = TRUE,
                        add = FALSE, control_convert = NULL){
  
  if(is.null(new_data)){

    #if new_data is not provided, but xlim and ylim provided, we use xlim and ylim to build mask instead of 
    #extracting mask from fit
    if(any(!is.null(xlim), !is.null(ylim))){
      if(!all(!is.null(xlim), !is.null(ylim))) stop('please provide both xlim and ylim.')
      x = seq(from = xlim[1], to = xlim[2], length.out = x_pixels)
      y = seq(from = ylim[1], to = ylim[2], length.out = y_pixels)
      mask = data.frame(x = rep(x, each = y_pixels), y = rep(y, x_pixels))
      is_new_mask = TRUE
    } else {
      mask = get_mask(fit)[[session]]
      is_new_mask = FALSE
    }
    
  } else {

    #if new_data is provided, then xlim and ylim will just do what they are supposed to do,
    #to trim the plot instead of building "mask"
    stopifnot(any(is(new_data, 'data.frame'), is(new_data, 'matrix')))
    stopifnot(all(c('x', 'y') %in% colnames(new_data)))
    mask = new_data[, c('x', 'y')]
    is_new_mask = TRUE
  }
  
  traps = get_trap(fit)[[session]]
  
  if('D' %in% get_par_extend_name(fit)){


    #firstly deal with all scenarios that there may be any new covariate provided 
    if(!is.null(D_cov) | !is.null(new_data)){

      if(!is.null(D_cov)){
        stopifnot(is(D_cov, 'list'))
        stopifnot(any(c('session', 'mask') %in% names(D_cov)))
      }

      #considering the possibility that user may only want to change part of covariates and remains other as the same
      #as the model, for example, 3 covariates related to D, 'weather' (session related), 'noise' and 'forest_type'(loc related)
      #and user only want to change weather, or to change noise, or to change noise and forest_type, and etc.

      #so we create 2 data frame, one for all new covariates, and one for old covariates. We must make sure these two data frame
      #contains the same mask points, and then replace the covariates in the "old" data frame with the same column in the "new",
      #if the same covariate appears in the "new" data frame.
      if(!is.nul(new_data)){
        #in new_data, user could include any location related covariates directly, and it also contains x and y,
        #so we could directly use it instead of mask
        new_covariates = as.data.frame(new_data)
      } else {
        new_covariates = as.data.frame(mask)
      }

      old_covariates = as.data.frame(mask)
      
      #build the new_covariates based on all information we could have
      if(!is.null(D_cov$mask)){
        if(is.null(control_convert)){
          control_convert = vector('list', 2)
          names(control_convert) = c('mask', 'loc_cov')
          control_convert$mask = list(mask)
          control_convert$loc_cov = D_cov$mask
        } else {
          control_convert$mask = list(mask)
          control_convert$loc_cov = D_cov$mask
        }
        
        cov_mask = do.call('location_cov_to_mask', control_convert)
        new_covariates = cbind(new_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      }
      

      if(!is.null(D_cov$session)){
        stopifnot(nrow(D_cov$session) == 1)
        for(i in colnames(D_cov$session)) new_covariates[[i]] = D_cov$session[,i]
        
      }
      
      #build the old_covariates
      ##for location related covariates, if user assigns new mask points by either new_data or xlim & ylim
      ##we need to do interpolation again, otherwise use the argument from the fit
      data_in_model = get_par_extend_data(fit)
      if(is_new_mask){
        old_loc_cov = get_loc_cov(fit)

        if(!is.null(old_loc_cov)){
          if(is.null(control_convert)){
            control_convert = vector('list', 2)
            names(control_convert) = c('mask', 'loc_cov')
            control_convert$mask = list(mask)
            control_convert$loc_cov = old_loc_cov
          } else {
            control_convert$mask = list(mask)
            control_convert$loc_cov = old_loc_cov
          }

          cov_mask = do.call('location_cov_to_mask', control_convert)
          old_covariates = cbind(old_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
        }


      } else if(!is.null(data_in_model$mask)){
        tem = data_in_model$mask
        col_names_tem = colnames(tem)
        if('session' %in% col_names_tem){
          tem = tem[tem$session == session, , drop = FALSE]
          tem = tem[, -which(col_names_tem == 'session'), drop = FALSE]
        }

        if('mask' %in% col_names_tem){
          tem = tem[order(tem$mask),]
          tem = tem[, -which(col_names_tem == 'mask'), drop = FALSE]
        }
        
        old_covariates = cbind(old_covariates, tem)
      }

      ##for session related covariates, it is simple, just take them from the input argument of fit
      
      if(!is.null(data_in_model$session)){
        tem = data_in_model$session[which(data_in_model$session == session), , drop = FALSE]
        tem = tem[, -which(colnames(tem) == 'session'), drop = FALSE]
        for(i in colnames(tem)) old_covariates[[i]] = tem[,i]
      }


      #update the old_covariates by the new_covariates

      for(i in colnames(old_covariates)){
        if(all(i != 'x', i != 'y', i %in% colnames(new_covariates))){
          old_covariates[[i]] = new_covariates[[i]]
        }
      }

      
      par_info = get_data_param(fit)
      par_info = subset(par_info, par == 'D')
      gam = get_gam(fit, 'D')
      values_link = as.vector(coef(fit, types = 'linked', pars = 'D'))
      tem = try({D.mask = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask, values_link, old_covariates)})
      if(is(tem, 'try-error')) stop('Please make sure all covariates needed for "D" are provided.')
      D.mask = unlink.fun(link = par_info$link, value = D.mask)
    } else {
      if(is_new_mask){
        #in this scenario, we need to extract the location related covariates (not mask, but the data provided by user)
        #and do the "location_cov_to_mask" again.
        old_covariates = as.data.frame(mask)
        old_loc_cov = get_loc_cov(fit)

        if(!is.null(old_loc_cov)){
          if(is.null(control_convert)){
            control_convert = vector('list', 2)
            names(control_convert) = c('mask', 'loc_cov')
            control_convert$mask = list(mask)
            control_convert$loc_cov = old_loc_cov
          } else {
            control_convert$mask = list(mask)
            control_convert$loc_cov = old_loc_cov
          }

          cov_mask = do.call('location_cov_to_mask', control_convert)
          old_covariates = cbind(old_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
        }

        data_in_model = get_par_extend_data(fit)
        if(!is.null(data_in_model$session)){
          tem = data_in_model$session[which(data_in_model$session == session), , drop = FALSE]
          tem = tem[, -which(colnames(tem) == 'session'), drop = FALSE]
          for(i in colnames(tem)) old_covariates[[i]] = tem[,i]
        }

        par_info = get_data_param(fit)
        par_info = subset(par_info, par == 'D')
        gam = get_gam(fit, 'D')
        values_link = as.vector(coef(fit, types = 'linked', pars = 'D'))
        tem = try({D.mask = get_extended_par_value(gam, par_info$n_col_full, par_info$n_col_mask, values_link, old_covariates)})
        if(is(tem, 'try-error')) stop('Please make sure all covariates needed for "D" are provided.')
        D.mask = unlink.fun(link = par_info$link, value = D.mask)


      } else {
        D.mask <- fit$D.mask[[session]]
      }
      
    }

  } else {
    #when D is not extended, it is literally a constant, just take the 1st estimated D from session 1
    #if D varies between sessions, it is extended, so there is no problem to just take session1
    D.mask <- rep(fit$D.mask[[1]][1], nrow(mask))

  }
  
  
  
  if (is.null(xlim)){
    xlim <- range(mask[, 'x'])
  }
  if (is.null(ylim)){
    ylim <- range(mask[, 'y'])
  }
  

  mask.keep <- xlim[1] <= mask[, 1] & xlim[2] >= mask[, 1] &
    ylim[1] <= mask[, 2] & ylim[2] >= mask[, 2]
  mask <- mask[mask.keep, ]
  unique.x <- sort(unique(mask[, 1]))
  unique.y <- sort(unique(mask[, 2]))
  z <- squarify(mask, D.mask[mask.keep])
  
  #not sure what "show.cv" does
  #if(!show.cv){
    z <- scale*z
  #}
  
  if (is.null(zlim)){
    zlim <- c(0, max(z, na.rm = TRUE))
  }
  z[z > zlim[2]] <- zlim[2]
  levels <- pretty(zlim, 10)
  if (!add){
    plot(mask, type = "n", asp = 1, xlab = "", ylab = "")
  }
  fields::image.plot(x = unique.x, y = unique.y, z = z, zlim = zlim, col = viridis::viridis(100), add = TRUE)
  points(traps, col = "black", pch = 4, lwd = 2)
  if (plot_contours){
    contour(x = unique.x, y = unique.y, z = z, levels = levels,
            drawlabels = TRUE, add = TRUE)
  }
}






