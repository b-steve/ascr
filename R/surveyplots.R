
#' Title
#'
#' @param fit 
#' @param new_covariates data.frame; contains any covariates will be used for all extended parameters (if not be skipped)
#' @param param_extend_skip character; skip extended parameter, for skipped extended parameters,
#' use its intercept as the value for this parameter
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
show_detfn <- function(fit, new_covariates = NULL, param_extend_skip = NULL, xlim = NULL, ylim = NULL,
                       main = NULL, xlab = NULL, ylab = NULL, col = NULL, add = FALSE, ...){

  det_fn = get_detfn(fit)
  
  if(is.null(xlab)) xlab = "Distance"
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
#' @param show.cv 
#' @param new_data 
#' @param D_cov 
#' @param xlim 
#' @param ylim 
#' @param x_pixels 
#' @param y_pixels 
#' @param zlim 
#' @param scale 
#' @param plot_contours 
#' @param add 
#' @param control_convert_loc2mask
#' @param arg_col 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
show_Dsurf <- function(fit, session = NULL, show.cv = FALSE, new_data = NULL, D_cov = NULL, xlim = NULL, ylim = NULL,
                        x_pixels = 50, y_pixels = 50, zlim = NULL, scale = 1, plot_contours = TRUE,
                        add = FALSE, control_convert_loc2mask= NULL, arg_col = list(n = 100), ...){
  extra_args = list(...)
  
  pred = predict_with_location(fit, session_select = ifelse(is.null(session), 1, session), 
                               new_data = new_data, D_cov = D_cov, xlim = xlim, ylim = ylim,
                                x_pixels = x_pixels, y_pixels = y_pixels, se_fit = show.cv,
                               control_convert_loc2mask= control_convert_loc2mask)
  
  mask = as.matrix(pred[, c('x', 'y')])
  if(!show.cv){
    D.mask = pred[,'est']
  } else {
    D.mask = pred[,'std'] / pred[,'est']
  }
  
  if (is.null(xlim)){
    xlim <- range(mask[, 'x'])
  }
  if (is.null(ylim)){
    ylim <- range(mask[, 'y'])
  }
  

  unique.x <- sort(unique(mask[, 1]))
  unique.y <- sort(unique(mask[, 2]))
  z <- squarify(mask, D.mask)
  

  if(!show.cv){
    z <- scale*z
  }
  
  if (is.null(zlim)){
    zlim <- c(0, max(z, na.rm = TRUE))
  }
  z[z > zlim[2]] <- zlim[2]
  levels <- pretty(zlim, 10)
  
  col_fn = viridis::viridis
  
  if (!add){
    #plot(mask, type = "n", asp = 1, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
    fields::image.plot(x = unique.x, y = unique.y, z = z, zlim = zlim, col = do.call("col_fn", arg_col), 
                       asp = 1, xlim = xlim, ylim = ylim, xlab = "", ylab = "")
  } else {
    fields::image.plot(x = unique.x, y = unique.y, z = z, zlim = zlim, col = do.call("col_fn", arg_col), add = TRUE)
  }
  
  
  if(!is.null(session)){
    traps = get_trap(fit)[[session]]
    trap_plot = extra_args$trap_plot
    if(is.null(trap_plot$col)) trap_col = 1
    if(is.null(trap_plot$pch)) trap_pch = 4
    if(is.null(trap_plot$lwd)) trap_lwd = 2
    points(traps, col = trap_col, pch = trap_pch, lwd = trap_lwd)
  }
  
  
  
  if (plot_contours){
    contour(x = unique.x, y = unique.y, z = z, levels = levels,
            drawlabels = TRUE, add = TRUE)
  }
  invisible(pred)
}




#' Title
#'
#' @param dat 
#' @param session 
#' @param types a character, "survey", "capt" or "covariates"
#' @param control a list,
#' @param anime 
#' @param ask 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.ascr_data <- function(dat, types = NULL, session = NULL, anime = FALSE, control = NULL,
                           ask = TRUE, ...){
  if(is.null(types)) stop('please specify the argument "types".')
  n.sessions = length(dat$traps)
  extra_args = list(...)
  
  ################################################################################################
  if(types == 'survey'){
    if(is.null(session)){
      session = 1
      cat("Argument 'session' is missing, the 1st session is plotted by default.\n")
    }
    
    stopifnot(session <= n.sessions)
    
    masks = get_mask_from_data(dat)[[session]]
    
    if(is.null(control$pch_mask)) pch_mask = "." else pch_mask = control$pch_mask
    if(is.null(control$pch_trap)) pch_trap = 16 else pch_trap = control$pch_trap
    if(is.null(control$col_mask)) col_mask = 8 else col_mask = control$col_mask
    if(is.null(control$col_trap)) col_trap = "red" else col_trap = control$col_trap
    if(is.null(control$cex_mask)) cex_mask = 2 else cex_mask = control$cex_mask
    if(is.null(control$cex_trap)) cex_trap = 1 else cex_trap = control$cex_trap
    
    plot(masks, pch = pch_mask, cex = cex_mask, asp = 1, col = col_mask, ...)
    points(get_trap_from_data(dat)[[session]], pch = pch_trap, col = col_trap, cex = cex_trap)
  }
  
  
  ################################################################################################
  if(types == 'capt'){
    #when session is not assigned, plot all capture history by default
    if(is.null(session)) session = 0
    
    stopifnot(session <= n.sessions)
    
    
    #in this "types", argument cex (cex_det) is used to control the traps' symbol
    #and the cex of the circle used to show the activated traps (cex_capt) should be always cex_det + 2
    #in order to cover the symbol of detectors
    if(is.null(extra_args$cex)){
      cex_det = 3
      cex_capt = 5
    } else {
      cex_det = extra_args$cex
      if(cex_det > 7){
        warning("The size of symbol of traps might be too big to show other information properly.")
      }
      cex_capt = extra_args$cex + 2
    } 
    
    if(is.null(control$fps)) fps = 2 else fps = control$fps
    
    
    #get capture history from data
    #be careful that, in this function, the animal_ID and ID will be
    #converted to natural successive numbers
    capt = get_capt_for_plot(dat)
    c_names = colnames(capt)
    if("bearing" %in% c_names) is.bearing = TRUE else is.bearing = FALSE
    if("dist" %in% c_names) is.dist = TRUE else is.dist = FALSE
    if("ss" %in% c_names) is.ss = TRUE else is.ss = FALSE
    if("animal_ID" %in% c_names) animal.model = TRUE else animal.model = FALSE
    
    if(is.ss & is.dist){
      cat("Extra information of 'signal strengh' and 'distance' are both probvided, 'distance' will be discarded.\n")
      is.dist = FALSE
      capt = capt[, -which(c_names == "dist"),drop = FALSE]
    }
    
    
    #when animal ID or call ID equals to zero, it means all detection.
    #currently, animation could only be applied for single session.
    if(any(!is.null(control$animal_ID), !is.null(control$ID), anime) & session == 0){
      session = 1
      cat("Argument 'session' is missing, the 1st session is plotted by default.\n")
    }
    
    if(animal.model){
      if(is.null(control$animal_ID)){
        a_id = 0
        #we cannot plot with call ID but without animal_ID
        if(!is.null(control$ID)) stop("Please provide information about animal_ID.")
      }  else {
        a_id = control$animal_ID
        if(anime){
          a_id = 0
          cat("Assigned animal_ID will be ignored as animated capture history will be generated.\n")
        }
      }
    }
    
    if(is.null(control$ID)){
      c_id = 0
    } else {
      c_id = control$ID
      if(anime){
        c_id = 0
        cat("Assigned ID will be ignored as animated capture history will be generated.\n")
      }
    } 
    

    t_list = get_trap_from_data(dat)
    m_list = get_mask_from_data(dat)
    ##################################################################################
    if(anime == TRUE){
      s = session
      traps = as.data.frame(t_list[[s]])
      colnames(traps) = c('trap_x', 'trap_y')
      traps$trap = seq(nrow(traps))
      
      masks = m_list[[s]]
      buffer = attr(masks, 'buffer')
      masks = as.data.frame(masks)
      
      capt_session = subset(capt, capt$session == s)
      #for dist model, re-scale ss and dist for plotting purpose 
      capt_session = dist_ss_rescale(capt_session, buffer, is.ss, is.dist, is.bearing)
      
      if(animal.model){
        capt_session$keys = paste(capt_session$animale_ID, capt_session$ID, sep = '-')
      } else {
        capt_session$keys = capt_session$ID
      }
      
      capt_session = merge(capt_session, traps, by = "trap", all.x = TRUE)
      capt_session = sort.data(capt_session, 'data.full')
      
      #plot basic info
      base_plot = ggplot(data = masks, mapping = aes(x = x, y = y)) + 
        coord_quickmap(xlim = range(masks$x), ylim = range(masks$y)) +
        theme_void()
      trap_plot = base_plot + 
        geom_point(data = traps, mapping = aes(x = trap_x, y = trap_y, shape = 4), size = cex_det) + 
        scale_shape_identity()
      
      
      #if is.bearing, plot arrow, otherwise plot circle or just point
      if(is.bearing){
        if(is.ss){
          #if is.dist or is.ss, use length of the arrow to indicate dist or -1*ss,
          #the information in column "dist" and "ss" are already been rescaled
          #in the previous step by function dist_ss_rescale(), so here we use them directly

          capt_session$arrow_len = capt_session$ss
          trap_plot = trap_plot + labs(caption = "arrow length shows re-scaled -1*ss.")
        } else if(is.dist){
          capt_session$arrow_len = capt_session$dist
          trap_plot = trap_plot + labs(caption = "arrow length shows re-scaled dist.")
        } else {
          capt_session$arrow_len = rep(0.5*buffer, nrow(capt_session))
        }
        
        capt_plot = trap_plot +
          geom_point(data = capt_session, mapping = aes(x = trap_x, y = trap_y, group = keys),
                     size = cex_capt, colour = "red") +
          geom_segment(data = capt_session, mapping = aes(x = trap_x, y = trap_y,
                                                          xend = trap_x + sin(bearing) * arrow_len,
                                                          yend = trap_y + cos(bearing) * arrow_len,
                                                          group = keys),
                       arrow = arrow(length = unit(0.02, "npc")), colour = "red")
        
 
        #end of bearing plot
      } else {
        #if not bearing, use circle if dist or ss is available, use literally point otherwise
        capt_plot = trap_plot +
          geom_point(data = capt_session, mapping = aes(x = trap_x, y = trap_y, group = keys),
                     size = cex_capt, colour = "red")
        

        #if is.dist or is.ss, use radius of circle to indicate scaled -1*dist or ss, since the values
        #have been re scaled, use them directly as radius
        if(is.ss){
          capt_plot = capt_plot +
            geom_point(data = capt_session, mapping = aes(x = trap_x, y = trap_y, group = keys),
                       size = capt_session$ss, colour = "red", shape = 1) +
            labs(caption = "radius shows re-scaled ss.")
        } else if(is.dist){
          capt_plot = capt_plot +
            geom_point(data = capt_session, mapping = aes(x = trap_x, y = trap_y, group = keys),
                       size = capt_session$dist, colour = "red", shape = 1) +
            labs(caption = "radius shows re-scaled dist.")
        }
        #end of non-bearing plot
      }
      
      
      anime_plot = capt_plot +
        gganimate::transition_time(keys)
      
      if(animal.model){
        anime_plot = anime_plot + 
          ggtitle(paste0('session ', s,', animal_ID - call_ID: {frame_time}'),
                  subtitle = 'Frame {frame} of {nframes}')
      } else {
        anime_plot = anime_plot + 
          ggtitle(paste0('session ', s, ', call_ID: {frame_time}'),
                  subtitle = 'Frame {frame} of {nframes}')
      }
      
      
      n_keys = length(unique(capt_session$keys))
      animate(anime_plot, nframes = n_keys, fps = fps)
      
      #end of animation plotting 
      
    } else {
      
      #start of statics plotting
      
      if(ask){
        ## Setting par(ask).
        ask.save <- par("ask")
        par(ask = TRUE)
        ## Making sure par is restored on function exit.
        on.exit(par(ask = ask.save))
      }
      
      
      #if session is 0, plot all sessions, create a vector for a sessions for plot
      if(session == 0) session = seq(n.sessions)
      
      for(s in session){
        traps = as.data.frame(t_list[[s]])
        masks = m_list[[s]]
        buffer = attr(masks, 'buffer')
        masks = as.data.frame(masks)
        capt_session = subset(capt, capt$session == s)
        #for dist model, re-scale ss and dist for plotting purpose 
        capt_session = dist_ss_rescale(capt_session, buffer, is.ss, is.dist, is.bearing)
        
        if(animal.model){
          if(a_id != 0){
            capt_session = subset(capt_session, animal_ID %in% a_id)
            if(c_id != 0){
              capt_session = subset(capt_session, ID %in% c_id)
            }
          }
          keys = paste(capt_session$animal_ID, capt_session$ID)
        } else {
          if(c_id != 0){
            capt_session = subset(capt_session, ID %in% c_id)
          }
          keys = capt_session$ID
        }
        
        u_keys = unique(keys)

        #plot basic info
        base_plot = ggplot(data = masks, mapping = aes(x = x, y = y)) + 
          coord_quickmap(xlim = range(masks$x), ylim = range(masks$y)) +
          theme_void()
        trap_plot = base_plot + 
          geom_point(data = traps, mapping = aes(x = x, y = y, shape = 4), size = cex_det) + 
          scale_shape_identity()
        
        
        if(is.bearing){
          #if bearing, use arrow
          #bearing: x1 = x0 + sin(bearing) * length, y1 = y0 + cos(bearing) * length
          
          if(is.ss){
            #if is.dist or is.ss, use length of the arrow to indicate dist or -1*ss,
            #the information in column "dist" and "ss" are already been rescaled
            #in the previous step by function dist_ss_rescale(), so here we use them directly
            len = capt_session$ss
            trap_plot = trap_plot + labs(caption = "arrow length shows re-scaled -1*ss")
          } else if(is.dist){
            len = capt_session$dist
            trap_plot = trap_plot + labs(caption = "arrow length shows re-scaled dist")
          } else {
            len = rep(0.5*buffer, nrow(capt_session))
          }
          
          for(k in u_keys){
            i_k = which(keys == k)
            one_call = capt_session[i_k,,drop = FALSE]
            actived_traps = traps[one_call$trap,]
            len_1call = len[i_k]
            if(animal.model){
              sub_title = paste0("session: ",s, ", animal ID: ", one_call$animal_ID[1],
                                 ", call ID: ", one_call$ID[1])
            } else {
              sub_title = paste0("session: ",s, ", call ID: ", one_call$ID[1])
            }

            plot_one_call = trap_plot + 
              geom_point(data = actived_traps, mapping = aes(x = x, y = y), size = cex_capt, colour = "red") +
              geom_segment(data = actived_traps, mapping = aes(x = x, y = y,
                                                               xend = x + sin(one_call$bearing) * len_1call,
                                                               yend = y + cos(one_call$bearing) * len_1call),
                           arrow = arrow(length = unit(0.02, "npc")), colour = "red") +
              labs(subtitle = sub_title)
            
            print(plot_one_call)
            
          }
          
          
          #end of if(is.bearing)
        } else {
          #if not bearing, use circle if dist or ss is available, use literally point otherwise
          
          #if is.dist or is.ss, use radius of circle to indicate scaled -1*dist or ss, since the values
          #have been re scaled, use them directly as radius
          if(is.ss){
            rad = capt_session$ss
            trap_plot = trap_plot + labs(caption = "radius shows re-scaled ss.")
          } else if(is.dist){
            rad = capt_session$dist
            trap_plot = trap_plot + labs(caption = "radius shows re-scaled dist.")
          }
            
          
          for(k in u_keys){
            i_k = which(keys == k)
            one_call = capt_session[i_k,,drop = FALSE]
            actived_traps = traps[one_call$trap,]
            
            if(animal.model){
              sub_title = paste0("session: ",s, ", animal ID: ", one_call$animal_ID[1],
                                 ", call ID: ", one_call$ID[1])
            } else {
              sub_title = paste0("session: ",s, ", call ID: ", one_call$ID[1])
            }
            
            plot_one_call = trap_plot + 
              geom_point(data = actived_traps, mapping = aes(x = x, y = y), size = cex_capt, colour = "red") +
              labs(subtitle = sub_title)
            
            #if is.dist or is.ss, add circles
            if(is.dist | is.ss){

              rad_1call = rad[i_k]
              plot_one_call = plot_one_call +
                geom_point(data = actived_traps, mapping = aes(x = x, y = y), size = rad_1call, shape = 1,
                           colour = "red")
            } 
            print(plot_one_call)
            
          }
          #end of !is.bearing plot
        }
        
       #end of session s 
      }
      
      #end of non-anime plot
    }
  
  #end of types "capt"
  }
  

  ################################################################################################

  if(types == 'covariates'){
    if(is.null(session)){
      session = 1
      cat("Argument 'session' is missing, the 1st session is plotted by default.\n")
    }
    
    masks = as.data.frame(get_mask_from_data(dat)[[session]])
    masks$mask = seq(nrow(masks))
    masks_mat = as.matrix(masks[, c('x', 'y'), drop = FALSE])
    D_cov_for_model = dat$par.extend$data$mask
    if(is.null(D_cov_for_model)){
      stop('There is no location related covariates, nothing to plot.')
    }
    
    D_cov_input = dat$arg.input$loc_cov
    
    
    D_cov_for_model = D_cov_for_model[D_cov_for_model$session == session,
                                      -which(colnames(D_cov_for_model) == 'session')]
    
    D_cov_for_model = merge(D_cov_for_model, masks, by = 'mask')
    D_cov_for_model = D_cov_for_model[, -which(colnames(D_cov_for_model) == 'mask')]
    
    cov_list = colnames(D_cov_for_model)
    cov_list = cov_list[-which(cov_list == 'x' | cov_list == 'y')]
    
    unique.x <- sort(unique(masks[, 'x']))
    unique.y <- sort(unique(masks[, 'y']))
    xlim <- range(masks[, 'x'])
    ylim <- range(masks[, 'y'])
    col_fn = viridis::viridis
    if(is.null(control$arg_col)) arg_col = list(n = 100, option = "D")
    if(is.null(control$plot_contours)) plot_contours = TRUE
    
    if(ask){
      ## Setting par(ask).
      ask.save <- par("ask")
      par(ask = TRUE)
      ## Making sure par is restored on function exit.
      on.exit(par(ask = ask.save))
    }
    
    
    for(i in cov_list){
      if(is(D_cov_for_model[[i]], 'numeric')){
        
        z <- squarify(masks_mat, D_cov_for_model[[i]])
        zlim <- c(0, max(z, na.rm = TRUE))
        fields::image.plot(x = unique.x, y = unique.y, z = z, zlim = zlim, col = do.call("col_fn", arg_col), 
                           asp = 1, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
                           main = paste0("Plot of covariate ", i, ", for session ", session))
        
        if(plot_contours){
          levels <- pretty(zlim, 10)
          contour(x = unique.x, y = unique.y, z = z, levels = levels, drawlabels = TRUE, add = TRUE)
        }
      } else {
        
        v = D_cov_for_model[[i]]
        v_num = as.numeric(as.factor(v))
        
        match_table = data.frame(code = v_num, name = v)
        match_table = match_table[!duplicated(match_table),,drop = FALSE]
        
        z <- squarify(masks_mat, v_num)
        zlim <- range(match_table$code)
        
        col = do.call("col_fn", arg_col)
        idx = round(scale_convert(match_table$code, seq(arg_col$n)),0)
        
        col_leng = col[idx]
        
        image(x = unique.x, y = unique.y, z = z, zlim = zlim, col = col, 
              asp = 1, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
              main = paste0("Plot of covariate ", i, ", for session ", session))
        legend(x = "topright", legend = match_table$name, fill = col_leng)
        
      }
    }
    #end of types == 'covariates'
  }
  
    
}


#' Title
#'
#' @param fit 
#' @param session 
#' @param types 
#' @param control 
#' @param main 
#' @param ... 
#'  
#'
#' @return
#' @export
#'
#' @examples
plot.ascr_tmb = function(fit, session = 1, types = NULL, control = NULL, main = NULL, ...){
  if(types %in% c('survey', 'capt', 'D_covariates')){
    plot.ascr_data(fit$args, session = session, types = types, control = control, ...)
  }
}
