#' Show demo
#'
#' @param data_name 
#' @param fit 
#' @param gradient_free 
#' @param sv_link 
#'
#' @return
#' @export
demo_fit = function(data_name, fit = TRUE, gradient_free = FALSE, sv_link = NULL){
  dat = get(data_name)
  
  par_extend_model = dat$par_extend_model
  sv = dat$sv
  fix = dat$fix
  bounds = dat$bounds
  detfn = dat$detfn
  ss_opts = dat$ss.opts
  
  dat[c('par_extend_model', 'sv', 'fix', 'bounds', 'detfn', 'ss.opts')] = NULL
  
  

  
  dat_model = list()
  
  dat_model$dat = do.call('read.ascr', dat)
  dat_model$par_extend_model = par_extend_model
  dat_model$detfn = detfn
  dat_model$sv = sv
  dat_model$fix = fix
  dat_model$bounds = bounds  
  dat_model$ss.opts = ss_opts
  dat_model$gr_skip = gradient_free
  dat_model$sv_link = sv_link
  
  
  
  if(fit){
    model_output = do.call('fit.ascr', dat_model)
    return(list(read_input = dat, read_output = dat_model$dat, dat_model = dat_model, fit = model_output))
  } else {
    return(list(read_input = dat, read_output = dat_model$dat, dat_model = dat_model, fit = NULL))
  }
}



#' Title
#'
#' @param table_return 
#'
#' @return
#' @export
#'
#' @examples
show_demo_options = function(table_return = TRUE){
  output = c('bearing_dist_hn', 'bearing_hn', 'dist_hn', 'ihd', 'ihd_ext', 'mul_ses', 'mul_ses_ext',
             'simple_hhn', 'hhn_cue', 'simple_hr', 'ss', 'ss_toa', 
             'ind_bearing_dist', 'ind_toa_hhn', 'ind_ss', 'ind_ss_log', 'ind_ss_sp')
  if(table_return){
    descriptions = matrix(c('hn' ,'bearing & dist' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'bearing' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'dist' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'NULL' ,'D' ,'FALSE' ,'1',
                            'hn' ,'NULL' ,'sigma & D' ,'FALSE' ,'1',
                            'hn' ,'bearing & dist' ,'NULL' ,'FALSE' ,'2',
                            'hn' ,'bearing & dist' ,'g0 & sigma' ,'FALSE' ,'2',
                            'hhn' ,'NULL' ,'NULL' ,'FALSE' ,'1',
                            'hhn' ,'cue_rates' ,'NULL' ,'FALSE' ,'1',
                            'hr' ,'NULL' ,'NULL' ,'FALSE' ,'1',
                            'ss' ,'ss' ,'NULL' ,'FALSE' ,'1',
                            'ss' ,'ss & toa' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'bearing & dist' ,'alpha & D' ,'TRUE' ,'3',
                            'hhn' ,'toa' ,'D' ,'TRUE' ,'2',
                            'ss' ,'ss' ,'b0.ss & D' ,'TRUE' ,'3',
                            'ss_log', 'ss', 'NULL', 'TRUE', '3',
                            'ss_spherical', 'ss', 'NULL', 'TRUE', '3'),
                          nrow = 17, byrow = T)
    
    descriptions = cbind(1:17, descriptions)
    colnames(descriptions) = c('index', 'det_fn', 'extra_info', 'extended_par', 'individual_id', 'n_sessions')
    rownames(descriptions) = output
    descriptions = as.data.frame(descriptions)

    return(descriptions)
  } else {
    return(output)
  }

}

