#' Show demo
#'
#' @param data_name 
#' @param fit 
#'
#' @return
#' @export
demo_fit = function(data_name, fit = TRUE){
  dat = get(data_name)
  
  dat_model = do.call('fit.data', dat)
  
  if(fit){
    model_output = fit.ascr(dat_model)
    return(list(data_input = dat, data_for_model = dat_model, model_output = model_output))
  } else {
    return(list(data = dat, data_for_model = dat_model, model_output = NULL))
  }
}



#' Title
#'
#' @param table_print 
#'
#' @return
#' @export
#'
#' @examples
show_demo_options = function(table_print = TRUE){
  output = c('bearing_dist_hn', 'bearing_hn', 'dist_hn', 'ihd', 'ihd_ext', 'mul_ses', 'mul_ses_ext',
             'simple_hhn', 'hhn_cue', 'simple_hr', 'ss', 'ss_toa', 
             'ind_bearing_dist', 'ind_toa_hhn', 'ind_ss', 'ind_ss_log', 'ind_ss_sp')
  if(table_print){
    descriptions = matrix(c('hn' ,'bearing & dist' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'bearing' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'dist' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'NULL' ,'D' ,'FALSE' ,'1',
                            'hn' ,'NULL' ,'sigma & D' ,'FALSE' ,'1',
                            'hn' ,'bearing & dist' ,'NULL' ,'FALSE' ,'2',
                            'hn' ,'bearing & dist' ,'g0 & sigma' ,'FALSE' ,'2',
                            'hhn' ,'NULL' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'cue_rates' ,'NULL' ,'FALSE' ,'1',
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
    print(descriptions)
  }
  
  invisible(output)
}

