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
  dat$gr_skip = gradient_free
  dat$sv_link = sv_link
  
  dat_model = do.call('read.ascr', dat)
  
  if(fit){
    model_output = fit.ascr(dat_model)
    return(list(data_input = dat, data_for_model = dat_model, model_output = model_output))
  } else {
    return(list(data_input = dat, data_for_model = dat_model, model_output = NULL))
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

