#' Show demo
#'
#' @param data_name 
#' @param fit 
#'
#' @return
#' @export
demo_fit = function(data_name, fit = TRUE, dev = FALSE){
  #data_name could be one of follows (also could be found in ./data/):
  #bearing_dist_hn, bearing_hn, dist_hn, ihd, ihd_ext, mul_ses, mul_ses_ext,
  #simple_hhn, simple_hn, simple_hr, ss, ss_toa, toa
  
  #for individual identification included model, the names of demo data are:
  #ind_bearing_dist, ind_toa_hhn, ind_ss
  
  dat = get(data_name)
  dat$dev = dev
  if(fit){
    model_output = with(dat, fit_og(capt = capt_input, traps = traps, mask = mask, sv = sv_input, fix = fix_input,
                                      detfn = detfn, local = local_input, bounds = bounds_input, cue.rates = cue.rates,
                                      sound.speed = sound.speed, ss.opts = ss.opts, survey.length = survey.length,
                                      par.extend = par.extend, dev = dev))
    return(list(data = dat, model_output = model_output))
  } else {
    return(list(data = dat, model_output = NULL))
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
             'simple_hhn', 'simple_hn_cue', 'simple_hr', 'ss', 'ss_toa', 'toa',
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
                            'hn' ,'toa' ,'NULL' ,'FALSE' ,'1',
                            'hn' ,'bearing & dist' ,'alpha & D' ,'TRUE' ,'3',
                            'hhn' ,'toa' ,'D' ,'TRUE' ,'2',
                            'ss' ,'ss' ,'b0.ss & D' ,'TRUE' ,'3',
                            'ss_log', 'ss', 'NULL', 'TRUE', '3',
                            'ss_spherical', 'ss', 'NULL', 'TRUE', '3'),
                          nrow = 18, byrow = T)
    
    descriptions = cbind(1:18, descriptions)
    colnames(descriptions) = c('index', 'det_fn', 'extra_info', 'extended_par', 'individual_id', 'n_sessions')
    rownames(descriptions) = output
    descriptions = as.data.frame(descriptions)
    print(descriptions)
  }
  
  invisible(output)
}

