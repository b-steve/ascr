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
