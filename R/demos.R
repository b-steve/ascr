demo_fit = function(data_name, fit = TRUE){
  #data_name could be one of follows (also could be found in ./data/):
  #bearing_dist_hn, bearing_hn, dist_hn, ihd, ihd_ext, mul_ses, mul_ses_ext,
  #simple_hhn, simple_hn, simple_hr, ss, ss_toa, toa
  dat = get(data_name)
  if(fit){
    model_output = with(dat, fit.ascr.tmb(capt = capt_input, traps = traps, mask = mask, sv = sv_input, fix = fix_input,
                                          detfn = detfn, local = local_input, bounds = bounds_input, cue.rates = cue.rates,
                                          sound.speed = sound.speed, ss.opts = ss.opts, survey.length = survey.length,
                                          par.extend = par.extend))
    return(list(data = dat, model_output = model_output))
  } else {
    return(list(data = dat, model_output = NULL))
  }
}
