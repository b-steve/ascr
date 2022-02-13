
test_that("joint bearing/dist fitting -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_print = F)[1])$model_output
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.04365338912, 1.68535571449, 3.93414147381, 1.61475318813, 7.72525085354)
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.0308861068845, 0.1338424162815, 0.0989981712497, 0.1020762055936)
  pars_link_names = c("g0_link", "sigma_link", "kappa_link", "alpha_link", "D_link")
  pars_names = c('g0', 'sigma', 'kappa', 'alpha', 'D')
  
  
  #test estimations as a whole, and one by one by selecting "par"
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  for(i in 1:5){
    relative.error = as.vector(abs((coef(fit, types = 'linked', pars = pars_names[i]) - pars_link_est_values[i])/pars_link_est_values[i]))
    expect_equal(relative.error, 0, tolerance = 1e-4)
  }
  
  #test std errors as a whole, and one by one
  pars_link_names_unfixed = pars_link_names[-1]
  pars_names_unfixed = pars_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  for(i in 1:4){
    relative.error = as.vector(abs((stdEr(fit, types = 'linked', pars = pars_names_unfixed[i]) - pars_link_std_values[i])/pars_link_std_values[i]))
    expect_equal(relative.error, 0, tolerance = 1e-4)
  }
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  #test confidence interval
  conf_95 = matrix(c(36.04365338912, 1.62482005737, 3.67181515830, 1.42072033795, 7.52518516690, 36.04365338912, 1.74589137161,
                     4.19646778933, 1.80878603831, 7.92531654018), ncol = 2)
  
  conf_95_unfixed = conf_95[-1,, drop = FALSE]
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  for(i in 1:4){
    relative.error = max(abs((confint(fit, types = 'linked', linked = FALSE, level = 0.95, pars = pars_names_unfixed[i]) - conf_95_unfixed[i,,drop = FALSE])/conf_95_unfixed[i,,drop = FALSE]))
    expect_equal(relative.error, 0, tolerance = 1e-4)
  }
  
  
  conf_90 = matrix(c(36.04365338912, 1.63455258956, 3.71399028995, 1.45191568709, 7.55735043654, 36.04365338912, 1.73615883942,
                     4.15429265767, 1.77759068917,  7.89315127054), ncol = 2)
  
  conf_90_unfixed = conf_90[-1,, drop = FALSE]
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  for(i in 1:4){
    relative.error = max(abs((confint(fit, types = 'linked', linked = FALSE, level = 0.9, pars = pars_names_unfixed[i]) - conf_90_unfixed[i,,drop = FALSE])/conf_90_unfixed[i,,drop = FALSE]))
    expect_equal(relative.error, 0, tolerance = 1e-4)
  }
  
  
})