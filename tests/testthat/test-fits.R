
test_that("joint bearing/dist fitting -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[1])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.04365, 1.607376, 3.876604, 1.892325, 7.835287)

  pars_link_names = c("g0_link", "sigma_link", "kappa_link", "alpha_link", "D_link")
  pars_names = c('g0', 'sigma', 'kappa', 'alpha', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)

  #test fitted estimations
  pars_est_values = c(1.000000, 4.989701, 48.260046, 6.634774, 2528.260710)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.03049075, 0.1397999, 0.1059173, 0.101818)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  pars_names_unfixed = pars_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  #test fitted std errors
  pars_std_values = c(0.1521397, 6.7467492, 0.7027375, 257.4224243)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names_unfixed] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.043653, 1.547615, 3.602601, 1.684730, 7.635727, 36.043653, 1.667137,
                     4.150607, 2.099919, 8.034846), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(36.043653, 1.557223, 3.646654, 1.718106, 7.667811, 36.043653, 1.657529,
                     4.106554, 2.066543,  8.002763), ncol = 2)

  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(1.000000, 4.691513, 35.036661, 5.257433, 2023.722030, 1.000000, 5.287890,
                            61.483431, 8.012114, 3032.799390), ncol = 2)
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(1.000000, 4.700248, 36.693560, 5.390997, 2070.876621, 1.000000, 5.296980,
                                   63.472502, 8.165506, 3086.664919), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density model with sigma extended -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[5])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.0436534, 2.0142606, -0.2688631, 5.7910711, 0.1208786)
  
  pars_link_names = c("g0_link", "sigma.(Intercept)_link", "sigma.brandsony_link", "D.(Intercept)_link", "D.noise_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.13479140, 0.06918262, 2.06221583, 0.21047749)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.0436534, 1.7500743, -0.4044586, 1.7492023, -0.2916497, 36.0436534, 2.2784469,
                     -0.1332677, 9.8329399, 0.5334069), ncol = 2)

  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(36.0436534, 1.7925485, -0.3826584, 2.3990279, -0.2253260, 36.0436534, 2.2359727,
                     -0.1550679, 9.1831143,  0.4670833), ncol = 2)

  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  

  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', noise = 7.7)
  pars_names_og = c('g0', 'sigma', 'D')
  pars_names_og_unfixed = pars_names_og[-1]
  
  expected_values = c(1, 5.728178, 830.341102)
  o = coef(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.9222248, 400.2406206)
  o = stdEr(fit, new_covariates = new_data)[pars_names_og_unfixed]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1, 3.92065, 45.88390, 1, 7.535705, 1614.798303),
                           ncol = 2)
  o = confint(fit, new_covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(6.721837, 0.4820195, 5.777096, 7.666577), nrow = 1)
  o = predict(fit, type = 'link', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(830.3411, 400.2406, 45.8839, 1614.798), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("joint bearing/dist model with 2 sessions and g0 & sigma extended -- half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[7])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.5378647, -1.0670451, 1.0895859, 0.1642574, 3.5215606, 1.0289551, 7.8714508)
  
  pars_link_names = c("g0.(Intercept)_link", "g0.weathersunny_link", "sigma.(Intercept)_link", "sigma.brandsony_link",
                      "kappa_link", "alpha_link", "D_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test linked std errors
  pars_link_std_values = c(1.58088441, 1.36630628, 0.07671678, 0.07736679, 0.51800547, 0.21100953, 0.17920083)
 
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(-1.56061182, -3.74495623, 0.93922376, 0.01262129, 2.50628848, 0.61538402, 7.52022367,
                     4.6363412, 1.6108660, 1.2399480, 0.3158935, 4.5368326, 1.4425262, 8.2226780), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(-1.06245877, -3.31441897, 0.96339802, 0.03700037, 2.66951737, 0.68187531, 7.57669171,
                     4.1381881, 1.1803287, 1.2157738, 0.2915144, 4.3736037, 1.3760349, 8.1662100), ncol = 2)

  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  #test fitted confident interval
  pars_names = c("g0.(Intercept)", "g0.weathersunny", "sigma.(Intercept)", "sigma.brandsony",
                 "kappa", "alpha", "D")
  
  conf_95_fitted = matrix(c(-1.56061182, -3.74495623, 0.93922376, 0.01262129, -0.51676398, 1.64091053, 1700.67103015,
                            4.6363412, 1.6108660, 1.2399480, 0.3158935, 68.1911480, 3.9553705, 3542.0610432), ncol = 2)
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(-1.56061182, -3.74495623, 0.93922376, 0.01262129, 12.25934474, 1.85036704, 1844.97990856,
                                   4.6363412, 1.6108660, 1.2399480, 0.3158935, 93.3945155, 4.2313715, 3724.4632673), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', weather = 'sunny')
  pars_names_og = c('g0', 'sigma', 'kappa', 'alpha', 'D')
  
  expected_values = c(0.6155777, 3.5037832, 33.8371920, 2.7981405, 2621.3660367)
  o = coef(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.1210091, 0.2763892, 17.5278506, 0.5904343, 469.7509821)
  o = stdEr(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(0.3784042, 2.9620703, -0.5167640, 1.6409105, 1700.6710302,
                             0.8527513, 4.0454960, 68.1911480, 3.9553705, 3542.0610432),
                           ncol = 2)
  o = confint(fit, new_covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("cue rate included -- hazard half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[9])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.405026, 1.519666, 6.203402)
  
  pars_link_names = c("sigma_link", "lambda0_link", "D_link")
  pars_names = c('sigma', 'lambda0', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(4.075632, 4.570699, 494.428405)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  #test linked std errors
  pars_link_std_values = c(0.07907883, 0.19694115, 0.13699999)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(0.3222962, 0.9001587, 67.7366874)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(1.250034, 1.133669, 5.934887, 1.560017, 1.905664, 6.471917), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(1.274953, 1.195727, 5.978057, 1.535099, 1.843605, 6.428747), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(3.443943, 2.806420, 361.666938, 4.707320, 6.334977, 627.189873), ncol = 2)
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(3.490462, 3.107034, 377.997397, 4.758904, 6.723868, 646.722570), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("Signal strength & toa model", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[12])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.497327, 1.356604, 2.228033, -6.398831, 7.812649)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "sigma.toa_link", "D_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'sigma.toa', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(8.977683e+01, 3.882984e+00, 9.281592e+00, 1.663501e-03, 2.471669e+03)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test linked std errors
  pars_link_std_values = c(0.01398757, 0.05076599, 0.04959745, 0.08643598, 0.09945806)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(1.255760e+00, 1.971236e-01, 4.603433e-01, 1.437863e-04, 2.458274e+02)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.469912, 1.257105, 2.130824, -6.568242, 7.617715,
                     4.524742, 1.456104, 2.325242, -6.229420, 8.007583), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.474319, 1.273101, 2.146453, -6.541006, 7.649055,
                     4.520334, 1.440107, 2.309614, -6.256657, 7.976243), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(8.731559e+01, 3.496629e+00, 8.379336e+00, 1.381685e-03, 1.989856e+03,
                            9.223808e+01, 4.269339e+00, 1.018385e+01, 1.945317e-03, 2.953482e+03), ncol = 2)
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(8.734902e+01, 3.515228e+00, 8.421802e+00, 1.404263e-03, 2.033909e+03,
                                   9.227213e+01, 4.289214e+00, 1.022916e+01, 1.970595e-03, 3.003649e+03), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density & toa model with individual identity -- hazard half normal", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[14])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(0.7713772, 1.8352614, -6.6387914, 10.6122325, -0.4882451, 2.1657845)
  
  pars_link_names = c("sigma_link", "lambda0_link", "sigma.toa_link", "D.(Intercept)_link", "D.noise_link", "mu_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test std error
  pars_link_std_values = c(0.07513875, 0.26995338, 0.10891509, 5.11553494, 0.47807457, 0.14269433)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(0.6241079, 1.3061625, -6.8522610, 0.5859682, -1.4252540, 1.8861088,
                     0.9186464, 2.3643603, -6.4253217, 20.6384967, 0.4487639, 2.4454603), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(0.6477849, 1.3912276, -6.8179407, 2.1979263, -1.2746077, 1.9310732,
                     0.8949694, 2.2792952, -6.4596420, 19.0265387, 0.2981176, 2.4004958), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(noise = 7.7)
  pars_names_og = c('sigma', 'lambda0', 'sigma.toa', 'D', 'mu')
  
  expected_values = c(2.162743e+00, 6.266772e+00, 1.308608e-03, 9.464759e+02, 8.721441e+00)
  o = coef(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(1.625058e-01, 1.691736e+00, 1.425272e-04, 1.374890e+03, 1.244500e+00)
  o = stdEr(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1.844237e+00, 2.951030e+00, 1.029260e-03, -1.748258e+03, 6.282266e+00,
                             2.481248e+00, 9.582514e+00, 1.587956e-03, 3.641210e+03, 1.116062e+01),
                           ncol = 2)
  o = confint(fit, new_covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(6.852746, 1.452641, 4.005622, 9.699869), nrow = 1)
  o = predict(fit, type = 'link', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(946.4759, 1374.89, -1748.258, 3641.21), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(946.4759, 54.90596, 16315.47), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, linked = T, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("signal strength model with log link and individual identity", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[16])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(1.676551, -2.532202, 2.287277, 4.529954, 2.162836)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(5.3470820, 0.0794838, 9.8480851, 92.7542590, 8.6957638)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.001742142, 0.015412422, 0.025178517, 0.174155514, 0.063878707)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(0.009315378, 0.001225038, 0.247960173, 16.153665696, 0.555474145)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(1.673136, -2.562410, 2.237928, 4.188615, 2.037636,
                     1.679966, -2.501994, 2.336626, 4.871292, 2.288036), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(1.673685, -2.557553, 2.245862, 4.243493, 2.057765,
                     1.679417, -2.506851, 2.328692, 4.816414, 2.267907), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(5.32882417, 0.07708277, 9.36209207, 61.09365602, 7.60705450,
                            5.36533978, 0.08188483, 10.33407809, 124.41486199, 9.78447314), ncol = 2)
  
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(5.32885530, 0.07711868, 9.37388885, 65.93141820, 7.67245027,
                                   5.36537098, 0.08192147, 10.34626944, 130.48942065, 9.85556188), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("signal strength model with spherical link and individual identity", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[17])$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.627961, -1.679443, 2.305293, 4.234750, 2.247667)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(102.3052278, 0.1864778, 10.0271191, 69.0444275, 9.4656292)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #g0 is fixed, so it has no std error estimation
  pars_link_std_values = c(0.02051216, 0.61022968, 0.02797164, 0.19999997, 0.06543341)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(2.0985008, 0.1137943, 0.2804750, 13.8088834, 0.6193684)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.587758, -2.875471, 2.250470, 3.842757, 2.119420,
                     4.6681639, -0.4834148, 2.3601167, 4.6267429, 2.3759144), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.594221, -2.683182, 2.259284, 3.905779, 2.140039,
                     4.6617003, -0.6757045, 2.3513026, 4.5637209, 2.3552956), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(98.19224173, -0.03655491, 9.47739826, 41.97951330, 8.25168944,
                            106.4182138, 0.4095105, 10.5768399, 96.1093417, 10.6795690), ncol = 2)
  
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(98.27382219, 0.05638956, 9.49219541, 46.65394249, 8.32630792,
                                   106.5020104, 0.6166739, 10.5921879, 102.1807100, 10.7608483), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density model with sigma extended -- half normal - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[5], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(36.0436534, 2.0142606, -0.2688632, 5.7910710, 0.1208786)
  
  pars_link_names = c("g0_link", "sigma.(Intercept)_link", "sigma.brandsony_link", "D.(Intercept)_link", "D.noise_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  pars_link_std_values = c(0.13479106, 0.06918256, 2.06221107, 0.21047694)
  
  #test linked std errors
  pars_link_names_unfixed = pars_link_names[-1]
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names_unfixed] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #expect an error here as g0 is fixed, no estimated std error
  expect_error(stdEr(fit, par = 'g0'))
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(36.0436534, 1.7500750, -0.4044585, 1.7492116, -0.2916486, 36.0436534, 2.2784462,
                     -0.1332678, 9.8329304, 0.5334059), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(36.0436534, 1.7925490, -0.3826583, 2.3990357, -0.2253251, 36.0436534, 2.2359722,
                     -0.1550680, 9.1831064,  0.4670824), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(brand = 'sony', noise = 7.7)
  pars_names_og = c('g0', 'sigma', 'D')
  pars_names_og_unfixed = pars_names_og[-1]
  
  expected_values = c(1, 5.728178, 830.341085)
  o = coef(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(0.9222219, 400.2402285)
  o = stdEr(fit, new_covariates = new_data)[pars_names_og_unfixed]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)

  expected_values = matrix(c(1, 3.920656, 45.884652, 1, 7.535699, 1614.797518),
                           ncol = 2)
  o = confint(fit, new_covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(6.721837, 0.4820191, 5.777097, 7.666577), nrow = 1)
  o = predict(fit, type = 'link', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(830.3411, 400.2402, 45.88465, 1614.798), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("Signal strength & toa model - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[12], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.497327, 1.356604, 2.228033, -6.398831, 7.812649)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "sigma.toa_link", "D_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'sigma.toa', 'D')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(8.977683e+01, 3.882984e+00, 9.281592e+00, 1.663501e-03, 2.471669e+03)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test linked std errors
  pars_link_std_values = c(0.01398754, 0.05076602, 0.04959744, 0.08643598, 0.09945805)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(1.255757e+00, 1.971236e-01, 4.603432e-01, 1.437863e-04, 2.458274e+02)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.469912, 1.257104, 2.130824, -6.568242, 7.617715,
                     4.524742, 1.456104, 2.325242, -6.229420, 8.007583), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.474320, 1.273101, 2.146453, -6.541006, 7.649055,
                     4.520334, 1.440107, 2.309614, -6.256656, 7.976243), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(8.731560e+01, 3.496629e+00, 8.379336e+00, 1.381685e-03, 1.989856e+03,
                            9.223807e+01, 4.269340e+00, 1.018385e+01, 1.945317e-03, 2.953482e+03), ncol = 2)
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(8.734903e+01, 3.515228e+00, 8.421803e+00, 1.404263e-03, 2.033909e+03,
                                   9.227212e+01, 4.289214e+00, 1.022916e+01, 1.970595e-03, 3.003649e+03), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})


test_that("heterogeneous density & toa model with individual identity -- hazard half normal - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[14], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(0.7713772, 1.8352614, -6.6387914, 10.6122319, -0.4882450, 2.1657846)
  
  pars_link_names = c("sigma_link", "lambda0_link", "sigma.toa_link", "D.(Intercept)_link", "D.noise_link", "mu_link")
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test std error
  pars_link_std_values = c(0.07513802, 0.26995232, 0.10891498, 5.11553391, 0.47807445, 0.14269430)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(0.6241094, 1.3061645, -6.8522608, 0.5859697, -1.4252537, 1.8861089,
                     0.9186450, 2.3643582, -6.4253220, 20.6384941, 0.4487637, 2.4454603), ncol = 2)
  
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  
  #test another confidence level
  conf_90 = matrix(c(0.6477861, 1.3912293, -6.8179406, 2.1979274, -1.2746075, 1.9310733,
                     0.8949682, 2.2792934, -6.4596422, 19.0265364, 0.2981175, 2.4004958), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  
  
  ############################################################################################################
  #test for new_data provided
  new_data = data.frame(noise = 7.7)
  pars_names_og = c('sigma', 'lambda0', 'sigma.toa', 'D', 'mu')
  
  expected_values = c(2.162743e+00, 6.266772e+00, 1.308608e-03, 9.464758e+02, 8.721442e+00)
  o = coef(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = c(1.625042e-01, 1.691730e+00, 1.425270e-04, 1.374889e+03, 1.244500e+00)
  o = stdEr(fit, new_covariates = new_data)[pars_names_og]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(1.844240e+00, 2.951043e+00, 1.029260e-03, -1.748257e+03, 6.282266e+00,
                             2.481245e+00, 9.582501e+00, 1.587956e-03, 3.641209e+03, 1.116062e+01),
                           ncol = 2)
  o = confint(fit, new_covariates = new_data)[pars_names_og,]
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  expected_values = matrix(c(6.852745, 1.452641, 4.005622, 9.699869), nrow = 1)
  o = predict(fit, type = 'link', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(946.4758, 1374.889, -1748.257, 3641.209), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  expected_values = matrix(c(946.4758, 54.90597, 16315.46), nrow = 1)
  o = predict(fit, type = 'response', newdata = new_data, linked = T, confidence = T)
  relative.error = max(abs((o - expected_values)/expected_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
})


test_that("signal strength model with spherical link and individual identity - no gradient", {
  #fit the model using 'demo_fit'
  fit = demo_fit(show_demo_options(table_return = F)[17], gradient_free = TRUE)$fit
  
  ##########################################################################################################
  #check coefficients estimations without back transformation
  
  pars_link_est_values = c(4.627962, -1.679417, 2.305293, 4.234735, 2.247669)
  
  pars_link_names = c("b0.ss_link", "b1.ss_link", "sigma.ss_link", "D_link", "mu_link")
  pars_names = c('b0.ss', 'b1.ss', 'sigma.ss', 'D', 'mu')
  
  
  #test linked estimations
  relative.error = max(abs((coef(fit, types = 'linked')[pars_link_names] - pars_link_est_values)/pars_link_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test fitted estimations
  pars_est_values = c(102.3053100, 0.1864827, 10.0271110, 69.0433944, 9.4656476)
  relative.error = max(abs((coef(fit, types = 'fitted')[pars_names] - pars_est_values)/pars_est_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  pars_link_std_values = c(0.02050358, 0.61003840, 0.02797149, 0.20000145, 0.06543333)
  
  
  relative.error = max(abs((stdEr(fit, types = 'linked')[pars_link_names] - pars_link_std_values)/pars_link_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  #test fitted std errors
  pars_std_values = c(2.0976247, 0.1137616, 0.2804732, 13.8087791, 0.6193689)
  relative.error = max(abs((stdEr(fit, types = 'fitted')[pars_names] - pars_std_values)/pars_std_values))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  
  ############################################################################################################
  
  #test confidence interval
  conf_95 = matrix(c(4.587775, -2.875070, 2.250469, 3.842740, 2.119422,
                     4.6681478, -0.4837633, 2.3601156, 4.6267309, 2.3759162), ncol = 2)
  
  #[pars_link_names, ] is use to make sure the order is correct
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.95)[pars_link_names,]
  
  relative.error = max(abs((o - conf_95)/conf_95))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('2.5 %', '97.5 %')))
  
  
  #test another confidence level
  conf_90 = matrix(c(4.594236, -2.682840, 2.259284, 3.905762, 2.140041,
                     4.6616870, -0.6759927, 2.3513015, 4.5637083, 2.3552975), ncol = 2)
  
  
  o = confint(fit, types = 'linked', linked = FALSE, level = 0.9)[pars_link_names,]
  
  relative.error = max(abs((o - conf_90)/conf_90))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  expect_true(all(colnames(o) == c('5 %', '95 %')))
  
  #test fitted confident interval
  conf_95_fitted = matrix(c(98.19404116, -0.03648596, 9.47739355, 41.97868473, 8.25170694,
                            106.4165789, 0.4094514, 10.5768284, 96.1081041, 10.6795883), ncol = 2)
  
  o = confint(fit, types = 'fitted', linked = FALSE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted)/conf_95_fitted))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  #test for argument "linked = T"
  conf_95_fitted_linked = matrix(c(98.2755539, 0.0564122, 9.4921905, 46.6531090, 8.3263254,
                                   106.5003050, 0.6164591, 10.5921762, 102.1794778, 10.7608676), ncol = 2)
  o = confint(fit, types = 'fitted', linked = TRUE, level = 0.95)[pars_names,]
  relative.error = max(abs((o - conf_95_fitted_linked)/conf_95_fitted_linked))
  expect_equal(relative.error, 0, tolerance = 1e-4)
  
  ############################################################################################################
  #since this model has no parameter been extended, skip the test for 
  
})

