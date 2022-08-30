#' Title
#'
#' @param sim_name 
#' @param n.rand 
#' @param fit 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
sim_study = function(sim_name, n.rand = 1, fit = FALSE, seed = 810){
  sim_args = sim_args_generator(sim_name)
  sim_args$n.rand = n.rand
  set.seed(seed)
  if(n.rand != 1) message("simulation progress:")
  simulated_capt = do.call('sim.capt', sim_args)
  output = list()
  output$sim_capt = simulated_capt$capt
  output$sim_args = sim_args
  
  if(fit){
    message("modelling started.")
    #sort out arguments required for model fitting first (excluding capt, as this could vary depends on n.rand)
    fit_args = fit_args_generator_from_sim(sim_name, simulated_capt$args)
    
    #if we only run 1 simulation and fit model to it, we could store the output of fit
    #if we run more simulations and fit models to all of them, we could not store that many
    #output objects, so we only plot the result of linked coefficients in a histogram with the
    #"true" value marked as a red vertical line
    if(n.rand == 1){
      fit_args$capt = create.capt(simulated_capt$capt, fit_args$traps)
      sim_fit = do.call('fit_og', fit_args)
      
      output$sim_fit_args = fit_args
      output$sim_fit = sim_fit
    } else {
      #generate the data frame which contains default link function for each parameter
      dat_par = default_df_link()
      
      #get the linked values for all parameters
      true_values = param_transform(sim_args$param, dat_par)
      est_values = vector('list', n.rand)
      fit_args$tracing = FALSE
      
      for(i in 1:n.rand){
        fit_args$capt = create.capt(simulated_capt$capt[[i]], fit_args$traps)
        sim_fit = do.call('fit_og', fit_args)
        est_values[[i]] = get_coef(sim_fit)
        message(paste0("finished: ", i, "/", n.rand))
        #write.csv(sim_fit$args$par.extend$data$mask, paste0('df_m_fit', i, '.csv'), row.names = F)
      }
      
      
      #remove the capture history and output the rest of arguments for model fitting
      fit_args$capt = NULL
      output$sim_fit_args = fit_args
      #we only output the linked scale of estimated coefficients
      output$sim_fit_coef_link = est_values
      
      
      #plot the est_values

      for(i in names(true_values)){
        for(j in 1:length(true_values[[i]])){
          #extract estimations as a vector and the true value
          est = sapply(est_values, function(x) x[[i]][j])
          tru = true_values[[i]][j]
          
          hist(est, xlim = range(est, tru),
               main = paste0(i, '[', j, ']_link: Sim vs. True'),
               xlab = 'Value')
          
          abline(v = tru, col = 2)
          legend('topright', 'true value', col = 2, lty = 1)
          box()
        }
      }
      
    }
  }
  
  return(output)
  
}
