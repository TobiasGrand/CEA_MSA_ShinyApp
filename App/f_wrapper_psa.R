
################################# PSA wrapper ###################################

# Create a wrapper function for the PSA model to allow for easy input of parameters and output of results
f_wrapper_psa <- function(
    ############ Adjustable USER INPUTS ############
    n_sim = 1000,
    d_c = 0.035,
    d_e = 0.035,
    rr_trt_mod_p = 0.5,
    rr_trt_mod_e = 0.5,
    c_Ci = 987.75,
    c_Md = 2000.19,
    c_Vd = 3679.37,
    c_Td = 3679.37,
    c_trt = 1000,
    progress = NULL)
 {
  
  
  ###### Unadjustable inputs
  
  base_params <- f_gen_param()
  
  #################################################################################
  
  ########## Run PSA ###############
  
  # Create inputs for PSA
  df_psa <- f_gen_psa(n_sim = n_sim, c_trt = c_trt, c_Ci = c_Ci, c_Md = c_Md, c_Vd = c_Vd, c_Td = c_Td)
  
  # Create matrix to store results 
  m_out <- matrix(data = NaN, 
                  nrow = n_sim,
                  ncol = 5,
                  dimnames = list(1:n_sim, c("Cost_BSC", "Cost_trt", "QALY_BSC", "QALY_trt", "ICER")))
  
  update_freq <- max(1, floor(n_sim/20)) # Update progress every 5% of iterations
  
  # Run model for each simulation
  for (i in 1:n_sim) {
   
     # Update progress every update_freq iterations
    if (i %% update_freq == 0 || i == n_sim) {
      # Using withProgress in the calling function. See action button where withProgress is located.
      shiny::setProgress(value = i/n_sim, 
                         detail = paste0(i, " of ", n_sim, " simulations (", 
                                         round(100*i/n_sim), "%)"))
    }
    
    # Flag that this is not deterministic, but PSA analysis. 
    
    psa_params <- as.list(df_psa[i,])
    psa_params$is_deterministic <- FALSE
    
    # User-specified parameters first that are not in the PSA data frame
            
    psa_params$c_trt <- c_trt  # Explicitly pass the user-specified treatment cost
    psa_params$rr_trt_mod_p <- rr_trt_mod_p  # Treatment effect on progression
    psa_params$rr_trt_mod_e <- rr_trt_mod_e  # Treatment effect on events
    psa_params$d_c <- d_c  # Discount rate for costs
    psa_params$d_e <- d_e  # Discount rate for effects
    
    # Add all missing parameters from base_params (only once per iteration)
    for (param_name in names(base_params)) {
      if (is.null(psa_params[[param_name]])) {
        psa_params[[param_name]] <- base_params[[param_name]]
      }
    }
    
    m_out[i,] <- f_MM_MSA(psa_params,
                          is_deterministic = FALSE)  
  } # End loop
  
  # Convert matrix to data frame for plotting
  df_out <- as.data.frame(m_out)
  
  # Return PSA results and the PSA input data
  return(list(
    results = df_out,
    psa_data = df_psa
  ))
  
}  # End of f_wrapper_psa function

