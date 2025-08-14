################################# PSA wrapper ###################################

# Create a wrapper function for t
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
  
  # Create a list to store all parameters for later download and debugging. Important that is a dataframe. 
  all_parameters <- list()
  
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
      if (!param_name %in% names(psa_params)) {
        psa_params[[param_name]] <- base_params[[param_name]]
      }
    }
    
    ############# Internal validation of parameters ############
    
    # Verify all model parameters are correctly transferred to PSA. We use the first iteration for this. 
    if (i == 1) {  
      missing_params <- setdiff(names(base_params), names(psa_params))
      param_valid <- length(missing_params) == 0
      expected_count <- length(names(base_params)) + 1  # base params + is_deterministic flag (FALSE)
      
      cat("\n=== PARAMETER VALIDATION (PSA) ===\n")
      cat("- Parameters count:", length(names(psa_params)), 
          "(Expected: 44 + is_deterministic flag (FALSE))\n")
      
      if (!param_valid || length(names(psa_params)) != expected_count) {
        warning("Unexpected parameter count in PSA analysis")
        cat("- Parameter validation:", "FAILED ✗\n")
        
        if (!param_valid) {
          cat("\nMISSING PARAMETERS:\n")
          cat(paste(missing_params, collapse=", "), "\n")
          warning("PSA may produce unreliable results due to missing parameters")
        }
      } else {
        cat("- Parameter validation:", "PASSED ✓\n")
      }
      cat("======================\n\n")
    }
    ########################
    
    # Store all parameters for debugging and download
    all_parameters[[i]] <- psa_params
    
    
    m_out[i,] <- f_MM_MSA(psa_params,
                          is_deterministic = FALSE)  
  } # End loop
  
  # Convert to a data frame
  df_all_params <- do.call(rbind, lapply(1:n_sim, function(i) {
    # Convert list to data frame row with iteration number
    df <- as.data.frame(t(unlist(all_parameters[[i]])))
    df$iteration <- i
    return(df)
  }))
  
  
  # Convert matrix to data frame for plotting
  df_out <- as.data.frame(m_out)
  
  # Return PSA results and the PSA input data
  return(list(
    results = df_out,
    psa_data = df_psa,
    all_parameters = df_all_params
  ))
  
}  # End of f_wrapper_psa function