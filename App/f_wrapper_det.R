#################### Deterministic wrapper ####################################

# Create wrapper for deterministic model without user specification of parameters for a fixed base case
f_wrapper_det <- function() {
  
# Get parameters which were defined earlier
  
  params <- f_gen_param()
  
  # Flag that this is a deterministic analysis
  
  params$is_deterministic <- TRUE
  
  # Run through model as deterministic
  results <- f_MM_MSA(
    params = params,  
    is_deterministic = TRUE)
  
  # Format results for naming in app
  if (is.list(results)) {
    results_det <- results$deterministic
    df_traces_melt_lim <- results$df_traces_melt_lim
    breakdown <- results$breakdown
  } else {
    results_det <- results
    df_traces_melt_lim <- NULL
    breakdown <- NULL
  }
  
  # Return results
  return(list(
    df_traces_melt_lim = df_traces_melt_lim,
    results_det = results_det,
    breakdown = breakdown
  ))
}
