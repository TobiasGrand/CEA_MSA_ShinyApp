# We define a function for parameters that need variation in Shiny App

f_gen_psa <- function(n_sim = 1000,
                      d_c = 0.035,
                      d_e = 0.035,
                      rr_trt_mod_p = 0.5,
                      rr_trt_mod_e = 0.5,
                      c_Ci = 987.75,
                      c_Md = 2000.19,
                      c_Vd = 3679.37,
                      c_Td = 3679.37,
                      c_trt = 1000){
  
  # We define a function for parameters and their distributions for inclusion in Probabilistic Sensitivity Analysis (PSA) 
  
  set.seed(11)
  
  
  df_psa <- data.frame(
    # Transitions probabilities
    
    ## Estimated distributions based on a sample size of 141 apart from cost data.
    
    p_CiMd = rbeta(n_sim, shape1 = 12.1, shape2 = 128.9), # alpha = n x mean and beta = n - alpha
    
    # Several transitions are not allowed in Q intensity matrix during multistate modelling, which is why they do not feature here. 
    
    p_MdVd = rbeta(n_sim, shape1 = 14.1, shape2 = 126.9),
    
    #p_MdD was kept constant because mean is close to zero and alpha zero. 
    
    p_VdTd = rbeta(n_sim, shape1 = 4.2, shape2 = 136.8),
    p_VdD = rbeta(n_sim, shape1 = 2.8, shape2 = 138.2),
    
    p_TdD = rbeta(n_sim, shape1 = 4.4, shape2 = 136.4),
    
    # Event probabilities
    
    p_UD =  rbeta(n_sim, shape1 = 72.1, shape2 = 68.9), 
    
    p_HD =  rbeta(n_sim, shape1 = 79.9, shape2 = 61.1),
    
    p_BD =  rbeta(n_sim, shape1 = 82.1, shape2 = 58.9),
    
    # Costs
    ## We use gamma distribution, which is routinely used for cost data to reflect the right skewed nature of costs.
    
    
    c_Ci = rgamma(n_sim, shape = c_Ci/50, scale = 50), # No data to inform distributions. 
    
    c_Md = rgamma(n_sim, shape = c_Md/50, scale = 50),
    
    c_Vd = rgamma(n_sim, shape = c_Vd/50, scale = 50),
    
    c_Td = rgamma(n_sim, shape = c_Td/50, scale = 50),
    
    
    # Utilities
    ## For sanity check: Run in console and check with min() if they are within plausible range of EQ-5D-3L (-0.594) in UK. 
    
    u_Ci = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "1 + 2"], shape2 = utilities$beta[utilities$UMSARS_IV == "1 + 2"] ), # Alpha and beta derived from patient level data.
    u_Md = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "3"], shape2 = utilities$beta[utilities$UMSARS_IV == "3"] ),   
    u_Vd = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "4"], shape2 = utilities$beta[utilities$UMSARS_IV == "4"] ),
    
    # Disutilities
    
    du_UD = rbeta(n_sim, shape1 = 0.7, shape2 = 140.3) 
    
    # du_HD was derived. No data to inform distributions.
    # du_BD was derived. No data to inform distributions.
  )
  
  return(df_psa)
}



