
# We define a function for the Markov model to which we pass either deterministic or probabilistic parameters.

f_MM_MSA <- function(params,
                     is_deterministic = FALSE) {
 
   
  
  with(as.list(params), {
    
    
    # Comprehensive debugging to check if parameters are passed correctly to model code, which the code indicates.

     # if(is_deterministic || (!is_deterministic && runif(1) < 0.1)) {
     # 
     # 
     #   # Show debug info for ~10% of runs
     #   cat("\n===== PARAMETER VALUE VERIFICATION =====\n")
     #   cat("Running mode:", ifelse(is_deterministic, "DETERMINISTIC", "PSA"), "\n\n")
     # 
     #   cat("===== KEY CALCULATION PARAMETERS =====\n")
     #   # Basic model parameters
     #   cat("n_cycles: from params=", params$n_cycles, "| direct=", n_cycles, "| USED=", n_cycles, "\n")
     #   cat("d_c: from params=", params$d_c, "| direct=", d_c, "| USED=", d_c, "\n")
     #   cat("d_e: from params=", params$d_e, "| direct=", d_e, "| USED=", d_e, "\n")
     # 
     #   # Cost parameters
     #   cat("\n===== COST PARAMETERS =====\n")
     #   cat("c_Ci: from params=", params$c_Ci, "| local=", c_Ci, "\n")
     #   cat("c_Md: from params=", params$c_Md, "| local=", c_Md, "\n")
     #   cat("c_Vd: from params=", params$c_Vd, "| local=", c_Vd, "\n")
     #   cat("c_Td: from params=", params$c_Td, "| local=", c_Td, "\n")
     #   cat("c_trt: from params=", params$c_trt, "| local=", c_trt, "\n")
     # 
     #   # Utility parameters
     #   cat("\n===== UTILITY PARAMETERS =====\n")
     #   cat("u_Ci: from params=", params$u_Ci, "| local=", u_Ci, "(annual equiv)\n")
     #   cat("u_Md: from params=", params$u_Md, "| local=", u_Md, "(annual equiv)\n")
     #   cat("u_Vd: from params=", params$u_Vd, "| local=", u_Vd, "(annual equiv)\n")
     #   cat("u_Td: from params=", params$u_Td, "| local=", u_Td, "(annual equiv)\n")
     # 
     #   # Treatment modifiers
     #   cat("\n===== TREATMENT MODIFIERS =====\n")
     #   cat("rr_trt_mod_p: from params=", params$rr_trt_mod_p, "| local=", rr_trt_mod_p, "\n")
     #   cat("rr_trt_mod_e: from params=", params$rr_trt_mod_e, "| local=", rr_trt_mod_e, "\n")
     # 
     #   cat("\n=====================================\n")
     # }
    

    # Flag that decides if deterministic or probabilistic analysis should run through model code. 
    
    params$is_deterministic <- is_deterministic
    
##########################  Parameter calculations  ###############################

# Adjust sampled parameters to be monthly. Check that only that the annual parameters are adjusted here.  

    # Convert annual utilities to monthly utilities
    
    u_Ci <- u_Ci / 12 # monthly utility for Ci
    u_Md <- u_Md / 12 # monthly utility for Md
    u_Vd <- u_Vd / 12 # monthly utility for Vd
    u_Td <- u_Td / 12 # monthly utility for Td
  
    # Convert annual disutilities to monthly disutilities
    
    du_UD <- du_UD / 12 # monthly disutility for urinary disorders
    du_HD <- du_HD / 12 # monthly disutility for hypotension disorders
    du_BD <- du_BD / 12 # monthly disutility for bowel disorders
    
    # Convert annual event probabilities to monthly probabilities
    
    p_UD <- 1 - (1 - p_UD) ^ (1/12) # monthly probability of urinary disorders
    p_HD <- 1 - (1 - p_HD) ^ (1/12) # monthly probability of hypotension disorders
    p_BD <- 1 - (1 - p_BD) ^ (1/12) # monthly probability of bowel disorders
  
# # Debugging statements to check the conversion of parameters from annual to monthly are working, which it is.
# 
#     print(paste("Original u_Ci:", params$u_Ci, "| Monthly u_Ci:", u_Ci, "\n"))
#     print(paste("Original p_UD:", params$p_UD, "| Monthly p_UD:", p_UD, "\n"))
#     print(paste("Original du_UD:", params$du_UD, "| Monthly du_UD:", du_UD, "\n"))

    # Calculate disutilities with and without event modifier
    
    disutilities <- c(p_UD * du_UD + p_HD * du_HD + p_BD * du_BD) # disutilities which will be applied to each health state in BSC arm
    disutilities_trt <- c(p_UD * du_UD * rr_trt_mod_e + p_HD * du_HD * rr_trt_mod_e + p_BD * du_BD * rr_trt_mod_e) # disutilities which will be applied to each health state in hypothetical intervention & BSC arm
    
    # Discount weights but first making them monthly
    
    d_c_m <- (1 + d_c) ^ (1/12) - 1 # discount rate for costs
    d_e_m <- (1 + d_e) ^ (1/12) - 1 # discount rate for effects
    
    
    d_wc <- 1 / (1 + d_c_m) ^ (0:(n_cycles)) # discount weight for costs
    d_we <- 1 / (1 + d_e_m) ^ (0:(n_cycles)) # discount weight for effects
    

################## Generate matrices for Markov model #########################  

# Generate transition matrices for BSC and hypothetical intervention & BSC    
    
# Create empty three dimensional array
    
    a_P <- array(0, dim = c(n_states, n_states, n_cycles),
                 dimnames = list(v_names_states, v_names_states, 0:(n_cycles-1)))
    
# Fill transition matrix with transition probabilities from generated parameters. 
    
   
    for(i in 1:n_cycles) {
      
      ### From Ci
      a_P["Ci","Ci",i] <- (1 - v_p_Dage[i]) * (1 - (p_CiMd + p_CiVd + p_CiTd + p_CiD))
      a_P["Ci","Md",i] <- (1 - v_p_Dage[i]) * p_CiMd
      a_P["Ci","Vd",i] <- (1 - v_p_Dage[i]) * p_CiVd  
      a_P["Ci","Td",i] <- (1 - v_p_Dage[i]) * p_CiTd
      a_P["Ci","D",i] <- (1 - v_p_Dage[i]) * p_CiD + v_p_Dage[i]
      
      ### From Md
      a_P["Md","Md",i] <- (1 - v_p_Dage[i]) * (1 - (p_MdVd + p_MdTd + p_MdD))
      a_P["Md","Vd",i] <- (1 - v_p_Dage[i]) * p_MdVd
      a_P["Md","Td",i] <- (1 - v_p_Dage[i]) * p_MdTd
      a_P["Md","D",i] <- (1 - v_p_Dage[i]) * p_MdD + v_p_Dage[i]
      
      ### From Vd
      a_P["Vd","Vd",i] <- (1 - v_p_Dage[i]) * (1 - (p_VdTd + p_VdD))
      a_P["Vd","Td",i] <- (1 - v_p_Dage[i]) * p_VdTd
      a_P["Vd","D",i] <- (1 - v_p_Dage[i]) * p_VdD + v_p_Dage[i]
      
      ### From Td
      a_P["Td","Td",i] <- (1 - v_p_Dage[i]) * (1 - p_TdD)
      a_P["Td","D",i] <- (1 - v_p_Dage[i]) * p_TdD + v_p_Dage[i]
      
      ### From D
      a_P["D","D",i] <- 1
    }
    
# Check if the number of cycles are equal and adds to 528
    # rowSums(a_P) 
    
# Create transition matrix for treatment modifier
    
    
    a_P_trt <- array(0, dim = c(n_states, n_states, n_cycles),
                     dimnames = list(v_names_states, v_names_states, 0:(n_cycles-1)))
    
# Fill transition matrix with transition probabilities with treatment modifier
    
    for(i in 1:n_cycles) {
      
      ### From Ci
      a_P_trt["Ci","Ci",i] <- (1 - v_p_Dage[i]) * (1 - (p_CiMd + p_CiVd + p_CiTd + p_CiD) * rr_trt_mod_p)
      a_P_trt["Ci","Md",i] <- (1 - v_p_Dage[i]) * p_CiMd * rr_trt_mod_p
      a_P_trt["Ci","Vd",i] <- (1 - v_p_Dage[i]) * p_CiVd * rr_trt_mod_p
      a_P_trt["Ci","Td",i] <- (1 - v_p_Dage[i]) * p_CiTd * rr_trt_mod_p
      a_P_trt["Ci","D",i] <- (1 - v_p_Dage[i]) * p_CiD * rr_trt_mod_p + v_p_Dage[i]
      
      ### From Md
      a_P_trt["Md","Md",i] <- (1 - v_p_Dage[i]) * (1 - (p_MdVd + p_MdTd + p_MdD) * rr_trt_mod_p)
      a_P_trt["Md","Vd",i] <- (1 - v_p_Dage[i]) * p_MdVd * rr_trt_mod_p
      a_P_trt["Md","Td",i] <- (1 - v_p_Dage[i]) * p_MdTd * rr_trt_mod_p
      a_P_trt["Md","D",i] <- (1 - v_p_Dage[i]) * p_MdD * rr_trt_mod_p + v_p_Dage[i]
      
      ### From Vd
      a_P_trt["Vd","Vd",i] <- (1 - v_p_Dage[i]) * (1 - (p_VdTd + p_VdD) * rr_trt_mod_p)
      a_P_trt["Vd","Td",i] <- (1 - v_p_Dage[i]) * p_VdTd * rr_trt_mod_p
      a_P_trt["Vd","D",i] <- (1 - v_p_Dage[i]) * p_VdD * rr_trt_mod_p + v_p_Dage[i]
      
      ### From Td
      a_P_trt["Td","Td",i] <- (1 - v_p_Dage[i]) * (1 - p_TdD)
      a_P_trt["Td","D",i] <- (1 - v_p_Dage[i]) * p_TdD + v_p_Dage[i]
      
      ### From D
      a_P_trt["D","D",i] <- 1
      
    }
    
# Check if the number of cycles are equal and adds to 528
    # rowSums(a_P_trt) 
    
# Create Markov trace matrix for BSC 
    
    m_TR <- matrix(data = NA, 
                   nrow = (n_cycles + 1),
                   ncol = n_states,
                   dimnames = list(0:n_cycles, v_names_states))
    
# Initial distribution across health states. Note this will be the same for both strategies.
    
    m_TR[1,] <- c("Ci" = 1, "Md" = 0, "Vd" = 0, "Td" = 0, D = 0) # Everyone starts in Ci state
    
# Check if first row in Markov has the above initial distribution. 
    
    # head(m_TR)
    
# Assigning initial distribution. 
    
    m_TR_trt <- m_TR
    
######################## Run Markov traces ####################################

#  Run Markov traces for BSC
    
    for (t in 1:n_cycles) {
      m_TR[t + 1,] <- m_TR[t,] %*% a_P[,,t]
    }
    
# Run Markov model for Hypothetical Intervention & BSC
    
    for (t in 1:n_cycles) {
      m_TR_trt[t + 1,] <- m_TR_trt[t,] %*% a_P_trt[,,t]
    }

################ Internal model validation ##########################
    
    # Function to validate both traces and transition matrices
    check_model_integrity <- function(m_TR, m_TR_trt, a_P, a_P_trt, tolerance = 1e-6) { # Adding small tolerance to account for rounding errors
      # --- Check Markov traces ---
      bsc_trace_sums <- rowSums(m_TR)
      bsc_trace_valid <- all(abs(bsc_trace_sums - 1) < tolerance)
      
      trt_trace_sums <- rowSums(m_TR_trt)
      trt_trace_valid <- all(abs(trt_trace_sums - 1) < tolerance)
      
      # --- Check transition matrices (sample cycle 1) ---
      bsc_matrix_sums <- rowSums(a_P[,,1])
      bsc_matrix_valid <- all(abs(bsc_matrix_sums - 1) < tolerance)
      
      trt_matrix_sums <- rowSums(a_P_trt[,,1])
      trt_matrix_valid <- all(abs(trt_matrix_sums - 1) < tolerance)
      
      # Only print if there's an issue or in deterministic mode
      if (!bsc_trace_valid || !trt_trace_valid || !bsc_matrix_valid || !trt_matrix_valid || is_deterministic) {
        cat("\n=== Internal MODEL VALIDATION ===\n")
        
        # Trace validation results
        cat("TRACE VALIDATION:\n")
        cat("- BSC trace sums to 1:", ifelse(bsc_trace_valid, "Yes ✓", "No ✗"), "\n")
        cat("- TRT trace sums to 1:", ifelse(trt_trace_valid, "Yes ✓", "No ✗"), "\n")
        
        # Transition matrix validation results
        cat("\nTRANSITION MATRIX VALIDATION:\n")
        cat("- BSC matrix rows sum to 1:", ifelse(bsc_matrix_valid, "Yes ✓", "No ✗"), "\n")
        cat("- TRT matrix rows sum to 1:", ifelse(trt_matrix_valid, "Yes ✓", "No ✗"), "\n")
        
        # Report specific issues if any
        if (!bsc_trace_valid) {
          problem_cycles <- which(abs(bsc_trace_sums - 1) >= tolerance)
          cat("\nBSC trace issues in cycles:", head(problem_cycles), "...\n")
        }
        
        if (!trt_trace_valid) {
          problem_cycles <- which(abs(trt_trace_sums - 1) >= tolerance)
          cat("TRT trace issues in cycles:", head(problem_cycles), "...\n")
        }
        
        if (!bsc_matrix_valid) {
          problem_rows <- which(abs(bsc_matrix_sums - 1) >= tolerance)
          cat("BSC matrix issues in rows:", problem_rows, 
              "(sums:", round(bsc_matrix_sums[problem_rows], 6), ")\n")
        }
        
        if (!trt_matrix_valid) {
          problem_rows <- which(abs(trt_matrix_sums - 1) >= tolerance)
          cat("TRT matrix issues in rows:", problem_rows,
              "(sums:", round(trt_matrix_sums[problem_rows], 6), ")\n")
        }
        
        cat("======================\n")
      }
      
      # Return overall validation status
      return(bsc_trace_valid && trt_trace_valid && bsc_matrix_valid && trt_matrix_valid)
    }
    
    # Execute validation after creating matrices and running traces
    model_valid <- check_model_integrity(m_TR, m_TR_trt, a_P, a_P_trt)
  
####################### Data for plot of Markov trace for deterministic model for use in app ##############
    
    # Compare Markov traces for BSC and Hypothetical Intervention & BSC by converting matrices to data frames for ggplot
    
    df_TR <- as.data.frame(m_TR)
    df_TR$Cycle <- 1:nrow(df_TR)
    df_TR$Treatment <- "BSC"
    
    df_TR_trt <- as.data.frame(m_TR_trt)
    df_TR_trt$Cycle <- 1:nrow(df_TR_trt)
    df_TR_trt$Treatment <- "Hypothetical Intervention & BSC"
    
    # Combine the data frames
    df_traces <- rbind(df_TR, df_TR_trt)
    
    # Combine data frames for ggplot
    df_traces_melt <- melt(df_traces, id.vars = c("Cycle", "Treatment"), variable.name = "HealthState", value.name = "Proportion")
    
    # We limit to 350 circles for better illustration.
    df_traces_melt_lim <- subset(df_traces_melt, Cycle <= 350)
    
    # Data that show the cumulative state occupancy over time.
    
    df_traces_melt_lim$HealthState <- factor(df_traces_melt_lim$HealthState, 
                                             levels = c("D", "Td", "Vd", "Md", "Ci"))
    
    
  ################### Calculate costs and QALYs #################################
  
# Calculate costs and QALYs
  
# Summarise costs and effects by vectors
  
  v_c_BSC <- c("Ci" = c_Ci, "Md" = c_Md, "Vd" = c_Vd, "Td" = c_Td, "D" = 0) # Costs for BSC
  v_u_BSC <- c("Ci" = u_Ci - disutilities, "Md" = u_Md - disutilities, "Vd" = u_Vd - disutilities, "Td" = u_Td - disutilities, "D" = 0) # QALYs for BSC
  
  v_c_trt <- c("Ci" = c_Ci + c_trt, "Md" = c_Md + c_trt, "Vd" = c_Vd + c_trt, "Td" = c_Td, "D" = 0) # Costs for Hypothetical Drug & BSC
  v_u_trt <- c("Ci" = u_Ci - disutilities_trt, "Md" = u_Md - disutilities_trt, "Vd" = u_Vd - disutilities_trt, "Td" = u_Td - disutilities, "D" = 0) # QALYs for Hypothetical Drug & BSC
  
  # Multiply vectors by the markov trace (proportion in each state by time)
  
  v_E_BSC <- (m_TR %*% v_u_BSC) # QALYs for BSC
  v_C_BSC <- (m_TR %*% v_c_BSC) # Costs for BSC
  
  v_E_trt <- (m_TR_trt %*% v_u_trt) # QALYs for Hypothetical Drug & BSC
  v_C_trt <- (m_TR_trt %*% v_c_trt) # Costs for Hypothetical Drug & BSC
  
  # Apply discount rates
  
  v_E_BSC_dis <- sum(v_E_BSC * d_we) # Discounted QALYs for BSC
  v_C_BSC_dis <- sum(v_C_BSC * d_wc) # Discounted costs for BSC
  
  v_E_trt_dis <- sum(v_E_trt * d_we) # Discounted QALYs for Hypothetical Drug & BSC
  v_C_trt_dis <- sum(v_C_trt * d_wc) # Discounted costs for Hypothetical Drug & BSC
  
  
  ################ Results #####################################
  
  # Calculate incremental costs and QALYs (Incremental cost effectiveness ratio)
  
  results <- c(
    "Costs_BSC" = v_C_BSC_dis,
    "Costs_trt" = v_C_trt_dis,
    "QALYs_BSC" = v_E_BSC_dis,
    "QALYs_trt" = v_E_trt_dis,
    "ICER" = (v_C_trt_dis - v_C_BSC_dis) / (v_E_trt_dis - v_E_BSC_dis)
  )
  
########### Create deterministic cost breakdown for use in app ##########################
  
# Define full labels. Instead of abbreviations from v_states_names
  
  state_labels <- c("Completely or not completely independent", 
                    "More dependent", 
                    "Very dependent", 
                    "Totally dependent", 
                    "Death")
  
  # For BSC strategy: calculate discounted costs and QALYs per health state
  cost_breakdown_BSC <- sapply(v_names_states, function(s) {
    sum(m_TR[, s] * v_c_BSC[s] * d_wc)
  })
  qalys_breakdown_BSC <- sapply(v_names_states, function(s) {
    sum(m_TR[, s] * v_u_BSC[s] * d_we)
  })
  
  # For Hypothetical Intervention & BSC strategy: calculate discounted costs and QALYs per health state
  cost_breakdown_trt <- sapply(v_names_states, function(s) {
    sum(m_TR_trt[, s] * v_c_trt[s] * d_wc)
  })
  qalys_breakdown_trt <- sapply(v_names_states, function(s) {
    sum(m_TR_trt[, s] * v_u_trt[s] * d_we)
  })
  
  # Create a data frame with the breakdown using the custom labels and updated column names.
  df_breakdown <- data.frame(
    `Health states` = state_labels,
    `Costs for BSC` = round(cost_breakdown_BSC, 0),
    `QALYs for BSC` = round(qalys_breakdown_BSC, 2),
    `Costs for hypothetical intervention and BSC` = round(cost_breakdown_trt, 0),
    `QALYs for hypothetical intervention and BSC` = round(qalys_breakdown_trt, 2),
    check.names = FALSE
  )
  
  if(is_deterministic) {
    # For deterministic analysis, return both trace data and results
    return(list(
      deterministic = results,
      df_traces_melt_lim = df_traces_melt_lim,
      breakdown = df_breakdown
    ))
  } else {
    # For PSA, just return results
    return(results)
  } # End of return(results)
  
  }) # End of with(as.list(params)
  
} # End of f_MM_MSA function
  
