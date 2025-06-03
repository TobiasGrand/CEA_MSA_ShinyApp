# We shape distributions first from available information. 

# We approximate the distribution of the event probabilities with beta distribution from the publication: http://dx.doi.org/10.1016/S1474-4422(12)70327-7. 
# For example, alpha = events + 1 and beta = no events + 1, which is a simple bayesian approach that treats the distribution as a binominal prior. 
# Also, referred to as the pseudo count rule or rule of succession.

alpha_p_UD <- 72 + 1 # events per study population
beta_p_UD <- 69 + 1 # no events per study population

alpha_p_HD <- 80 + 1 # events per study population
beta_p_HD <- 61 + 1 # no events per study population

alpha_p_BD <- 82 + 1 # events per study population
beta_p_BD <- 59 + 1 # no events per study population

# Example of how to derive alpha and beta values for du_UD from publication: https://doi.org/10.1177/0272989X11401031. 

mean_du_UD <- 0.0054 # From reference paper ": Standard error = 0.0072, Mean = -0.0054.
se_du_UD <- 0.0072

alpha_du_UD <- ((1 - mean_du_UD) / se_du_UD^2 - 1 / mean_du_UD) * mean_du_UD^2
beta_du_UD <- alpha_du_UD * (1 / mean_du_UD - 1)

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
    
    p_CiMd = rbeta(n_sim, shape1 = 0.0515, shape2 = 0.5458), # Alpha and beta values calculated from multi state modelling for those transitions allowed in Q matrix. 
   
    # Several transitions are not allowed in Q intensity matrix during multistate modelling, which is why they do not feature here. 
   
    p_MdVd = rbeta(n_sim, shape1 = 0.1114, shape2 = 1.045),

   #p_MdD was kept constant because mean is close to zero and alpha zero. 
   
    p_VdTd = rbeta(n_sim, shape1 = 0.1102, shape2 = 3.5233),
    p_VdD = rbeta(n_sim, shape1 = 0.07, shape2 = 3.4874),
   
    p_TdD = rbeta(n_sim, shape1 = 0.0343, shape2 = 1.059),
    
    # Event probabilities

    p_UD =  rbeta(n_sim, shape1 = alpha_p_UD, shape2 = beta_p_UD), 
   
    p_HD =  rbeta(n_sim, shape1 = alpha_p_HD, shape2 = beta_p_HD),
   
    p_BD =  rbeta(n_sim, shape1 = alpha_p_BD, shape2 = beta_p_BD),
    
    # Costs
    ## We use gamma distribution, which is routinely used for cost data to reflect the right skewed nature of costs. We use shape parameter 4 and scale parameter mean/4 to maintain the mean while reducing extreme values.
    
    
    c_Ci = rgamma(n_sim, shape = 4, scale = c_Ci/4), # No data to inform distributions. 
    
    c_Md = rgamma(n_sim, shape = 4, scale = c_Md/4),
    
    c_Vd = rgamma(n_sim, shape = 4, scale = c_Vd/4),
    
    c_Td = rgamma(n_sim, shape = 4, scale = c_Td/4),
    
    
    # Utilities
      ## For sanity check: Run in console and check with min() if they are within plausible range of EQ-5D-3L (-0.594) in UK. 
    
    u_Ci = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "1 + 2"], shape2 = utilities$beta[utilities$UMSARS_IV == "1 + 2"] ), # Alpha and beta derived from patient level data.
    u_Md = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "3"], shape2 = utilities$beta[utilities$UMSARS_IV == "3"] ),   
    u_Vd = rbeta(n = n_sim, shape1 = utilities$alpha[utilities$UMSARS_IV == "4"], shape2 = utilities$beta[utilities$UMSARS_IV == "4"] ),
    
    # Disutilities
    
    du_UD = rbeta(n_sim, shape1 = alpha_du_UD, shape2 = beta_du_UD) 
    
   # du_HD was derived. No data to inform distributions.
   # du_BD was derived. No data to inform distributions.
  )
  
  return(df_psa)
}



# # Code for manual checking of visuals included in suppl. materials.
# 
# visual_psa <- f_gen_psa(
#   n_sim = 1000,  # or whatever number you prefer
#   d_c = 0.035,
#   d_e = 0.035,
#   rr_trt_mod_p = 0.5,
#   rr_trt_mod_e = 0.5,
#   c_Ci = 987.75,
#   c_Md = 2000.19,
#   c_Vd = 3679.37,
#   c_Td = 3679.37,
#   c_trt = 1000
# )
# 
# ## For transition probabilities
# 
# # p_CiMd
# plot(density(visual_psa$p_CiMd), main = "p_CiMd")
# visual_psa %>%
#   summarise(mean = mean(p_CiMd),
#             sd = sd(p_CiMd),
#             min = min(p_CiMd),
#             max = max(p_CiMd))
# p_CiMd <- mean(visual_psa$p_CiMd)
# print(p_CiMd)
# 
# # p_MdVd
# plot(density(visual_psa$p_MdVd), main = "p_MdVd")
# visual_psa %>%
#   summarise(mean = mean(p_MdVd),
#             sd = sd(p_MdVd),
#             min = min(p_MdVd),
#             max = max(p_MdVd))
# p_MdVd <- mean(visual_psa$p_MdVd)
# print(p_MdVd)
# 
# # p_VdTd
# plot(density(visual_psa$p_VdTd), main = "p_VdTd")
# visual_psa %>%
#   summarise(mean = mean(p_VdTd),
#             sd = sd(p_VdTd),
#             min = min(p_VdTd),
#             max = max(p_VdTd))
# p_VdTd <- mean(visual_psa$p_VdTd)
# print(p_VdTd)
# 
# # p_VdD
# plot(density(visual_psa$p_VdD), main = "p_VdD")
# visual_psa %>%
#   summarise(mean = mean(p_VdD),
#             sd = sd(p_VdD),
#             min = min(p_VdD),
#             max = max(p_VdD))
# p_VdD <- mean(visual_psa$p_VdD)
# print(p_VdD)
# 
# # p_TdD
# plot(density(visual_psa$p_TdD), main = "p_TdD")
# visual_psa %>%
#   summarise(mean = mean(p_TdD),
#             sd = sd(p_TdD),
#             min = min(p_TdD),
#             max = max(p_TdD))
# p_TdD <- mean(visual_psa$p_TdD)
# print(p_TdD)
# 
# ## For event probabilities
# 
# # p_UD
# plot(density(visual_psa$p_UD), main = "p_UD")
# visual_psa %>%
#   summarise(mean = mean(p_UD),
#             sd = sd(p_UD),
#             min = min(p_UD),
#             max = max(p_UD))
# p_UD <- mean(visual_psa$p_UD)
# print(p_UD)
# 
# # p_HD
# plot(density(visual_psa$p_HD), main = "p_HD")
# visual_psa %>%
#   summarise(mean = mean(p_HD),
#             sd = sd(p_HD),
#             min = min(p_HD),
#             max = max(p_HD))
# p_HD <- mean(visual_psa$p_HD)
# print(p_HD)
# 
# # p_BD
# plot(density(visual_psa$p_BD), main = "p_BD")
# visual_psa %>%
#   summarise(mean = mean(p_BD),
#             sd = sd(p_BD),
#             min = min(p_BD),
#             max = max(p_BD))
# p_BD <- mean(visual_psa$p_BD)
# print(p_BD)
# 
# ## For costs
# 
# # c_Ci
# plot(density(visual_psa$c_Ci), main = "c_Ci")
# visual_psa %>%
#   summarise(mean = mean(c_Ci),
#             sd = sd(c_Ci),
#             min = min(c_Ci),
#             max = max(c_Ci))
# c_Ci <- mean(visual_psa$c_Ci)
# print(c_Ci)
# 
# # c_Md
# plot(density(visual_psa$c_Md), main = "c_Md")
# visual_psa %>%
#   summarise(mean = mean(c_Md),
#             sd = sd(c_Md),
#             min = min(c_Md),
#             max = max(c_Md))
# c_Md <- mean(visual_psa$c_Md)
# print(c_Md)
# 
# # c_Vd
# plot(density(visual_psa$c_Vd), main = "c_Vd")
# visual_psa %>%
#   summarise(mean = mean(c_Vd),
#             sd = sd(c_Vd),
#             min = min(c_Vd),
#             max = max(c_Vd))
# c_Vd <- mean(visual_psa$c_Vd)
# print(c_Vd)
# 
# # c_Td
# plot(density(visual_psa$c_Td), main = "c_Td")
# visual_psa %>%
#   summarise(mean = mean(c_Td),
#             sd = sd(c_Td),
#             min = min(c_Td),
#             max = max(c_Td))
# c_Td <- mean(visual_psa$c_Td)
# print(c_Td)
# 
# ## For utilities
# 
# # u_Ci
# plot(density(visual_psa$u_Ci), main = "u_Ci")
# visual_psa %>%
#   summarise(mean = mean(u_Ci),
#             sd = sd(u_Ci),
#             min = min(u_Ci),
#             max = max(u_Ci))
# u_Ci <- mean(visual_psa$u_Ci)
# print(u_Ci)
# 
# # u_Md
# plot(density(visual_psa$u_Md), main = "u_Md")
# visual_psa %>%
#   summarise(mean = mean(u_Md),
#             sd = sd(u_Md),
#             min = min(u_Md),
#             max = max(u_Md))
# u_Md <- mean(visual_psa$u_Md)
# print(u_Md)
# 
# # u_Vd
# plot(density(visual_psa$u_Vd), main = "u_Vd")
# visual_psa %>%
#   summarise(mean = mean(u_Vd),
#             sd = sd(u_Vd),
#             min = min(u_Vd),
#             max = max(u_Vd))
# u_Vd <- mean(visual_psa$u_Vd)
# print(u_Vd)
# 
# ## For disutilities
# 
# # du_UD
# plot(density(visual_psa$du_UD), main = "du_UD")
# visual_psa %>%
#   summarise(mean = mean(du_UD),
#             sd = sd(du_UD),
#             min = min(du_UD),
#             max = max(du_UD))
# du_UD <- mean(visual_psa$du_UD)
# print(du_UD)
# 
# 


