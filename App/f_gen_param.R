

###################### Parameter generation #######################################

# Load csv files with cost, utility, transition data and life tables to make available in environment
costs <- read.csv("Data/costs.csv", header = TRUE)

utilities <- read.csv("Data/utilities.csv", header = TRUE)

transitions <- read.csv("Data/transitions.csv", header = TRUE)
colnames(transitions) <- c("From", "mUMSARS 1 + 2", "3", "4", "5", "6 (Death)")
rownames(transitions) <- transitions$From
transitions$From <- NULL

# Load and process life table data
lifetable <- read.csv("Data/lifetable.csv", header = TRUE)

lifetable <- lifetable %>%
  rename_with(
    ~ str_replace(.x, "X80\\.\\.years\\.old", "X80.100.years.old"),
    starts_with("X80..years.old"))

lifetable <- lifetable %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "age_col",
    values_to = "hazard_rate"
  ) %>%
  mutate(
    # Remove leading "X" and trailing ".years.old"
    age_col = str_remove(age_col, "^X"),
    age_col = str_remove(age_col, "\\.years\\.old$")
  ) %>%
  separate(
    age_col,
    into = c("start_age", "end_age"),
    sep = "\\.",
    fill = "right",
    extra = "drop"
  ) %>%
  mutate(
    start_age = as.numeric(start_age),
    end_age   = as.numeric(end_age),
    # Treat missing end_age as start_age so there's at least one valid value
    end_age   = ifelse(is.na(end_age), start_age, end_age)
  ) %>%
  filter(is.finite(start_age), is.finite(end_age)) %>%
  rowwise() %>%
  mutate(age = list(seq(start_age, end_age))) %>%
  unnest(age)

lifetable <- lifetable[ , !(names(lifetable) %in% c("Code", "start_age", "end_age")) ]

# Convert death rates to monthly probabilities
lifetable <- lifetable %>%
  mutate(
    p_Dage = 1 - (1 - (hazard_rate / 1000))^(1/12)
  )

################## Define all model parameters ##############################
f_gen_param <- function() {
  # This function generates all the parameters needed for the model
  
# Define constant parameters
n_age_init <- 56 # age at baseline
n_age_end <- 100 # age at end of model
n_cycles <- (n_age_end - n_age_init) * 12 # monthly cycles
v_names_states <- c("Ci", "Md", "Vd", "Td", "D") # names of health states
n_states <- length(v_names_states) # number of health states
d_c <- 0.035 # discount rate for costs
d_e <- 0.035 # discount rate for effects (QALYs)
n_sim <- 1000 # number of simulations for PSA

# Define strategy names
strategies <- c("BSC", "Hypothetical Intervention & BSC")

# Generate all-cause mortality vector
v_p_Dage <- rep(lifetable$p_Dage[lifetable$age >= n_age_init & lifetable$age < n_age_end], each = 12)[1:n_cycles]

# Extract transition probabilities
p_CiCi <- transitions["mUMSARS 1 + 2", "mUMSARS 1 + 2"]
p_CiMd <- transitions["mUMSARS 1 + 2", "3"]
p_CiVd <- transitions["mUMSARS 1 + 2", "4"]
p_CiTd <- transitions["mUMSARS 1 + 2", "5"]
p_CiD <- transitions["mUMSARS 1 + 2", "6 (Death)"]

p_MdMd <- transitions["3", "3"]
p_MdVd <- transitions["3", "4"]
p_MdTd <- transitions["3", "5"]
p_MdD <- transitions["3", "6 (Death)"]

p_VdVd <- transitions["4", "4"]
p_VdTd <- transitions["4", "5"]
p_VdD <- transitions["4", "6 (Death)"]

p_TdTd <- transitions["5", "5"]
p_TdD <- transitions["5", "6 (Death)"]

p_DD <- transitions["6 (Death)", "6 (Death)"]

# Define treatment modifier
rr_trt_mod_p <- 0.5 # monthly relative risk for treatment on transitions

# Extract cost parameters
c_Ci <- round(costs$Total[costs$UMSARS == "UMSARS_IV 1 + 2"], 2)
c_Md <- round(costs$Total[costs$UMSARS == "UMSARS_IV 3"], 2)
c_Vd <- round(costs$Total[costs$UMSARS == "UMSARS_IV 4"], 2)
c_Td <- round(costs$Total[costs$UMSARS == "UMSARS_IV 5"], 2)
c_D <- 0
c_trt <- 1000 # monthly treatment cost

# Extract utility parameters
u_Ci <- utilities$Mean_Index_Score[utilities$UMSARS_IV == "1 + 2"]
u_Md <- utilities$Mean_Index_Score[utilities$UMSARS_IV == "3"]
u_Vd <- utilities$Mean_Index_Score[utilities$UMSARS_IV == "4"]
u_Td <- utilities$Mean_Index_Score[utilities$UMSARS_IV == "5"]
u_D <- 0

# Define disutilities
du_UD <- 0.0054 # disutility for urinary disorders
du_HD <- 0.067  # disutility for hypotension disorders
du_BD <- 0.194  # disutility for bowel disorders

# Define event probabilities
p_UD <- 0.511 # annual probability of urinary disorders
p_HD <- 0.567 # annual probability of hypotension disorders
p_BD <- 0.582 # annual probability of bowel disorders

# Define treatment modifier for events
rr_trt_mod_e <- 0.5 # monthly relative risk for treatment on events

# Create a list of parameters 
params <- list(
  strategies = strategies,
  n_age_init = n_age_init,
  n_age_end = n_age_end,
  n_cycles = n_cycles,
  v_names_states = v_names_states,
  n_states = n_states,
  d_c = d_c,
  d_e = d_e,
  n_sim = n_sim,
  v_p_Dage = v_p_Dage,
  p_CiCi = p_CiCi,
  p_CiMd = p_CiMd,
  p_CiVd = p_CiVd,
  p_CiTd = p_CiTd,
  p_CiD = p_CiD,
  p_MdMd = p_MdMd,
  p_MdVd = p_MdVd,
  p_MdTd = p_MdTd,
  p_MdD = p_MdD,
  p_VdVd = p_VdVd,
  p_VdTd = p_VdTd,
  p_VdD = p_VdD,
  p_TdTd = p_TdTd,
  p_TdD = p_TdD,
  p_DD = p_DD,
  rr_trt_mod_p = rr_trt_mod_p,
  c_Ci = c_Ci,
  c_Md = c_Md,
  c_Vd = c_Vd,
  c_Td = c_Td,
  c_D = c_D,
  c_trt = c_trt,
  u_Ci = u_Ci,
  u_Md = u_Md,
  u_Vd = u_Vd,
  u_Td = u_Td,
  u_D = u_D,
  du_UD = du_UD,
  du_HD = du_HD,
  du_BD = du_BD,
  p_UD = p_UD,
  p_HD = p_HD,
  p_BD = p_BD,
  rr_trt_mod_e = rr_trt_mod_e
)

return(params)

} # End of function

