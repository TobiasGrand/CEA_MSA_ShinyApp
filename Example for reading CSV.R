# Code to read a CSV file (PSA_all_parameters) in R

df_all_params <- read.csv(
  file = "PSA_all_parameters_20250813.csv", # Modify the file path to match yours
  header = TRUE,
  sep = ",",
  dec = "."
)

head(df_all_params)  # Display the first few rows of the dataframe

