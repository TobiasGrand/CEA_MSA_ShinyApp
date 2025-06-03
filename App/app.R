        # Shiny app for a Markov Model in Multiple System Atrophy 
                
                      # Coded by Tobias Grand

####################### Load packages #########################################

# For externals - remember to change working directory to your local drive if running code. It is disabled now to allow the upload to shinyapp.io.

# setwd("C:/Users/togn/OneDrive - H. Lundbeck A S/University of Sheffield/Part III Cost-effectiveness modelling MSA/Markov Model MSA/Markov Model MSA Clean (wrappers)/App") # Set working directory

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(reshape2)
library(shiny)
library(bslib)
library(shinythemes)

###################### Source R scripts #######################################

# Then load the model functions
source("f_gen_param.R")
source("f_gen_psa.R")  
source("f_MM_MSA.R")
source("f_wrapper_det.R")
source("f_wrapper_psa.R")


########################### Load the data ##################################

deterministic_results <- f_wrapper_det()  # Generates deterministic results (the wrapper calls f_gen_param() internally)

########################## Define UI ######################################

ui <- navbarPage( # Opens the page navbar 
  # title 
  title = "Markov Model in Multiple System Atrophy",
  theme = bs_theme(bootswatch = "lumen"), # Theme of the app
  
  sidebar = bslib::sidebar( # Opens the sidebar
    width = 450, #adjusting width of sidebar, 
    
    # Heading run model 
    
    tags$h6(tags$strong("Run Model")),
    
    # Run model 
    actionButton(inputId = "run_model", label = "Run Model"),
    
    # space between inputs and action button
    br(),
    
    # Insert heading
    tags$h6(tags$strong("Model Inputs")),
    
    # Define number of simulations 
    sliderInput(inputId = "SI_n_sim", label = "Number simulations", min = 100, max = 10000, value = 1000),
    sliderInput(inputId = "SI_d_c", label = "Discount rate for costs", min = 0, max = 1, value = 0.035, step = 0.001),
    sliderInput(inputId = "SI_d_e", label = "Discount rate for effects", min = 0, max = 1, value = 0.035, step = 0.001),
    sliderInput(inputId = "SI_rr_trt_mod_p", label = "Slowing of disease progression: treatment modifier applied to transition probabilities", min = 0.01, max = 1, value = 0.5),
    sliderInput(inputId = "SI_rr_trt_mod_e", label = "Reduce clinical events: treatment modifier applied to event probabilities", min = 0.01, max = 1, value = 0.5),
    
    # cost heading
    tags$h6(tags$em("Monthly costs")), 
    
    sliderInput(inputId = "SI_c_Ci", label = "Cost of Ci (£)", min = 0, max = 100000, value = 987.75, step = 10),
    sliderInput(inputId = "SI_c_Md", label = "Cost of Md (£)", min = 0, max = 100000, value = 2000.19, step = 10),
    sliderInput(inputId = "SI_c_Vd", label = "Cost of Vd (£)", min = 0, max = 100000, value = 3679.37, step = 10),
    sliderInput(inputId = "SI_c_Td", label = "Cost of Td (£)", min = 0, max = 100000, value = 3679.37, step = 10),
    sliderInput(inputId = "SI_c_trt", label = "Cost of hypothetical intervention (£)", min = 0, max = 10000, value = 1000, step = 10)
 
     ), # Closes the sidebar
  
  nav_spacer(), # Adds space between the sidebar and the main panel
  
  
    # Tab 1: Model Structure
    nav_panel("Model Structure", tags$img(src = "ModelStructure.JPG", width = 1600, height = 1100),
              br(),
              tags$p("The image shows a Markov model structure for multiple system atrophy (MSA).",
              "The model has five health states: Ci (completely or not completely independent) Md (more dependent)",
              " Vd (very dependent) , Td (totally dependent), and D(death).",
              "Arrows indicate possible transitions between health states.")
              
              ), # Closing nav_panel
    
    # Tab 2: Deterministic analysis
    nav_panel ("Deterministic Analysis",
               tags$h4("Fixed results"),
               tags$p("Results for base-case analysis reported in publication"),
               uiOutput("SO_icer_Table_des"),
               br(), 
               tags$p("Proportion in each state over model cycles"),
               plotOutput("SO_markov_trace"),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               tags$p("Breakdown of costs and QALYs by health state"),
               uiOutput("SO_breakdown_Table")
               
               ), # Closing nav_panel
    
    # Tab 3: Sensitivity analysis
    nav_panel ("Sensitivity Analysis",
               tags$h4("Changeable results"),
               tags$p("The cost-effectiveness results for PSA simulations"), 
               uiOutput("SO_icer_Table"), 
               br(),
               tags$p("Cost-effectiveness plane"),
               plotOutput("SO_CE_plane"),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               br(),
               tags$p("Comparison of sampled (PSA) versus deterministic Values"),
               uiOutput("SO_PSA_table"),
               tags$p("The PSA data is available for download as csv file"),
               downloadButton("DownloadPSA", "Download PSA data")
               
               
               ), # Closing nav_panel
  
  # Tab 4: Modelling analysis plan
  nav_panel("Modelling Analysis Plan",
            tags$p("Here is the modelling analysis plan"),
            tags$iframe(style = "height:2000px; width:100%; scrolling = yes", 
                        src = "MAP_clean.pdf")
            
            ), # Closing nav_panel
  
  # Tab intro
  nav_panel("Introduction for Non-Technical Users",
            tags$h2("Welcome! This is an introduction for non-technical users."),
            br(),
            br(),
            tags$h4("Model Structure"),
            tags$p("We have created a model that reflects functional progression of Multiple System Atrophy (MSA). People in the model start as being completely or not completely independent, but over time their disease progress and they move to the health state more dependent. However, they can also die, both from natural causes and disease. The arrows indicate where people can move. For example, in any state people can die, but also remain in their state. The model change every month, meaning that there is a chance of moving or remaining at each monthly interval."),
            tags$h4("Deterministic analysis"),
            tags$p("The results in this section are fixed and will NOT change if model settings are changed. It is because this analysis reflects the analysis from our paper. The way the results are calculated is by adding together the quality of life and costs from each disease state over monthly periods. For example, in more severe disease states patients will use more healthcare resources and have worse quality of life. 
A treatment that reduces progression of disease will allow patients to stay in disease states with higher quality of life and lower costs for longer periods. However, a treatment does also have a cost, which means that the price of treatment is also accounted for in the model. The results table shows the costs of treatment and best supportive care (best available treatment before the introduction of a hypothetical intervention). Moreover, the quality of life is also captured, as a metric called QALY, which is a year in perfect health. It is now possible to calculate costs-effectiveness. Generally, a new treatment is more expensive, but it provides more quality of life to the patient. For example, by reducing progression of disease. 
If the incremental cost-effectiveness ratio (ICER) is between or below 20,000 to 30,000 pounds sterling it is regarded as cost-effective in the United Kingdom. However, in some scenarios, for example when costly time in hospital is avoided because of a new treatment, it turns out that the new interventions is both cheaper and provides more health. We usually then say that the new treatment is dominating existing treatments.
The illustration below the table with results shows how patients move and occupy the disease states over time. For example, if a treatment reduce progression to death, one will see that the coloured areas e.g., in the earlier health states will be larger, because they are occupied for longer before reaching death. 
"),
            tags$h4("Sensitivity Analysis"),
            tags$p("The sensitivity analysis tab shows results that can be changed by specifying parameters on the left-hand side. It is complicated to fully understand how all parameters play together without carefully examining the code or modelling analysis plan. Remember to press the button “Run Model” when something has been changed in the model. We suggest that non-health economist mostly vary costs or risk modifiers.
For example, increase or decrease costs associated with disease states or treatment. There is not a lot of data available for MSA, and data on costs of Vd and Td were aggregated. So, one could decide to increase cost of Td or the costs of the intervention, and then press “Run Model”. For slowing of disease progression e.g., by 0.5 means a factor of 0.5 is multiplied with the probability of moving, which means that the chance of progressing to a more severe disease state is halfed. Thus, if you set both reduction scrollers to 1, which means equal rate of transition and probability of having an event, the QALYs should be equal. 
We allow all parameters to vary across the chosen number of simulations, which is why results deviate from the deterministic results. We can see the different ICERs on the cost-effectiveness plane. Larger scattering of the dots (ICERS) illustrates larger uncertainty or variation in the analysis. To conclude, it may be better just to look at the table, because even for trained health economists, it may be difficult to comprehend without sufficient insights into the modelling code. 
"),
            tags$h4("Modelling Analysis Plan"),
            tags$p("We have provided a modelling analysis plan for those who want deeper insights into the parameters of the model. It is mostly, for people who are deeply interested in understanding how the modelling components play together and why certain variables have their assigned value."),
            tags$h4("Other resources"),
            tags$p("We have provided a link to our paper that discusses issues, challenges and opportunities for economic evaluations of orphan drugs in rare diseases. It is a good introduction to the topic of health econonomics for rare diseases"),
            
  ), # Closing nav_panel

  nav_menu( "Other Resources",
    nav_item(tags$a(href = "https://link.springer.com/article/10.1007/s40273-024-01370-2", "📄 Grand et al 2024: Issues, Challenges and Opportunities for Economic Evaluations of Orphan Drugs in Rare Diseases: An Umbrella Review"))
  ) # Closing nav_item
  
  
) # Close the navbarPage



########################## Define Server ##################################

server <- function(input, output, session) {
  
  
  ################################################################ Deterministic ###################################################################

  
  # Retrieve cost-effectiveness table (deterministic)
  output$SO_icer_Table_des <- renderUI({
   
    # Get results
    results_format <- deterministic_results$results_det
    
    # Format results to match required structure
    formatted_results <- data.frame(
      Strategy = c("Hypothetical Intervention & BSC", "BSC"),
      Costs = c(results_format["Costs_trt"], results_format["Costs_BSC"]),
      QALYs = c(results_format["QALYs_trt"], results_format["QALYs_BSC"]),
      ICER = c(NA, results_format["ICER"])
    )
    
    rownames(formatted_results) <- NULL
    
    # Round the results
    formatted_results$Costs <- round(formatted_results$Costs, 0)
    formatted_results$QALYs <- round(formatted_results$QALYs, 2)
    formatted_results$ICER <- round(formatted_results$ICER, 0)
    
    # Create the table
    shiny::HTML(
      kable(formatted_results, 
            caption = paste("ICER(£ per QALY) = ", format(round(formatted_results$ICER[2], 2), big.mark=",")), 
            align = "c", 
            format = "html") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = T)
    )
  })
  
  # Retrieve Markov trace 
  output$SO_markov_trace <- renderPlot({

    # Get the Markov trace data
    df_traces_melt_lim <- deterministic_results$df_traces_melt_lim
    
    # Create the plot
    ggplot(df_traces_melt_lim, aes(x = Cycle, y = Proportion, fill = HealthState)) +
      geom_area(color = "white", linewidth = 0.2) +
      facet_wrap(~ Treatment) +
      labs(title = "",
           x = "Cycle number",
           y = "Proportion of Patients in Health State") +
      scale_x_continuous(limits = c(0, 350)) +
      scale_fill_manual(
        breaks = c("Ci", "Md", "Vd", "Td", "D"),
        labels = c("Completely or not completely independent", 
                   "More dependent", 
                   "Very dependent", 
                   "Totally dependent", 
                   "Death"),
        values = c("Ci" = "#74c476",  # green
                   "Md" = "#a1d99b",  # lighter green
                   "Vd" = "#fd8d3c",  # orange
                   "Td" = "#e34a33",  # red
                   "D"  = "#969696")  # grey
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.position = "bottom")
  }, height = 800) # Adjusts height and closes the render plot function
  
  # Retrieve costs and QALYs by health state
  
  output$SO_breakdown_Table <- renderUI ({
    
    df_breakdown <- deterministic_results$breakdown
    
    rownames(df_breakdown) <- NULL
    
    ## Display the breakdown table using kable
    shiny::HTML(
      kable(df_breakdown, caption = "*Costs and QALYs are discounted. *Values are rounded", align = "c", format = "html") %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = T)
    )
    
  }) # Closes the render UI function
  
  ############################################################### PSA #############################################################################
  
  #Action button to run model
 observeEvent(input$run_model, ignoreNULL = F, {
   
   withProgress(message = "Running PSA simulations", value = 0, {
  
   
   # Run model function with Shiny Inputs
   
   results <- f_wrapper_psa(n_sim = input$SI_n_sim,
                            d_c = input$SI_d_c,
                            d_e = input$SI_d_e,
                            rr_trt_mod_p = input$SI_rr_trt_mod_p,
                            rr_trt_mod_e = input$SI_rr_trt_mod_e,
                            c_Ci = input$SI_c_Ci,
                            c_Md = input$SI_c_Md,
                            c_Vd = input$SI_c_Vd,
                            c_Td = input$SI_c_Td,
                            c_trt = input$SI_c_trt)
 
   #Store results
   df_model_res <- results$results
   df_model_psa <- results$psa_data
   
   # Download as csv for checking results and double coding in Excel. PSA results are also available from App/Data folder. 
   
   output$DownloadPSA <- downloadHandler(
     filename = function() {
       paste("PSA_results_", format(Sys.time(), "%Y%m%d"), ".csv", sep="")
     },
     content = function(file) {
       write.csv(df_model_psa, file, row.names = FALSE)
     }
   )
   
   ################################# CEA TABLE ###############################
   
   # Create cost-effectiveness table (mean of simulations)
   
   output$SO_icer_Table <- renderUI({
     
     # Create a table with the cost-effectiveness results
     df_res_table <- data.frame( 
       Strategy = c("Hypothetical Intervention & BSC", "BSC"),
       
       Costs = c(mean(df_model_res$Cost_trt), mean(df_model_res$Cost_BSC)),
       
       QALYs = c(mean(df_model_res$QALY_trt), mean(df_model_res$QALY_BSC)),
       
       
       ICER = c(NA, (mean(df_model_res$Cost_trt) - mean(df_model_res$Cost_BSC)) / (mean(df_model_res$QALY_trt) - mean(df_model_res$QALY_BSC)))
       
     ) # Closes dataframe for results table 
     
     # round dataframe results table 
     
     df_res_table$QALYs <- round(df_res_table$QALYs, 2)
     df_res_table$Costs <- round(df_res_table$Costs, 0)
     df_res_table$ICER <- round(df_res_table$ICER, 0)
     
     
     # Create a kable table with the results
     shiny::HTML( 
       kable(df_res_table, 
             caption = paste("ICER(£ per QALY) = ", format(round(df_res_table$ICER[2], 2), big.mark=",")), 
             align = "c", 
             format = "html") %>%
         kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                       full_width = T, 
                       position = "center"))
     
   }) # Closes the render table function

    ############################################################# Cost-effectiveness plane ########################################################
     
   # Create cost effectiveness plane 
   
   output$SO_CE_plane <- renderPlot ({
   
      # incremental cost and QALYs
     
     df_model_res$inc_C <- ((df_model_res$Cost_trt) - (df_model_res$Cost_BSC))
     
     df_model_res$inc_E <- ((df_model_res$QALY_trt) - (df_model_res$QALY_BSC))
     
     # Calculate symmetric limits around zero
     max_inc_E <- max(abs(df_model_res$inc_E))
     max_inc_C <- max(abs(df_model_res$inc_C))
     
     ggplot(df_model_res, aes(x = inc_E, y = inc_C)) +
       geom_point() +
       geom_hline(yintercept = 0, linetype = "dashed") +
       geom_vline(xintercept = 0, linetype = "dashed") +
       geom_abline(intercept = 0, slope = 30000, linetype = "dotted", color = "red") +
       annotate("text", x = -max_inc_E * 1, y = max_inc_C * -0.5, label = "£ 30,000 WTP threshold", color = "red", hjust = 0, vjust = 1, size = 5) +
       labs(
         title = "Cost-effectiveness plane",
         x = "Incremental QALYs",
         y = "Incremental Costs"
       ) +
       scale_x_continuous(
         labels = function(x) format(round(x, 2), scientific = FALSE),
         breaks = scales::pretty_breaks(n = 5)
       ) +
       scale_y_continuous(
         labels = function(x) paste0("£ ", format(round(x, 0), big.mark = ",", scientific = FALSE)),
         limits = c(-max_inc_C, max_inc_C),
         breaks = scales::pretty_breaks(n = 10)
       ) +
       theme_minimal(base_size = 12) +
       theme(
         legend.title = element_blank(),
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         axis.title = element_text(size = 12),
         axis.text = element_text(size = 14),
         legend.text = element_text(size = 14),
         legend.position = "bottom"
       )

   }, height = 800) # Adjust height and closes the render plot function
   
########################## PSA Table #########################################
   
   # Create table for PSA results
   
   output$SO_PSA_table <- renderUI({
     
     
     # Tabulate results to compare psa vs. determinisic values
     
     ## Summaries from df_psa
     vars <- c("p_CiMd", "p_MdVd", "p_VdTd", "p_VdD", "p_TdD",
               "p_UD", "p_HD", "p_BD",
               "c_Ci", "c_Md", "c_Vd", "c_Td",
               "u_Ci", "u_Md", "u_Vd",
               "du_UD")
     
     df_psa_summary <- do.call(rbind, lapply(vars, function(v) {
       data.frame(
         Parameter = v,
         PSA_Mean  = mean(df_model_psa[[v]]),
         PSA_SD    = sd(df_model_psa[[v]]),
         PSA_Min   = min(df_model_psa[[v]]),
         PSA_Max   = max(df_model_psa[[v]])
       )
     }))
     
     ## Original (non-PSA) values. Before conversion to monthly
     
     params <- f_gen_param()
     
     df_non_psa <- data.frame(
       Parameter  = c("p_CiMd", "p_MdVd", "p_VdTd", "p_VdD", "p_TdD",
                      "p_UD", "p_HD", "p_BD",
                      "c_Ci", "c_Md", "c_Vd", "c_Td",
                      "u_Ci", "u_Md", "u_Vd",
                      "du_UD"),
       Deterministic_Mean = c(params$p_CiMd,
                              params$p_MdVd,
                              params$p_VdTd,
                              params$p_VdD,
                              params$p_TdD,
                              params$p_UD,
                              params$p_HD,
                              params$p_BD,
                              params$c_Ci,
                              params$c_Md,
                              params$c_Vd,
                              params$c_Td,
                              params$u_Ci,
                              params$u_Md,
                              params$u_Vd,
                              params$du_UD))
     
     ## Merge tables
     tabulate <- merge(df_psa_summary, df_non_psa, by = "Parameter", all.x = TRUE)
     
     ## Print as a single table
     psa_dists <- kable(tabulate,
                        digits = 3,
                        caption = "See supplementary materials for visualisation of distributions",
                        col.names = c("Parameter", "PSA Mean", "PSA SD", "PSA Min", "PSA Max", "Deterministic Mean"),
                        format = "html",
                        table.attr = 'style="width: 100%;"',
                        escape = FALSE,
                        align = rep("c", ncol(tabulate))) %>%
       kable_styling(full_width = TRUE, position = "center") %>%
       column_spec(1, width = "20em", extra_css = "text-align: center;") %>%
       column_spec(2:5, width = "5em", extra_css = "text-align: center;") %>%
       column_spec(6, width = "10em", extra_css = "text-align: center;")
     
     HTML(kable_styling(psa_dists,
                        bootstrap_options = c("striped", "hover", "condensed"),
                        full_width = TRUE) %>%
            row_spec(0, bold = TRUE, extra_css = "text-align: center;") %>%
            row_spec(1:nrow(tabulate), background = ifelse(1:nrow(tabulate) %% 2 == 0, "#f9f9f9", "#ffffff")) %>%
            column_spec(1, bold = TRUE, extra_css = "text-align: center;") %>%
            column_spec(2:5, background = "#f9f9f9", extra_css = "text-align: center;") %>%
            column_spec(6, background = "#ffffff", extra_css = "text-align: center;"))
     
     
     
   }) # Closes the function for PSA table
   
   }) # Closes the withProgress function
   
 }) # Closes the observe event
  
} # Closes the server function

########################## Run the App ####################################

shinyApp(ui = ui, server = server)


