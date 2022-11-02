# Extracts stds from a plate and calculates std curves for various targets and saves it to google sheet
# Prashant; 23/April/21
# Adapted from COVID qPCR/WW_master branch


# Trigger into this function
  # if it is a standard curve holding file (Stdx), call standard curve processor
  # if(str_detect(flname, 'Std[:digit:]*')) process_standard_curve(flname)



process_standard_curve <- function(flnm, .dat_pol, dilutions_to_truncate = 0)
{ # note: Quantity of copies/well must be written in the template sheet for the standards
  
  # Preliminary naming ----
  
  # Extract a short name for the standard curve from the file name provide. short name has Stdxx_qxx
  fl_namer <- c('Std[:alnum:]*', '^q[:alnum:]*') %>% 
    map_chr(~str_match(flnm, .)) %>% 
    str_c(collapse = '_')
  
  if(is.na(fl_namer)) 
  {
    stop(str_c('Filename:', flnm,  ' - input to the standard curve function, is either missing the standard curve ID (ex: Std25), the run ID (Ex: q06) \n
         Example plate name : q07_16sRAM RNAtest-5_Std7'))
  }
  
  title_name <- str_c('Standard curve: ', fl_namer) # title name for plots
  
  
  # Data input ----
  
  std_results <- .dat_pol %>% 
    
    filter(str_detect(Sample_category, 'Std') | 
             str_detect(assay_variable, 'NTC|ntc')) %>% # only retain standards or NTCs
    mutate(Quantity = str_replace(assay_variable, 'NTC|ntc', '0') %>% as.numeric) %>%  # replacing NTC with 0 and make numeric
    
    filter(Target_name %in% # only retain the targets which have standards
             (filter(., 
                     str_detect(Sample_category, 'Std')) %>% # pulls only the targets which have standards
                pull(Target_name)) 
           ) %>% 
    # optional filtering to remove low concentration points in standard curve
    filter(Quantity > 1| Quantity == 0) %>%  # filtering only standard curve with realistic quantity
    
    # # Remove the bottom x dilutions by user input : to fine tune the std curve
    remove_last_n_dilutions(dilutions_to_truncate = dilutions_to_truncate)
  
  # plotting ----
  
  plt <- std_results %>%
    plotstdcurve(title_name, 'log(Copies/ul template)') # plot standard curve
  
  # # Extract the names of the targets in use
  # targets_used <- fl$Results %>% filter(Task == 'STANDARD') %>% pull(`Target_name`) %>% unique(.)  
  
  # Isolating standard curve variables (Quantity,CT) of the different targets into groups
  standard_curve_vars <- std_results %>% 
    filter(str_detect(Sample_category, 'Std') ) %>% # Only standards, NTCs are removed : old: Task == 'STANDARD' 
    select(Quantity, CT, Target_name) %>% 
    group_by(Target_name) # select required columns and group
  
  # Apply linear regression and find the model fitting results (equation and slope, R2 values) for each target
  std_table <- standard_curve_vars %>% 
    do(., equation = lm_std_curve(.), 
       params = lm_std_curve(., trig = 'coeff'), 
       dat = .[1,] ) # "do" applies functions to each group of the data
  
  # unnest might work here?
  std_table$params %<>% bind_rows() # Convert parameters and data into tibbles : "do" function makes htem lists
  std_table$dat %<>% bind_rows()  
  
  std_table$dat$CT <- max(standard_curve_vars$CT, na.rm = T) - 2 * seq_along(std_table$Target_name) + 2 # manual numbering for neat labelling with geom_text
  
  # Add labels to plot - linear regression equation
  plt.with.eqn <- plt + geom_text(data = std_table$dat, label = std_table$equation, parse = TRUE, show.legend = F, hjust = 'inward', nudge_x = 0, force = 10)
  print(plt.with.eqn)
  
  # Let the user approve the plot (in case some standards need to be excluded/ incorrect standards concentration order)
  proceed_with_standards <- menu(c('Yes', 'No'), title = paste("Check the standard curve plot:", 
                                                               fl_namer, 
                                                               "on the right side in Rstudio. 
                                                               
   Do you wish to continue with saving the standard curve parameters? 
   Select NO if you wish to truncate any dilutions (using 'dilutions_to_truncate') or change something else and re-run the script", sep=" "))
  
  if (proceed_with_standards == 2){ # if standards are rejected -- 
    
    # save rejected standards into separate folder 
    ggsave(str_c('qPCR analysis/Standards/bad standards/', fl_namer , '.png'), width = 5, height = 4)
    
    # Stop workflow with an error message
    stop("Cancel selected, script aborted. Standards not saved.")
  }
  
  
  # Save plot
  ggsave(str_c('qPCR analysis/Standards/', fl_namer , '.png'), width = 5, height = 4)
  
  # Data output ----
  
  
  # processing linear regression output
  efficiency_table <- tibble(Slope = std_table$params %>% 
                               pull(slope), y_intercept = std_table$params %>% 
                               pull(y_intercept) , Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = std_table$params %>% 
                               pull(r_square) %>% round(3)
  ) %>% 
    mutate(Target_name = std_table$dat$`Target_name`) %>% 
    select(Target_name, everything()) %>% # bring Target_name to the first column position 
    mutate(ID = fl_namer, .before = 1) %>%  # add the run ID
    mutate('master mix' = '') # add an empty entry for master mix -- record manually
  
  # formatting std curve points for writing (archive for meta-analysis of freeze thaw effects)
  clean_std_data <- std_results %>% # reorder columns 
    select(Target_name, CT, Quantity, biological_replicates, Sample_category, everything()) %>% 
    mutate(ID = fl_namer, .before = 1) # add the run ID
  
  # Writing data : google sheet or excel file 
  
  if(template_source == 'googlesheet') # Add parameters to standard curves
    sheet_append(sheeturls$plate_layouts_PK, efficiency_table, sheet = 'qPCR Std curves') else { # into google sheet
      write_csv(efficiency_table, path = 'qPCR analysis/Standards/qPCR_Std_curve_parameters.csv', append = TRUE)} # into excel file

  write_csv(clean_std_data, 'qPCR analysis/Standards/all_std_data_points.csv',append = TRUE) # Writing full std data for archive
  
  # return the table with std curve parameters
  return(efficiency_table)
}