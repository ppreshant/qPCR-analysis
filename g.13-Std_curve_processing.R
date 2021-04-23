# Extracts stds from a plate and calculates std curves for various targets and saves it to google sheet
# Prashant; 23/April/21
# Adapted from COVID qPCR/WW_master branch

process_standard_curve <- function(flnm)
{ # note: Quantity of copies/well must be written in the template sheet for the standards
  
  # Preliminary naming ----
  
  # file path and google sheet urls
  flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
  
  # Extract a short name for the standard curve from the file name provide. short name has Stdxx_Target_WWxx
  fl_namer <- c('Std[:alnum:]*', 'BCoV|N1|N2|BRSV|pMMoV', 'WW[:alnum:]*') %>% 
    map_chr(~str_match(flnm, .)) %>% 
    str_c(collapse = '_')
  
  if(is.na(fl_namer)) 
  {
    stop(str_c('Filename:', flnm,  ' - input to the standard curve function, is either missing the standard curve ID (ex: Std25), the WW ID (Ex: WW61) or the Target name (Ex: BCoV) \n
         Check README for proper naming convention: https://github.com/ppreshant/WW-CoV2-project/'))
  }
  
  title_name <- str_c('Standard curve: ', fl_namer ,' - Fastvirus 4x') # title name for plots
  
  
  # Data input ----
  
  fl <- readqpcr(flpath) # read file
  
  # Bring sample names from template google sheet
  plate_template <- get_template_for(flnm, sheeturls$templates)
  
  # this gives a vector to order the samples columnwise in the PCR plate or strip 
  # (by default : data is shown row-wise) => This command will enable plotting column wise order
  sample_order = columnwise_index(fl) 
  
  bring_results <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`, Task) %>% rename(Target = `Target Name`) %>%  .[sample_order,] %>%  
    
    # select only the results used for plotting, calculations etc. and arrange them according to sample order
    select(-`Sample Name`) %>% right_join(plate_template, by = 'Well Position') %>%  # Incorporate samples names from the google sheet by matching well position
    separate(Sample_name, c(NA, 'Category', 'Quantity'), sep = '-|_') %>% mutate_at('Quantity', ~ replace_na(as.numeric(.), 0)) %>% 
    filter(!is.na(Target))
  
  # optional filtering to remove low concentration points in standard curve
  bring_results %<>% filter(Quantity > 1| Quantity == 0) # filtering only standard curve within the linear range
  
  # plotting ----
  
  plt <- bring_results %>% filter(str_detect(Category, 'NTC|Std')) %>%  plotstdcurve(title_name, 'log(Copy #)') # plot standard curve
  
  # # Extract the names of the targets in use
  # targets_used <- fl$Results %>% filter(Task == 'STANDARD') %>% pull(`Target Name`) %>% unique(.)  
  
  # Isolating standard curve variables (Quantity,CT) of the different targets into groups
  standard_curve_vars <- bring_results %>% filter(Task == 'STANDARD')  %>% select(Quantity, CT, Target) %>% group_by(Target) # select required columns and group
  
  # Apply linear regression and find the model fitting results (equation and slope, R2 values) for each target
  std_table <- standard_curve_vars %>% do(., equation = lm_std_curve(.), params = lm_std_curve(., trig = 'coeff'), dat = .[1,] ) # "do" applies functions to each group of the data
  std_table$params %<>% bind_rows() # Convert parameters and data into tibbles : "do" function makes htem lists
  std_table$dat %<>% bind_rows()  
  
  std_table$dat$CT <- max(standard_curve_vars$CT, na.rm = T) - 2 * seq_along(std_table$Target) + 2 # manual numbering for neat labelling with geom_text
  
  # Add labels to plot - linear regression equation
  plt.with.eqn <- plt + geom_text(data = std_table$dat, label = std_table$equation, parse = TRUE, show.legend = F, hjust = 'inward', nudge_x = 0, force = 10)
  print(plt.with.eqn)
  
  # Let the user approve the plot (in case some standards need to be excluded/ incorrect standards concentration order)
  proceed_with_standards <- menu(c('Yes', 'No'), title = paste("Check the standard curve plot:", 
                                                               fl_namer, 
                                                               "on the right side in Rstudio. 
   Do you wish to continue with saving the standard curve parameters? Select NO if you wish to change something and re-run the script", sep=" "))
  
  if (proceed_with_standards == 2){
    stop("Cancel selected, script aborted.")
  }
  
  
  # Save plot
  ggsave(str_c('qPCR analysis/Standards/', fl_namer , '.png'), width = 5, height = 4)
  
  # Data output ----
  
  
  # processing linear regression out
  efficiency_table <- tibble(Slope = std_table$params %>% 
                               pull(slope), y_intercept = std_table$params %>% 
                               pull(y_intercept) , Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = std_table$params %>% 
                               pull(r_square) %>% round(3)
  ) %>% 
    mutate(Target = std_table$dat$`Target`) %>% 
    select(Target, everything()) %>% 
    mutate(ID = fl_namer, .before = 1)
  
  # Writing data
  sheet_append(sheeturls$data_dump, efficiency_table, sheet = 'Standard curves') # Add parameters to standard curves in data dump google sheet
  
}