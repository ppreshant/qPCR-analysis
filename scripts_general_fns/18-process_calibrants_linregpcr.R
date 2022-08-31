# 18-process_calibrants_linregpcr.R

# Identify Stds from linregpcr processing files and process and save data + plot after user approval

process_calibrants_linregpcr <- function(forplotting_cq.dat)
{
  
  # Extract a short name for the calibrant standard set from the file name provided. short name has qxx
  fl_namer <- str_match(base_title_name,'^q[:alnum:]*') %>% 
    str_c(., '_linregStd')
  
  calibrant.data <- 
    filter(forplotting_cq.dat, str_detect(Sample_category, 'Std')) %>% # select only standards
    select(-contains('var.')) %>%  # remove irrelevant label columns
    arrange(Target_name) %>% # arrange by target to keep replicates together
    
    mutate(std_concentration = as.numeric(assay_variable), .after = 1) %>% # make the concentrations numeric
    mutate(N0.per.copy = N0/std_concentration, .after = 1) %>%  # normalize the standards
    mutate(mean_N0.per.copy = mean_N0/std_concentration, .after = 1) %>% # normalize the standards - avg
    
    mutate(ID = str_replace(fl_namer, '_.*', ''), .before = 1) # add the qxx ID to the dataset
  
  
  # Plotting for visual validation
  
  calibrant_plt <- ggplot(calibrant.data, aes(x = Target_name, y = N0.per.copy)) + geom_point() +  # visual validation
    geom_point(aes(y = mean_N0.per.copy), size = 8, shape = '-', colour = 'red') + # show the mean
    ggtitle(str_c(fl_namer, ': N0/conc.'))
  log_cal_plt <- format_logscale_y(calibrant_plt) + ggtitle('Logscale') # logscale plot
  
  patchwork::wrap_plots(list(calibrant_plt, log_cal_plt), nrow = 2) %>% print() # combine plots
  
  
  # Let the user approve the calibrants (in case some standards need to be excluded/ --AddFeature)
  # --AddFeature: excluding targets directly from user input -- might need alternative to 'menu' function
  
  proceed_with_standards <- menu(c('Yes', 'No'), title = paste("Check the calibrant standards plots:", 
                                                               fl_namer, 
                                                               "on the right side in Rstudio. 
                                                               
   Do you wish to continue with saving the calibrant standards? 
   Select NO if you wish to remove any targets (using '?? option') or change something else and re-run the script", sep=" "))
  
  if (proceed_with_standards == 2){ # if standards are rejected -- 
    
    # save rejected standards into separate folder 
    ggsave(str_c('qPCR analysis/Standards/bad standards/', fl_namer , '.png'), width = 7, height = 6)
    
    # Stop workflow with an error message
    stop("Cancel selected, script aborted. Standards not saved.")
  }
  
  # If calibrants are accepted ----
  
  # Save plot
  ggsave(str_c('qPCR analysis/Standards/', fl_namer , '.png'), width = 7, height = 6)
  
  # Data output ----
  
  # Select the mean calibration data 
  mean_calibrant.data <- ungroup(calibrant.data) %>% select(1,2,3,5, mean_N0) %>% unique()
  
  # Writing data ----
  
  # write mean data
  # sheet_append(sheeturls$plate_layouts_PK, mean_calibrant.data, sheet = 'Linregpcr_calibrants') # Add parameters to standard curves into google sheet
  write_csv(mean_calibrant.data, 'qPCR analysis/Standards/linregpcr_calibrants.csv', append = TRUE)
  
  # Raw data of all replicates
  write_csv(calibrant.data, 'qPCR analysis/Standards/linregpcr_all_stds.csv', append = TRUE) # Writing full std data for archive
  
  
}