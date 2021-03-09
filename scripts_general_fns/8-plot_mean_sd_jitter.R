
# Plotting mean, sd and individual replicates in a jitter form
# Yet to be fully generalized function. But should work for any data frame by adding dummy columns. 
# Default parameters correspond to COVID-qPCR files

plot_mean_sd_jitter <- function(.data_list = long_processed_minimal, 
                                long_format = TRUE, 
                                
                                measure_var = 'Copy #', 
                                colour_var = Target, x_var = assay_variable, y_var = `Copy #`, 
                                facet_var = `Sample_name`, 
                                
                                sample_var = '.*', exclude_sample = F, 
                                WWTP_var = wwtp_manhole_names, exclude_WWTP = F, 
                                target_filter = '.*', 
                                
                                ascending_order = FALSE, 
                                title_text = title_name, ylabel = 'Genome copies/ul RNA', xlabel = plot_assay_variable, 
                                facet_style = 'grid')

{ # Convenient handle for repetitive plotting in the same format; Specify data format: long vs wide (specify in long_format = TRUE or FALSE)
  
  .dat_filtered <- .data_list %>% map( ~ filter(.x, 
                                                if('Sample_name' %in% colnames(.x)) str_detect(`Sample_name`, sample_var, negate = exclude_sample) else TRUE, 
                                                if('WWTP' %in% colnames(.x)) str_detect(WWTP, WWTP_var, negate = exclude_WWTP) else TRUE, 
                                                str_detect(Target, target_filter))
  )
  
  # filtering data to be plotted by user inputs
  if(long_format) # use long format if not plotting Copy #s - ex. Recovery, % recovery etc.
  {
    
    .data_to_plot <- .dat_filtered %>% map(filter,
                                           Measurement == measure_var)
    if(ascending_order) .data_to_plot$summ.dat %<>% mutate_at('WWTP', as.character) %>% 
      arrange(`mean`) %>% 
      mutate_at('WWTP', as_factor)
    
    y_var <- sym('value') # default y variable is value
    
    summ_actual_spike_in <- .dat_filtered$summ.dat %>% filter(str_detect(Measurement,'Actual'))
    
  } else
    
  {
    .data_to_plot <- .dat_filtered 
    
    if(ascending_order) .data_to_plot$summ.dat %<>% mutate_at('WWTP', as.character) %>% 
      arrange(`mean`) %>% 
      mutate_at('WWTP', as_factor)
    
  }
  
  # Exit with a useful message if data is empty
  if(.data_to_plot %>% map_lgl(plyr::empty) %>% any()) return('Data does not exist')  
  
  # plotting
  plt1 <- .data_to_plot$summ.dat %>% ggplot(aes(x = {{x_var}}, y = mean, colour = {{colour_var}})) +
    geom_point(size = 2) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) +
    
    # Individual data points
    geom_jitter(data = .data_to_plot$raw.dat, aes(y = {{y_var}}, alpha = map_chr({{y_var}}, ~. == 0), size = map_chr({{y_var}}, ~. == 0)), width = .2, show.legend = F ) +
    scale_alpha_manual(values = c(.3, 1)) + scale_size_manual(values = c(1, 2)) + # manual scale for emphasizing unamplified samples
    
    # Plotting actual spike ins (only for Recovery plot'; only with long format data )
    { if(measure_var == 'Recovered') list(geom_point(data = summ_actual_spike_in, colour = 'black', shape = 21), 
                                          geom_line(data = summ_actual_spike_in, aes(group = {{colour_var}})))
    } +
    
    # Facetting
    facet_grid(cols = vars({{facet_var}}), scales = 'free_x', space = 'free_x') +
    
    # experimental - conditional facetting (doesn't work for unknown reasons) : Just facet_grid the output to remove facets or add new!
    # { if (facet_style == 'grid') list(facet_grid(cols = vars({{facet_var}}), scales = 'free_x', space = 'free_x'))
    #   if (facet_style == 'wrap free') list(facet_wrap(facets =  vars({{facet_var}}), scales = 'free')) 
    #   else NULL
    # } +
    
    # Labelling
    ggtitle(title_text) + ylab(ylabel) #+ xlab(xlabel)
  
  plt1.formatted <- plt1 %>% format_classic() # clean formatting
  
}