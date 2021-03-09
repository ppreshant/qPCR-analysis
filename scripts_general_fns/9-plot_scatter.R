
# Scatter plot with a linear regression fit and equation
# Almost generalized: Only need to make it work with grouping_var - accounting for NULL and presence of a variable...

plot_scatter <- function(.data = processed_quant_data, 
                         x_var = N1_multiplex, y_var = N2_multiplex, colour_var = NULL, shape_var = NULL,
                         grouping_var = NULL, # CURRENTLY ONLY WORKS FOR NULL GROUP
                         already_pivoted_data = 'no',
                         title_text = title_name,
                         
                         show_y.equals.x_line = 'yes',
                         text_for_equation = 'Rsquare', # choice: "Rsquare" or "full equation"
                         measure_var = 'Copy #',
                         text_cols = minimal_label_columns,
                         sample_var = str_c(extra_categories, '|NTC|Vaccine'), 
                         exclude_sample = T)

{ # Convenient handle for repetitive plotting in the same format; 
  
  
  # Preliminaries
  
  
  # convert column names (target names) into string
  checkx <- strx <- paste(substitute(x_var)) 
  checky <- stry <- paste(substitute(y_var))
  
  modx <- function(x) x # unimplemented feature: modifier in case this function gets transformed variables
  
  if(already_pivoted_data == 'no')
  {  # filtering data for plotting according to function inputs
    .data_for_plot <- .data %>% 
      select(all_of(text_cols), biological_replicates, all_of(measure_var)) %>% 
      filter(str_detect(`Sample_name`, sample_var, negate = exclude_sample)) %>% # account for missing var
      pivot_wider(names_from = 'Target', values_from = all_of(measure_var)) %>% # Target can be generalized?
      ungroup() # why did you ungroup - for the lm ..?
  } else .data_for_plot <- .data #%>%  # direct carrying of data to next steps
  # {if(!is.null(grouping_var)) group_by(., {{grouping_var}}) else . } # DISABLED FOR CHECKING IF PLOTLY RUNS (GROUPS NOT WORKING RIGHT NOW
  
  
  
  # error handling
  
  # If plotting transformations of variables : Not fully implemented yet
  if(enexpr(x_var) %>% is.call()) {checkx <-  enexpr(x_var)[2] %>% paste(); }# modx <- eval(enexpr(x_var)[1])}
  if(enexpr(y_var) %>% is.call()) checky <-  enexpr(y_var)[2] %>% paste()
  
  # If X/Y vars are not present in the data, stop running and print a message
  if(.data_for_plot %>% names() %>% 
     {checkx %in% . & checky %in% . } %>% 
     !.) return('X or Y axis values for this scatterplot are not present in the data frame provided. 
Check if x_var and y_var are present in .data')
  
  # If duplicates of data exist in the data, stop running and print a message
  if((.data_for_plot %>% select(all_of(c(checkx, checky))) %>% 
      map_lgl(~class(.) == 'numeric') %>% 
      sum()) < 2) 
  { duplicated_data_points <- .data_for_plot %>% 
    filter(map({{x_var}}, length) > 1)
  print(duplicated_data_points)
  
  return('Repeated data instances for the same WWTP found in this scatterplot')
  }
  
  
  
  # For linear regression data
  
  # Making linear regression formula (source: https://stackoverflow.com/a/50054285/9049673)
  fmla <- new_formula(enexpr(y_var), enexpr(x_var))
  
  # Max and ranges for plotting
  xyeq <- .data_for_plot %>% summarise(across(where(is.numeric), max, na.rm = T)) %>% select(all_of(c(checkx, checky))) %>% min() %>% {.*0.9} %>% modx()
  
  # linear regression equation
  lin_reg_eqn <- .data_for_plot %>% mutate(across(all_of(c(checkx, checky)), ~if_else(.x == 0, NaN, .x))) %>% 
    lm(fmla, data = ., na.action = na.exclude) %>% lm_eqn(., trig = text_for_equation)
  
  
  
  
  # plotting part
  
  plt1 <- .data_for_plot %>% 
    ggplot(aes(x = {{x_var}}, y =  {{y_var}} )) +
    geom_point(size = 2, mapping = aes(colour = {{colour_var}}, shape = {{shape_var}})) +
    
    # linear regression
    geom_smooth(method = 'lm') + # .. , mapping = aes(group = {{grouping_var}})  # DISABLED GROUPS 
    geom_text(data = . %>% summarise(across(where(is.numeric), max, na.rm = T) ),
              # mapping = aes(group = {{grouping_var}}), # DISABLED FOR CHECKING IF PLOTLY RUNS (GROUPS NOT WORKING RIGHT NOW)
              label = lin_reg_eqn, parse = TRUE, show.legend = F, hjust = 'inward', nudge_x = -5) +
    
    # Dummy y = x line
    {if(show_y.equals.x_line == 'yes') list(geom_abline(slope = 1, intercept = 0, alpha = .4),
                                            annotate(geom = 'text', x = xyeq, y = xyeq, label = 'y = x', alpha = .3))
    } +
    
    # Labeling
    ggtitle(title_text, subtitle = measure_var)
}

