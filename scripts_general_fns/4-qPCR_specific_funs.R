# qPCR related functions (good to copy to the general qPCR project)

# qPCR related  ----

# lookup table of primer pairs and respective targets - not useful as same information is in `Target` already
primer_table <- c('q1-3' = 'Flipped', 'q4-5' = 'Flipped', 
                  'q5-11' = 'Unflipped', 'q1-2' = 'Unflipped', 'q9-10' = 'Unflipped', 'q4-2' = 'Unflipped', 'q2-4' = 'Unflipped', 'q5-11' = 'Unflipped', 'q6-7' = 'Unflipped',
                  'q12-13' = 'Backbone')



#' Calculate absolute copies with qPCR data and standard curve parameters
calculate_absolute_quant <- function(.df = forplotting_cq.dat, std_parameters = std_par)
{
  .df %>%
    select(-Copies_proportional) %>% # remove the dummy copies data
    
    group_by(Target_name) %>%
    nest() %>% # create a new column with data frames for each target
    
    # calculate copy number for each dataset, and the mean for replicates
    summarize(w.copy.data = map2(data, Target_name,  
                                 ~ absolute_calculation_within_target(.x, .y, std_parameters) ) 
    ) %>% 
    unnest(cols = c(w.copy.data)) %>%  # expand the absolute copy number data list
    
    # append the ID of the std curve and equation used
    mutate('std_curve id' = std_to_retrieve)
  
}


# function to back-calculate CT using standard curve parameters
absolute_calculation_within_target <- function(.df, .target_current, std_par)
{
  std_current <- std_par %>% filter(Target_name == .target_current)
  
  if(!plyr::empty(std_current)) { # if data for the current target exists
    # Back-calculate the unknown Copies from the CT using the linear standard curve
    mutate(.df, `Copies.per.ul.template` = 10^( (CT - std_current$y_intercept)/std_current$Slope) ) %>% 
      
      group_by(Sample_category, assay_variable) %>%  # grouping by everything except replicates
      mutate(mean_Copies.per.ul.template = mean(Copies.per.ul.template, na.rm = TRUE)) %>% # find mean of the replicates
      
      # record the std curve equation as text
      mutate('std_curve equation' = glue::glue('Cq = {std_current$Slope} x log10(copies) + {std_current$y_intercept}') )

    
  } else .df
}


# Plot Standard curve
plotstdcurve <- function(results_qpcr, plttitle, xlabel)
{
  plt <- results_qpcr %>% 
    ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target_name`) + geom_point() +
    theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle(plttitle, subtitle = expression(C[q])) + xlab(xlabel) + ylab('') +
    stat_smooth(data = filter(results_qpcr, str_detect(Sample_category, 'Std')), method ="lm", se = F) # plots linear regression line
}


# output the difference between consecutive CT values
optional1 <- function()
{
  tsumrev <- trev %>% group_by(`Sample_name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
  diff(tsumrev$CT) %>% round(2)
}


# Tm plots ----


# plotting functions for Melting temperature
# Obsolete: Were used in small_scale mode and old_assay mode

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample_name`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample_name`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample_name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample_name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

