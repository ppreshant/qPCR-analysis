# S071_q48_memory_timeseries.R

# copied elements from S050_q37_41 ; subsumed: S067_q45_plots.R --> change inputs to do this


# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script

# override filename and title name from user inputs 
title_name <- 'S071_q48_pilot 1'
flnms <- c('q48a_S071_22-4-23', 'q48b_S071_23-4-23')


# For S067 use these instead of above
# flnm <- 'q45_S067_top3 till d2_29-3-23'
# title_name <- 'q45_S067_top3'


# translate plasmid names (assay_variable) into meaningful names
plasmid_translation <- c('pPK' = '',
                         '143' = 'No memory',
                         '144|148' = 'Silent',
                         '147|150' = 'Frugal',
                         'pSS079|79' = 'Fluorescent')


# Load data ----

processed_data <- get_processed_datasets(flnms) %>% # clean_up_water_wells() %>% # for q45_S067
  
  # split columns
  separate(assay_variable, into = c('plasmid', 'organism'), sep = '/', remove = FALSE) %>% 
  separate(Sample_category, into = c('day', 'AHL (uM)'), sep = '/', remove = FALSE) %>% 
  
  # translate assay_variable for presentation plot
  mutate(across(plasmid, ~ str_replace_all(.x, plasmid_translation))) %>% 

  # change types
  mutate(across(`AHL (uM)`, 
                ~ if_else(str_detect(., 'I'), 10, 0) %>% # Induced = 10, uninduced
                  replace_na(0) %>% as.character %>% fct_inorder),   # NA is 0
         across(organism, ~ replace_na(., 'control')),
         across(day, ~ str_replace_all(., c('control' = 'd-1', '^d' = ''))))

# TODO : made day numeric for proper plotting

# Ratios ----

# take ratio to backbone
ratio_data <- select(processed_data, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)

# Summary plot ----

# Small multiples-timeseries plot : Modified from plate reader

timeseries <- 
  {filter(ratio_data, organism != 'control', plasmid != 'No memory') %>% # remove empty data
      # filter(forplotting_cq.dat, Target_name == 'backbone') %>% 
      
      ggplot(aes(day, flipped_fraction, colour = plasmid, shape = `AHL (uM)`)) + 
      
      geom_point() +
      # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
      scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
      scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
      
      
      # line like a dumbell plot
      geom_line(aes(group = interaction(plasmid, `AHL (uM)`, biological_replicates), 
                    alpha = `AHL (uM)`)) +  
      
      facet_wrap(vars(organism), scale = 'free_y') +
      
      # legend
      theme(legend.position = 'top') +
      guides(shape = guide_legend(nrow = 2)) +
      
      # labels
      ggtitle('Memory in wastewater', subtitle = title_name)} %>% 
  
  print()

# TODO : once day is numeric, add ticks only where data is present (from S050_q37_41)

ggplotly(timeseries)
ggsave(plot_as(title_name, '-fraction'), width = 5, height = 5)


# individual target plot ----

plot_timeseries_target <- function(filter_target = 'flipped')
{
  filter(forplotting_cq.dat, organism != 'control', plasmid != '143', # remove empty data
         Target_name == filter_target) %>% # filter specific target
    
    ggplot(aes(day, Copies_proportional, colour = plasmid, shape = `AHL (uM)`)) + 
    
    geom_point(size = 2) +
    # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
    scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
    scale_alpha_discrete(guide = 'none', range = c(0.5, 1)) + # control line transparency
    
    # line connect data
    geom_line(aes(group = interaction(plasmid, `AHL (uM)`), 
                  alpha = `AHL (uM)`)) +  
    
    facet_wrap(vars(organism), scale = 'free_y') +
    theme(legend.position = 'top') +
    
    # labels
    ggtitle(str_c('Memory in wastewater : ', filter_target), subtitle = title_name)
}

copies_flipped <- plot_timeseries_target() %>% print
ggsave(plot_as(title_name, '-copy-flip'), copies_flipped, width = 5, height = 5)

copies_bb <- plot_timeseries_target('backbone') %>% print
ggsave(plot_as(title_name, '-copy-bb'), width = 5, height = 5)

copies_bb <- plot_timeseries_target('chromosome') %>% print
ggsave(plot_as(title_name, '-copy-chr'), width = 5, height = 5)

