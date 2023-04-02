# S067_q45_plots.R

# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 

  clean_up_water_wells() %>% 
  
  # split columns
  separate(assay_variable, into = c('plasmid', 'organism'), sep = '/', remove = FALSE) %>% 
  separate(Sample_category, into = c('day', 'AHL (uM)'), sep = '/', remove = FALSE) %>% 
  
  # change types
  mutate(across(`AHL (uM)`, 
                ~ if_else(str_detect(., 'I'), 10, 0) %>% # Induced = 10, uninduced
                  replace_na(0) %>% as.character %>% fct_inorder),   # NA is 0
         across(organism, ~ replace_na(., 'control')),
         across(day, ~ str_replace_all(., c('control' = 'd-1', '^d' = ''))))

# TODO : take ratio of flipped / backbone ; plot ratio in a timecourse ; 


# Ratios ----

# take ratio to backbone
ratio_data <- select(forplotting_cq.dat, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)



# Plotting ----


# timeseries plot from plate reader


# Small multiples-timeseries plot
timeseries <- 
  {filter(ratio_data, organism != 'control', plasmid != '143') %>% # remove empty data
      # filter(forplotting_cq.dat, Target_name == 'backbone') %>% 
      
      ggplot(aes(day, flipped_fraction, colour = plasmid, shape = `AHL (uM)`)) + 
      
      geom_point(size = 2) +
      # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
      scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
      scale_alpha_discrete(guide = 'none', range = c(0.5, 1)) + # control line transparency
      
      
      # line like a dumbell plot
      geom_line(aes(group = interaction(plasmid, `AHL (uM)`, biological_replicates), 
                    alpha = `AHL (uM)`)) +  
      
      facet_wrap(vars(organism), scale = 'free_y') +
      theme(legend.position = 'top') +
      
      # labels
      ggtitle('Memory in wastewater', subtitle = title_name)} %>% 
  
  print()

ggplotly(timeseries)
ggsave(plot_as(title_name, '-timeseries'), width = 5, height = 5)

