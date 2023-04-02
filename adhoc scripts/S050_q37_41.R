# S050_q37 processing

source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- 'Fluorescent vs burdenless memory: S050'


flnms <- c('q37_S050_pilot-2_4-10-22', 'q41_S050_pilot-3_pSS079_1-11-22')

# translate plasmid names (assay_variable) into meaningful names
assay_variable_translation <- c('pPK' = '',
                                '144|148' = 'Silent',
                                '147|150' = 'Frugal',
                                'pSS079' = 'Fluorescent')


processed_data <- get_processed_datasets(flnms) %>% 
  
  # process day from sample_category
  separate(Sample_category, # separate induction and day
           into = c('Induction', 'day')) %>% 
  
  mutate(across(day, # convert day column into number
                ~ str_replace_all(.x, c('minus' = '-', 'd' = '')) %>%
                  replace_na('-1') %>% 
                  as.numeric)) %>% 
  
  # translate assay_variable for presentation plot
  mutate(across(assay_variable, ~ str_replace_all(.x, assay_variable_translation)))
  

# take ratio to backbone
ratio_data <- select(processed_data, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)
  


# plots ----

# plot all data (except NTC)
ggplot(filter(processed_data, assay_variable != 'ntc'),
       aes(day, 40 - CT, colour = assay_variable, shape = Induction, # 40-CT
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates),
                alpha = if_else(str_detect(Induction, 'Induced'), 1, 0.5) # emphasize data
                )) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
  scale_alpha_continuous(guide = 'none', range = c(0.5, 1)) + # control line transparency
  
  ggtitle(title_name) + 
  facet_wrap(facets = 'Target_name')

# interactive
ggplotly(dynamicTicks = T)

ggsave(plot_as('q41_S050_79+silents'), width = 5, height = 3)  


# plot only flipped data (except NTC)
plt_flip <- ggplot(filter(processed_data, assay_variable != 'ntc', 
                                  Target_name == 'flipped'),
       aes(day, 40 - CT, colour = assay_variable, shape = Induction, # 40-CT
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates))) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
  
  ggtitle(title_name) + 
  facet_wrap(facets = 'Target_name')

ggsave(plot_as('q41_S050_79+silents_flipped'), plt_flip, width = 4, height = 3)


# Plot ratios -- copy number

plt_copy_number <- {ggplot(filter(ratio_data, assay_variable != 'ntc'),
                   aes(day, plasmid_copy_number, colour = assay_variable, shape = Induction, 
                       label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates))) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
  
  ggtitle(title_name)} %>% 
  
  format_logscale_y()


# plot ratio -- flip fraction

plt_flip_fraction <- {ggplot(filter(ratio_data, assay_variable != 'ntc'),
                           aes(day, flipped_fraction, colour = assay_variable, shape = Induction, 
                               label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates),
                  alpha = if_else(str_detect(Induction, 'Induced'), 1, 0.5) # emphasize data
                  )) +
    
    scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
    scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
    scale_alpha_continuous(guide = 'none', range = c(0.5, 1)) + # control line transparency
    
    ggtitle(title_name)} %>% 
  
  print 

format_logscale_y(plt_flip_fraction)

ggsave(plot_as('q41_S050_79+silents_flipped_fraction'), plt_flip_fraction, width = 4, height = 3)

