# S050_q37 processing
# Run analysis.R till line 245

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


# plot ----

# plot all data (except NTC)
ggplot(filter(processed_data, assay_variable != 'ntc'),
       aes(day, Copies_proportional, colour = assay_variable, shape = Induction, # 40-CT
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates))) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
  
  ggtitle('q37_S050_pilot2') + 
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
  
  ggtitle('q37_S050_pilot2,3') + 
  facet_wrap(facets = 'Target_name')

ggsave(plot_as('q41_S050_flipped'), plt_flip, width = 4, height = 3)
