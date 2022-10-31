# S050_q37 processing
# Run analysis.R till line 245

processed_data <- 
  separate(polished_cq.dat, Sample_category, # separate induction and day
           into = c('Induction', 'day')) %>% 
  
  mutate(across(day, # convert day column into number
                ~ str_replace_all(.x, c('minus' = '-', 'd' = '')) %>%
                  replace_na('-1') %>% 
                  as.numeric))

ggplot(filter(processed_data, assay_variable != 'ntc'),
       aes(day, 40-CT, colour = assay_variable, shape = Induction, 
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates))) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 1, 8)) + # simplify x axis ticks
  
  ggtitle('q37_S050_pilot2') + 
  facet_wrap(facets = 'Target_name')

ggplotly(dynamicTicks = T)

ggsave(plot_as('q37_S050_pilot2'), width = 5, height = 3)  
