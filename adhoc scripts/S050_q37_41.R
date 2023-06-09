# S050_q37 processing

source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- 'Fluorescent vs burdenless memory: S050'

# user inputs  ----
flnms <- c('q37_S050_pilot-2_4-10-22', 'q41_S050_pilot-3_pSS079_1-11-22')

# select the column name to use for copies -- extrapolated as 2^ (40-Cq) vs absolute quantification
column_for_copies <- quo(Copies.per.ul.template)  # or Copies_proportional for uncalibrated data

# specify targets - as expressions : required for calculate_memory_ratios()
flipped <- expr(flipped)
backbone <- expr(backbone)
chromosome <- expr(chromosome)

# translate plasmid names (assay_variable) into meaningful names
assay_variable_translation <- c('pPK' = '',
                                '144|148' = 'Silent',
                                '147|150' = 'Frugal',
                                'pSS079' = 'Fluorescent')


# processing ----
processed_data <- get_processed_datasets(flnms) %>% 
  
  # process day from sample_category
  separate(Sample_category, # separate induction and day
           into = c('Induction', 'day')) %>% 
  
  mutate(across(day, # convert day column into number
                ~ str_replace_all(.x, c('minus' = '-', 'd' = '')) %>%
                  replace_na('-1') %>% 
                  as.numeric)) %>% 
  
  # translate assay_variable for presentation plot
  mutate(across(assay_variable, ~ str_replace_all(.x, assay_variable_translation))) %>% 
  
  # change order for plotting
  mutate(across(assay_variable, ~ fct_relevel(.x, c('Fluorescent', 'Silent', 'Frugal'))))

  

# ratios ----

grouping_vars_for_ratio <- c('assay_variable', 'Induction', 'day') # to group after pivoting by target, to take medians etc.

# take ratio to backbone
ratio_data <- calculate_memory_ratios(processed_data)


# plots ----

# plot Cq : all data (except NTC)
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


# Ratio plot ----
# flip fraction

plt_flip_fraction <- {ggplot(filter(ratio_data, Induction != 'Control'), # !str_detect(assay_variable, 'ntc|MG1655')
                           aes(day, flipped_fraction, colour = assay_variable, shape = Induction, 
                               label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(assay_variable, Induction, biological_replicates),
                  alpha = if_else(str_detect(Induction, 'Induced'), 1, 0.5) # emphasize data
                  )) +
    
    scale_shape_manual(values = c(19, 1)) + # determine shapes # c(4, 19, 1)
    scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
    scale_alpha_continuous(guide = 'none', range = c(0.5, 1)) + # control line transparency
    
    # induction window
    annotate('rect', xmin = -1, ymin = 0, xmax = 0, ymax = Inf, alpha = .2) +
    
    ggtitle(title_name)} %>% 
  
  print 

format_logscale_y(plt_flip_fraction)
ggsave(plot_as('q41_S050_79+silents_flipped_fraction'), plt_flip_fraction, width = 4, height = 3)

# presentable plot
present_flip_fraction <- 
  {plt_flip_fraction + ggtitle(NULL) + # remove title
  ylab('ON state fraction of plasmid') + guides(colour = guide_legend('Designs'))} %>% print

ggsave('qPCR analysis/Archive/q41_S050_flip_fraction.png', width = 4, height = 3)
ggsave('qPCR analysis/q41_S050_flip_fraction.pdf', width = 4, height = 3)
