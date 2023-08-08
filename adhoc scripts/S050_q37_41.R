# S050_q37_41 processing

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

# inducer translation : placeholders 
inducer_translation <- c('Induced' = 'I',
                         'Uninduced' = '0')

# processing ----
processed_data <- get_processed_datasets(flnms) %>% 
  
  # process day from sample_category
  separate(Sample_category, # separate Inducer and day
           into = c('Inducer', 'day')) %>% 
  
  mutate(across(day, # convert day column into number
                ~ str_replace_all(.x, c('minus' = '-', 'd' = '')) %>%
                  replace_na('-1') %>% 
                  as.numeric)) %>% 
  
  # translate assay_variable for presentation plot
  mutate(across(assay_variable, 
                ~ str_replace_all(.x, assay_variable_translation) %>% 
                  fct_relevel(c('Fluorescent', 'Silent', 'Frugal'))  # change order for plotting
                  )) %>% 
  
  # rename Inducers with shorter placeholders
  mutate(across(Inducer, ~ str_replace_all(.x, inducer_translation) %>% 
                  fct_relevel('0', 'I') # reorder
                  ))
  

# ratios ----

grouping_vars_for_ratio <- c('assay_variable', 'Inducer', 'day') # to group after pivoting by target, to take medians etc.

# take ratio to backbone
source('scripts_general_fns/22-memory_wrappers_ratio.R')
ratio_data <- calculate_memory_ratios(processed_data)


# Cq plots ----

# plot Cq : all data (except NTC)
ggplot(filter(processed_data, assay_variable != 'ntc'),
       aes(day, 40 - CT, colour = assay_variable, shape = Inducer, # 40-CT
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Inducer, biological_replicates),
                alpha = if_else(str_detect(Inducer, 'Induced'), 1, 0.5) # emphasize data
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
       aes(day, 40 - CT, colour = assay_variable, shape = Inducer, # 40-CT
           label = biological_replicates)) + 
  geom_point() + 
  geom_line(aes(group = interaction(assay_variable, Inducer, biological_replicates))) + 
  
  scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
  scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
  
  ggtitle(title_name) + 
  facet_wrap(facets = 'Target_name')

ggsave(plot_as('q41_S050_79+silents_flipped'), plt_flip, width = 4, height = 3)


# Copies plots ----

# plotting Copies of all targets -- potential S10 fig

plt_copies_all <- 
  {ggplot(
    mutate(processed_data, across(Target_name, ~ fct_relevel(.x, c('flipped', 'backbone', 'chromosome')))),
    
          #filter(processed_data, Inducer != 'Control'), # !str_detect(assay_variable, 'ntc|MG1655')
          
          aes(day, !!column_for_copies, colour = assay_variable, shape = Inducer, 
              label = biological_replicates)) + 
      geom_point() + 
      geom_line(aes(group = interaction(assay_variable, Inducer, biological_replicates),
                    # alpha = if_else(str_detect(Inducer, 'Induced'), 1, 0.5) # emphasize data
                    alpha = Inducer
      )) +
      
      # Aesthetics
      scale_shape_manual(values = c(1, 16, 6)) + # determine shapes # c(4, 19, 1)
      scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
      scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5, 0.2)) + # control line transparency
      
      # induction window
      annotate('rect', xmin = -1, ymin = 0, xmax = 0, ymax = Inf, alpha = .2) +
      
      # Layout : legend position
      theme(legend.position = 'top', 
            legend.box = 'vertical', # Split into 2 lines
            legend.spacing = unit(0, 'mm'), # minimize spacing
            legend.box.just = 'left') + # left justified
      
      facet_grid(rows = vars(Target_name), scales = 'free') +
      ggtitle(title_name)} %>% 
  
  print 

ggsave(plot_as('q41,37_S050-all copies'), width = 6, height = 8)

# logscale - different orientation
plt_copies_log <- 
  {plt_copies_all + facet_grid(cols = vars(Target_name)) + ggtitle(NULL)} %>% format_logscale_y()

ggsave(plot_as('q41,37_S050-log copies'), plt_copies_log, width = 6, height = 5)
ggsave('qPCR analysis/Archive/q41,37_S050-log copies.pdf', plt_copies_log, width = 6, height = 5)


# Ratio plot ----

# flipped fraction

plt_flip_fraction <- {ggplot(filter(ratio_data, Inducer != 'Control'), # !str_detect(assay_variable, 'ntc|MG1655')
                           aes(day, flipped_fraction, colour = assay_variable, shape = Inducer, 
                               label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(assay_variable, Inducer, biological_replicates),
                  # alpha = if_else(str_detect(Inducer, 'Induced'), 1, 0.5) # emphasize data
                  alpha = Inducer
                  )) +
    
    # Aesthetics
    scale_shape_manual(values = c(1, 16)) + # determine shapes # c(4, 19, 1)
    scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
    scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
    
    # induction window
    annotate('rect', xmin = -1, ymin = 0, xmax = 0, ymax = Inf, alpha = .2) +
    
    # Layout : legend position
    theme(legend.position = 'top', 
          legend.box = 'vertical', # Split into 2 lines
          legend.spacing = unit(0, 'mm'), # minimize spacing
          legend.box.just = 'left') + # left justified
    ggtitle(title_name)} %>% 
  
  print 

format_logscale_y(plt_flip_fraction)
ggsave(plot_as('q41_S050_79+silents_flipped_fraction'), plt_flip_fraction, width = 4, height = 3)


# presentable plot

# colours for memory picked by SS : coolors.co/
SS_colourscheme <- c('#9E2A2B', '#73956F', '#003844', '#F3C677', '#003844') # 'red', light green, light yellow, blue

present_flip_fraction <- 
  {plt_flip_fraction + ggtitle(NULL) + # remove title
      ylab('ON state fraction of plasmid') + guides(colour = "none") + # guide_legend('Designs')
      scale_color_manual(values = SS_colourscheme)} %>% print

ggsave('qPCR analysis/Archive/q41_S050_flip_fraction.png', width = 4, height = 5)
ggsave('qPCR analysis/q41_S050_flip_fraction-v2.1.pdf', width = 3, height = 3)


# Calculations / stats ---- 

# show the fold changes and % differences from d1 to d8
end_points <- filter(ratio_data, 
                     day == 1 | day == 8, Inducer =='I') %>% # select end points / induced only
  select(median_flipped_fraction) %>% unique %>% # only the summary value
  
  pivot_wider(names_from = day, values_from = median_flipped_fraction, # move each day to a column
              names_prefix = 'frac_d') %>% 
  mutate(differ_from_d1 = (frac_d1 - frac_d8)/frac_d1, # percent change from d1 to d8
         fold_from_d1 = frac_d1 / frac_d8) %>%  # fold change from d1 to d8
  print
# Conclusion : Silent and frugal change by < .2 fold or 30%


# T tests

# make nested data
condensed_data <- ratio_data %>% 
  group_by(assay_variable, Inducer) %>% 
  filter(str_detect(day, '^1|8')) %>% # select only d1 (first day measured for sil, fli) and d8
  nest() %>% 
  view


stat_data <- 
  mutate(condensed_data, 
    ttest = map(data, 
                ~ t.test(flipped_fraction ~ day, alternative = 'greater', paired = T,
                         data = .x)),
    
    pval = map(ttest, ~ .x$p.value)
    
    ) %>% view

stat_data$ttest[[6]]

# t.test for d1 d8 differences? -- paired data
# t.test(frac_d1, frac_d8, data = ..)

# testing subset
cs <- condensed_data$data[[6]]
t.test(flipped_fraction ~ day, alternative = 'greater', paired = T, data = cs)
# plasmid copy number ----

# Plot ratios -- copy number

plt_copy_number <- {ggplot(filter(ratio_data, assay_variable != 'ntc'),
                           aes(day, plasmid_copy_number, colour = assay_variable, shape = Inducer, 
                               label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(assay_variable, Inducer, biological_replicates))) + 
    
    scale_shape_manual(values = c(4, 19, 1)) + # determine shapes
    scale_x_continuous(breaks = c(-1, 0, 1, 7, 8)) + # simplify x axis ticks
    
    ggtitle(title_name)} %>% print

format_logscale_y(plt_copy_number)

