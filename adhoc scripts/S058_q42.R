# plotting for S058_q42

# Run general functions from analysis.R
source('./0-general_functions_main.R') # Source the general_functions file before running this

# User-inputs ----
flnm <- 'q42_wild memory-prelim_11-Jan-23' # mention the file name without the "-linreg" or "-processed" suffixes
title_name <- "q42_memory-wild"

# translate plasmid names (assay_variable) into meaningful names
assay_variable_translation <- c('pPK' = '',
                                '144|148' = 'Silent',
                                '147|150' = 'Frugal',
                                'pSS79' = 'Fluorescent')

# load data, process ----
processed_data <- 
  get_processed_datasets(flnm) %>% 
  separate(col = Sample_category, into = c('plasmid', 'organism'), remove = F) %>% 
  separate(assay_variable, into = c('day', 'AHL_uM'), sep = '/') %>% replace_na(list('AHL_uM' = '0'))

# more refining (from S050_q37_41 script)  
proc_data2 <-     
  mutate(processed_data, across(day, # convert day column into number
                ~ str_replace_all(.x, c('minus' = '-', 'd' = '')) %>%
                  replace_na('-1') %>% 
                  as.numeric)) %>% 
  
  # translate assay_variable for presentation plot
  mutate(across(plasmid, ~ str_replace_all(.x, assay_variable_translation)))

# take ratio to backbone
ratio_data <- select(proc_data2, -CT) %>% # remove the non unique columns
  filter(plasmid != 'Negative') %>% # remove ntcs
  
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)

# Plotting ----

plt_panels <-
  ggplot(filter(processed_data, Target_name == 'flipped'), # select flipped
                     aes(x = day, y = 40 - CT, colour = AHL_uM)) + 
  geom_point() + geom_line(aes(group = AHL_uM)) + 
  
  ggtitle(title_name) + 
  theme(legend.position = 'top') + 
  facet_wrap(facets = vars(Sample_category), ncol = 4, scales = 'free_x')

# save plot 
ggsave(plot_as(title_name), plt_panels, width = 4, height = 6)


# combine plot
# plot_condensed <- 
  ggplot(filter(processed_data, Target_name == 'flipped', plasmid != 'Negative'), # select flipped, no ntcs
         aes(x = day, y = 40 - CT, colour = organism, shape = AHL_uM)) + 
  geom_point() + geom_line(aes(group = interaction(AHL_uM, organism))) +
  scale_shape_manual(values = c(1, 16)) + 
  
  ggtitle(title_name) + 
  # theme(legend.position = 'top') + 
  facet_wrap(facets = vars(plasmid), scales = 'free')


# save plot 
ggsave(plot_as(title_name, '-condensed'), plot_condensed, width = 4, height = 4)



# More plotting ----

# Plot ratios -- copy number

plt_copy_number <- {ggplot(ratio_data,
                           aes(day, plasmid_copy_number, colour = plasmid, shape = AHL_uM, 
                               label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(plasmid, AHL_uM, biological_replicates))) + 
    
    scale_shape_manual(values = c(1, 19)) + # determine shapes
    scale_x_continuous(breaks = c(-1, 0, 1)) + # simplify x axis ticks
    
    theme(legend.position = 'top', legend.box = 'vertical') + 
    facet_wrap(facets = vars(organism), scales = 'free_x') +
  
    ggtitle(title_name)} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, '-copy_number'), plt_copy_number, width = 5.5, height = 6)

# plot ratio -- flip fraction

plt_flip_fraction <- {ggplot(ratio_data,
                             aes(day, flipped_fraction, colour = plasmid, shape = AHL_uM, 
                                 label = biological_replicates)) + 
    geom_point() + 
    geom_line(aes(group = interaction(plasmid, AHL_uM, biological_replicates),
                  alpha = if_else(str_detect(AHL_uM, 'Induced'), 1, 0.5) # emphasize data
    )) +
    
    scale_shape_manual(values = c(1, 19)) + # determine shapes
    scale_x_continuous(breaks = c(-1, 0, 1)) + # simplify x axis ticks
    scale_alpha_continuous(guide = 'none', range = c(0.5, 1)) + # control line transparency
    
    theme(legend.position = 'top', legend.box = 'vertical', legend.margin = margin()) + # tighter legends
    facet_wrap(facets = vars(organism), scales = 'free_x') +
    
    ggtitle(title_name)} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, '-flipped_fraction'), plt_flip_fraction, width = 4, height = 6)
