# S071_q48_memory_timeseries.R

# copied elements from S050_q37_41 ; subsumed: S067_q45_plots.R --> change inputs to do this


# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script

# override filename and title name from user inputs 
title_name <- 'S071_q48_pilot 1'
flnms <- c('q48a_S071_22-4-23', 
           'q48b_S071_23-4-23',
           'q48c_S071_d2_26-4-23',
           'q48d_S071_d8_26-4-23')


# For S067 use these instead of above
# flnm <- 'q45_S067_top3 till d2_29-3-23'
# title_name <- 'q45_S067_top3'


# translate plasmid names (assay_variable) into meaningful names
plasmid_translation <- c('pPK' = '',
                         '143' = 'No memory',
                         '144|148' = 'Silent',
                         '147|150' = 'Frugal',
                         'pSS079|79' = 'Fluorescent')

metadata_columns <- c("plasmid", "organism", "day", "AHL (uM)")

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
         across(day, ~ str_replace_all(., c('control' = 'd-1', '^d' = ''))),
         across(day, as.numeric)) # convert day to numeric

# take means for plotting : something wrong here, NAs not being ignored
processed_mean <- 
  reframe(processed_data,
          .by = all_of(c(metadata_columns, 'Target_name')),
          across(where(is.numeric), ~ mean(.x, na.rm = T)))

processed_mean_toplt <- filter(processed_mean, organism != 'control', plasmid != 'No memory')

# Ratios ----

# take ratio to backbone
ratio_data <- select(processed_data, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)


# mean data 

ratio_mean <- select(ratio_data, -biological_replicates) %>% 
  reframe(.by = all_of(metadata_columns),
          across(where(is.numeric), ~ mean(.x, na.rm = T)))

ratio_mean_toplt <- filter(ratio_mean, organism != 'control', plasmid != 'No memory')


# Normalizations ----

# divide by mean on d-1..

# verified that flipped and bb are not NaNs on d-1
# processed_mean %>% filter(is.na(Copies_proportional), day == -1) %>% view

normalized_data <- 
  processed_data %>% 
  group_by(across(all_of(c("plasmid", "organism", 'Target_name')))) %>% # all except days, AHL
  
  nest() %>% # nest all the other variables : for doing normalization and fitting exponential
  
  # bring out the mean and time 0
  mutate(initial_mean_Copies =  # get copies of no AHL at d-1 to normalize by
           map_dbl(data, 
                   ~ filter(., day == -1, `AHL (uM)` == 0) %>% {mean(.$Copies_proportional, na.rm = T)}),
         
         # Normalization : within each data frame in the nest
         data = map(data, 
                    ~ mutate(., normalized_Copies = 
                               .$Copies_proportional / initial_mean_Copies) %>%  # divide such that mean starts at 1
                      
                      group_by(day, `AHL (uM)`) %>% # group to calculate mean
                      
                      # calculating mean of the normalized
                      mutate(mean_normalized_Copies = mean(normalized_Copies, na.rm = TRUE))
         ) 
         
  ) %>%  
  
  unnest(data) # show full data 
  
  

# Summary plot ----

# Small multiples-timeseries plot : Modified from plate reader
xjitter <- position_jitter(width = 0.2, height = 0, seed = 1)

timeseries <- 
  {filter(ratio_data, organism != 'control', plasmid != 'No memory') %>% # remove empty data
      # filter(`AHL (uM)` == 0) %>% # filter uninduced only
      # filter(forplotting_cq.dat, Target_name == 'backbone') %>% 
      
      ggplot(aes(day, flipped_fraction, colour = plasmid, shape = `AHL (uM)`)) + 
      
      geom_point() + # points with jitter # position = xjitter
      scale_x_continuous(breaks = c(-1, 0, 2, 5, 8)) + # simplify x axis ticks
      scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
      scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
      
      
      # line join means / replicates
      # geom_line(aes(group = interaction(plasmid, `AHL (uM)`, biological_replicates), # join replicates
      #               alpha = `AHL (uM)`)) +
      geom_line(aes(alpha = `AHL (uM)`),
                data = ratio_mean_toplt) + # join means
      
      # facet_wrap(vars(organism), scale = 'free_y') +
      facet_grid(vars(organism), vars(plasmid), scale = 'free_y') +
      
      # legend
      theme(legend.position = 'top') +
      # guides(shape = guide_legend(nrow = 2)) + # split into 2 rows
      
      # labels
      ggtitle('Memory in wastewater', subtitle = title_name)} %>% 
  
  print()

ggplotly(timeseries)
ggsave(plot_as(title_name, '-fraction-facets'), width = 7, height = 5)
# ggsave(plot_as(title_name, '-mean-fraction'), width = 5, height = 5)
# ggsave(plot_as(title_name, '-fraction-uninduced'), width = 5, height = 5)



# individual target plot ----

plot_timeseries_target <- function(filter_target = 'flipped', .connect = 'mean',
                                   .data = processed_data, .yvar = Copies_proportional)
{
  filter(.data, organism != 'control', plasmid != 'No memory', # remove empty data
         Target_name == filter_target) %>% # filter specific target
    
    ggplot(aes(day, {{.yvar}}, colour = plasmid, shape = `AHL (uM)`)) + 
    
    geom_point(size = 2) +
    # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
    scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
    scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
    scale_x_continuous(breaks = c(-1, 0, 2, 5, 8)) + # simplify x axis ticks
    
    # line connect data / means
    {if(.connect == 'mean')
    {geom_line(aes(alpha = `AHL (uM)`),
               data = filter(processed_mean_toplt, Target_name == filter_target)) # join means
    } else 
    {geom_line(aes(group = interaction(plasmid, `AHL (uM)`, biological_replicates),
                   alpha = `AHL (uM)`)) }
      } +

    
    
    # facet_wrap(vars(organism), scale = 'free_y') +
    facet_grid(vars(organism), vars(plasmid), scale = 'free') +
    
    theme(legend.position = 'top') +
    
    # labels
    ggtitle(str_c('Memory in wastewater : ', filter_target), subtitle = title_name)
}

copies_flipped <- plot_timeseries_target(.connect = 'ind') %>% print
ggsave(plot_as(title_name, '-copy-flip'), copies_flipped, width = 7, height = 5)

# join means 
copies_flipped_mean <- plot_timeseries_target() %>% print #format_logscale_y()
ggsave(plot_as(title_name, '-copy-flip-mean'), copies_flipped_mean, width = 7, height = 5)
ggplotly(copies_flipped_mean)

# OTHER TARGETS
copies_bb <- plot_timeseries_target('backbone') %>% print
ggsave(plot_as(title_name, '-copy-bb'), width = 7, height = 5)

copies_chr <- plot_timeseries_target('chromosome') %>% print
ggsave(plot_as(title_name, '-copy-chr'), width = 7, height = 5)
ggplotly(copies_chr)

# Processed plots ---- 

# logscale
copies_flipped_mean %>% format_logscale_y()

# plotting only uninduced
copies_flipped_unind <- filter(processed_data, `AHL (uM)` == 0) %>% 
  {plot_timeseries_target(.data = ., .connect = 'ind')} %>% print

format_logscale_y(copies_flipped_unind) # logscale
ggplotly(copies_flipped_unind, dynamicTicks = T) # interactive

ggsave(plot_as(title_name, '-copy-flip-unind'), width = 7, height = 5)

# backbone, only uninduced
copies_bb_unind <- filter(processed_data, `AHL (uM)` == 0) %>% 
  {plot_timeseries_target('backbone', .data = ., .connect = 'ind')} %>% print


# Normalized plot ----

plt_normalized <- plot_timeseries_target(.connect = 'ind', 
                                         .data = normalized_data, .yvar = normalized_Copies) %>% print

normal_zoom <- {plt_normalized + facet_wrap(vars(plasmid, organism), scales = 'free')} %>% print
ggsave(plot_as(title_name, '-flipped_normalized'), width = 8, height = 7)

# interactive
ggplotly(plt_normalized, dynamicTicks = T)
