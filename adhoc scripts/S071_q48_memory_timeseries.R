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
         across(day, as.numeric)) %>%  # convert day to numeric
  
  # change order for plotting
  mutate(across(plasmid, ~ fct_relevel(.x, c('Fluorescent', 'Silent', 'Frugal'))))


# Ratios ----

# take ratio to backbone
ratio_data <- select(processed_data, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone)


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

  
# timeseries plotting function ---- 

plot_timeseries_target <- function(filter_target = 'flipped', .connect = 'mean',
                                   .data = processed_data, .yvar = Copies_proportional)
{
  data_to_plot <- filter(.data, organism != 'control', plasmid != 'No memory', # remove empty data
         if_any(matches('Target_name'), ~ .x == filter_target)) 
  
  
  mean_data <- reframe(data_to_plot,
                       .by = any_of(c(metadata_columns, 'Target_name')),
                       across(where(is.numeric), ~ mean(.x, na.rm = T)))
  
  data_to_plot %>% # filter specific target
    
    ggplot(aes(day, {{.yvar}}, colour = plasmid, shape = `AHL (uM)`)) + 
    
    geom_point(size = 2) +
    # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
    scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
    scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
    scale_x_continuous(breaks = c(-1, 0, 2, 5, 8)) + # simplify x axis ticks
    
    # line connect data / means
    {if(.connect == 'mean')
    {geom_line(aes(alpha = `AHL (uM)`),
               data = mean_data) # join means
               # data = filter(processed_mean_toplt, Target_name == filter_target)) # join means
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


# individual target plot ----

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


# Ratio plot ----

plt_ratio_mean <- plot_timeseries_target(.data = ratio_data, .yvar = flipped_fraction,
                                         filter_target = 'ratio') %>% print
ggsave(plot_as(title_name, '-fraction-facets'), width = 7, height = 5)


# More plots ---- 

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


# selected samples ----

copies_flipped_sel <- filter(processed_data, str_detect(organism, 'Ec|w4')) %>% 
  {plot_timeseries_target(.data = ., .connect = 'ind')} %>% print

# ggsave(plot_as(title_name, '-copy-flip-selected'), width = 7, height = 3)

ratio_sel <- filter(ratio_data, str_detect(organism, 'Ec|w4')) %>% 
  {plot_timeseries_target(.data = ., .yvar = flipped_fraction,
                          filter_target = 'ratio')} %>% print
  
ggsave(plot_as(title_name, '-ratio-selected'), width = 7, height = 4)



# Normalized plot ----

plt_normalized <- plot_timeseries_target(.connect = 'ind', 
                                         .data = normalized_data, .yvar = normalized_Copies) %>% print

normal_zoom <- {plt_normalized + facet_wrap(vars(plasmid, organism), scales = 'free')} %>% print
ggsave(plot_as(title_name, '-flipped_normalized'), width = 8, height = 7)

# interactive
ggplotly(plt_normalized, dynamicTicks = T)
