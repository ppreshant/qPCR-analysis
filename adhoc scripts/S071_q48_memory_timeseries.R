# S071_q48_memory_timeseries.R

# copied elements from S050_q37_41 ; subsumed: S067_q45_plots.R --> change inputs to do this


# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script


# User inputs -----
# override filename and title name from user inputs 
title_name <- 'q48_S071'

flnms <- c('q48a_S071_22-4-23', 
           'q48b_S071_23-4-23',
           'q48c_S071_d2_26-4-23',
           'q48d_S071_d8_26-4-23')


# select the column name to use for copies -- extrapolated as 2^ (40-Cq) vs absolute quantification
column_for_copies <- quo(Copies.per.ul.template)  # or Copies_proportional for uncalibrated data

# specify targets - as expressions : required for calculate_memory_ratios()
flipped <- expr(flipped)
backbone <- expr(backbone)
chromosome <- expr(chromosome)

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

grouping_vars_for_ratio <- metadata_columns # to group after pivoting by target, to take medians etc.

source('scripts_general_fns/22-memory_wrappers_ratio.R')
ratio_data <- calculate_memory_ratios(processed_data)


# # take ratio to backbone
# ratio_data <- select(processed_data, -CT) %>% # remove the non unique columns
#   pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
#   
#   mutate(plasmid_copy_number = backbone/chromosome,
#          flipped_fraction = flipped/backbone)



# timeseries plotting function ---- 

#' plots timeseries (by day) of signal for memory constructs
plot_timeseries_target <- function(filter_target = 'flipped', .connect = 'mean', point_style = 'straight',
                                   .data = processed_data, .yvar = Copies_proportional)
{
  data_to_plot <- filter(.data, organism != 'control', plasmid != 'No memory', # remove empty data
         if_any(matches('Target_name'), ~ .x == filter_target)) %>% 
    
    ungroup() # ungroup if it was grouped..
  
  
  mean_data <- reframe(data_to_plot,
                       .by = any_of(c(metadata_columns, 'Target_name')),
                       across(where(is.numeric), ~ mean(.x, na.rm = T)))
  
  # points jitter
  if(point_style == 'jitter')
  {pos_random <- position_jitter(width = 0.3, height = 0, seed = 1)
  
  } else pos_random <- position_identity()
  
  
  data_to_plot %>% # filter specific target
    
    ggplot(aes(day, {{.yvar}}, colour = plasmid, shape = `AHL (uM)`)) + 
    
    geom_point(size = .7, position = pos_random) +
    # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
    scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
    scale_alpha_discrete(guide = 'none', range = c(0.2, 0.5)) + # control line transparency
    scale_x_continuous(breaks = c(-1, 0, 2, 5, 8)) + # simplify x axis ticks
    
    # line connect data / means
    {if(.connect == 'mean')
    {geom_line(aes(alpha = `AHL (uM)`),
               data = mean_data, # join means
               position = pos_random)
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


# Ratio plot ----

plt_ratio_mean <- plot_timeseries_target(.data = ratio_data, .yvar = flipped_fraction,
                                         filter_target = 'ratio') %>% print
ggsave(plot_as(title_name, '-fraction-facets'), width = 7, height = 5)

present_ratio_mean <- 
  {plt_ratio_mean + ggtitle(NULL, subtitle = NULL) + # remove title
      ylab('ON state fraction of plasmid') + guides(colour = guide_legend('Designs')) + 
      scale_colour_manual(values = SS_colourscheme)
  } %>% print

ggsave('qPCR analysis/q48_S071_-all_flip_fraction.pdf', width = 6, height = 8)


# selected samples ----

# plotting only E. coli and W4 for fig 5 - clearer story
ratio_sel <- filter(ratio_data, str_detect(organism, 'Ec|w4')) %>% 
  {plot_timeseries_target(.data = ., .yvar = flipped_fraction,
                          filter_target = 'ratio')} %>% print


# presentable plot - clean up

# colours for memory picked by SS : coolors.co/
SS_colourscheme <- c('#9E2A2B', '#73956F', '#003844', '#F3C677', '#003844') # 'red', light green, light yellow, blue, magenta

present_flip_fraction <- 
  {ratio_sel + ggtitle(NULL, subtitle = NULL) + # remove title
      ylab('ON state fraction of plasmid') + guides(colour = guide_legend('Designs')) + 
      scale_colour_manual(values = SS_colourscheme)
      } %>% print

ggsave(plot_as(title_name, '-flip_fraction'), width = 7, height = 4)
# ggsave('qPCR analysis/Archive/q41_S050_flip_fraction.png', width = 7, height = 4)
ggsave('qPCR analysis/q48_S071_flip_fraction.pdf', width = 7, height = 4)



# plotting copies of flipped, maybe for suppl.?
copies_flipped_sel <- filter(processed_data, str_detect(organism, 'Ec|w4')) %>% 
  {plot_timeseries_target(.data = ., .connect = 'ind')} %>% print

# ggsave(plot_as(title_name, '-copy-flip-selected'), width = 7, height = 3)



# Statistics -----

# gp_noahl <- metadata_columns[!str_detect(metadata_columns, 'AHL')]
                             
condensed_data <- ratio_data %>% 
  group_by(across(all_of(metadata_columns[!str_detect(metadata_columns, 'AHL')]))) %>% 
  filter(str_detect(day, '^1|8'), str_detect(plasmid, 'Frugal'), str_detect(organism, 'Ec|w4')) %>% # select only d1 (first day measured for sil, fli) and d8
  nest() %>% print
  # view

stat_data <- 
  mutate(condensed_data, 
         ttest = map(data, 
                     ~ t.test(flipped_fraction ~ `AHL (uM)`, alternative = 'greater', paired = T,
                              data = .x)),
         
         pval = map_dbl(ttest, ~ .x$p.value)
         
  ) %>% print



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



# Normalized plot ----

plt_normalized <- plot_timeseries_target(.connect = 'ind', 
                                         .data = normalized_data, .yvar = normalized_Copies) %>% print

normal_zoom <- {plt_normalized + facet_wrap(vars(plasmid, organism), scales = 'free')} %>% print
ggsave(plot_as(title_name, '-flipped_normalized'), width = 8, height = 7)

# interactive
ggplotly(plt_normalized, dynamicTicks = T)
