# S067_q45_plots.R

# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script


# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 

  clean_up_water_wells() %>% 
  
  # split columns
  separate(assay_variable, into = c('plasmid', 'organism'), sep = '/', remove = FALSE) %>% 
  separate(Sample_category, into = c('day', 'induction'), sep = '/', remove = FALSE) %>% 
  
  # change types
  mutate(across(induction, 
                ~ if_else(str_detect(., 'I'), 10, 0) %>% # Induced = 10, uninduced
                  replace_na(0)))  # NA is 0

# TODO : take ratio of flipped / backbone ; plot ratio in a timecourse ; 


# Plotting ----


# timeseries plot from plate reader


# Small multiples-timeseries plot
processed.data %<>% mutate(across(Samples, ~ fct_relevel(., rev(sample_order)))) # order in desc of max

title_name <- 'S067b wastewater memory II'

filter(processed.data, category == 'Memory') %>% 
  
  ggplot(aes(Day, `RFP/OD`, colour = category, shape = `AHL (uM)`)) + 
  
  geom_point(size = 2) +
  # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
  scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
  
  # line like a dumbell plot
  geom_line(aes(group = interaction(`AHL (uM)`, category)), colour = 'black', alpha = 0.3) +  
  
  facet_wrap(vars(Samples), scale = 'free_y') +
  theme(legend.position = 'top') +
  
  # labels
  ggtitle('Memory in wastewater', subtitle = title_name)

ggsave(plot_as(title_name, '-timeseries-memory'), width = 5, height = 5)


# overview plot
ggplot(S67b2_order,
       aes(`RFP/OD`, Samples, colour = Day, shape = `AHL (uM)`)) + 
  
  geom_point(size = 2) +
  # scale_colour_brewer(palette = 'Dark2', direction = -1) + # change the values - orange for uninduced/0
  scale_shape_manual(values = c(1, 16)) + # shape : open and closed circles
  
  # line like a dumbell plot
  geom_line(aes(group = Samples), colour = 'black', alpha = 0.3) + 
  
  # facet_grid(vars(category), scale = 'free_y', space = 'free_y') +
  theme(legend.position = 'top') +
  
  # labels
  ggtitle('Memory in wastewater', subtitle = title_name)

ggsave(plot_as(title_name), width = 4, height = 6)



# Doesn't work
horz.cq <- {plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable) + 
    guides(colour = guide_legend(title = '', nrow = 2))} %>%  # remove legend title
  
  print()


horz.cq + ggrepel::geom_text_repel(aes(label = biological_replicates), show.legend = F,
                                   data = filter(forplotting_cq.dat, assay_variable == 'water'))
ggsave(plot_as(title_name, '-overview'))

# interactive  
plotly::ggplotly(horz.cq)



