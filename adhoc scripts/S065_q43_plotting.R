# S065_custom_plotting.R
# works for q43 and q44 too


# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- 'q44_S068 frugal lysis protocol'

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 
  filter(str_detect(Target_name, 'flipped|backbone|chromosome')) %>% # select only the triplex samples
  
  mutate_cond(str_detect(Sample_category, 'Water'),  # clean up names for 'Water' samples
              across(assay_variable, ~ str_replace(., 'spillH4', '11')), # change the name to a number
              
              biological_replicates = as.numeric(assay_variable), 
              across(contains('assay_var'), ~ 'water'))

# Plotting ----

horz.cq <- {plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                               .xaxis.label.custom = axislabel.assay_variable) + 
    guides(colour = guide_legend(title = '', nrow = 2))} %>%  # remove legend title
  
  print()


horz.cq + ggrepel::geom_text_repel(aes(label = biological_replicates), show.legend = F,
                                    data = filter(forplotting_cq.dat, assay_variable == 'water'))


# interactive  
plotly::ggplotly(horz.cq)

# Save plot ----

ggsave(plot_as(title_name), width = 3.5, height = 5)


# Plotting heatmap ----
source('scripts_general_fns/19-qPCR_data_to_heatmap.R')

# single target
qPCR_data_to_heatmap(.target = 'backbone')

ggsave(plot_as(title_name, '-cross contam-backbone'), width = 5, height = 5)


# all targets ~ heatmap + barcharts
qPCR_multitarget_heatmap()

ggsave(plot_as(title_name, '-cross contam'), width = 5, height = 5)
