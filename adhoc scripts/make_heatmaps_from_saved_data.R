# make_heatmaps_from_saved_data.R

source('./0-general_functions_main.R') # Source the general_functions file before running this


# Source user inputs  ----
source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm)

# Plotting heatmap ----
source('scripts_general_fns/19-qPCR_data_to_heatmap.R')


# all targets ~ heatmap + barcharts
qPCR_multitarget_heatmap()

ggsave(plot_as(title_name, '-cross contam'), width = 8, height = 5)


# q040 e : join by replicates
plot_facetted_assay()


# single target
qPCR_data_to_heatmap()

ggsave(plot_as(title_name, '-cross contam-flipped'), width = 5, height = 5)

