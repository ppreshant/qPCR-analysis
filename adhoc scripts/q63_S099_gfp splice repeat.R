# q63_S099_gfp splice repeat.R
# plots modified versions of the typical plot. 
# It's a simple code but saving it here for consistency/future ease of access

# comparison to 328 data ; work in progress

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name


# read processed data ----
forplotting_cq.dat <- get_processed_datasets(flnm)

# plot 40 - Cq ----
plt.cq <- plot_facetted_assay(.yvar_plot = 40-CT, 
                              .xaxis.label.custom = "Template name", 
                              flipped_plot = FALSE)

# save plot with proper sizing (custom to the number of panels in the data)
ggsave(plot_as(title_name), height = 3, width = 6)