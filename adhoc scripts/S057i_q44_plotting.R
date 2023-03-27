# S057i_q44_plotting.R


# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- 'S057i_marine lysate qPCR_q44'

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 
  
  filter(str_detect(Target_name,'genomic|plasmid|x$')) %>%  # select samples for the dataset
  mutate(across(Target_name, ~ str_replace(., 'Target-x', 'genomic')))

# plotting ----

horz.cq <- plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                               .xaxis.label.custom = axislabel.assay_variable)

ggsave(plot_as(title_name), width = 5, height = 5)


# Analysis ----

inhibition_analysis <- filter(forplotting_cq.dat, !str_detect(assay_variable, 'gDNA|nolysis|sybr')) %>% 
  
  mutate(across(assay_variable, ~ str_replace(., '1x', '1e0')),
         dilution = as.numeric(assay_variable),
         
         original_copies = Copies_proportional * dilution)

plt_inh <- plot_facetted_assay(.data = inhibition_analysis,
                    .yvar_plot = original_copies, .xvar_plot = assay_var.horz_label, 
                    .xaxis.label.custom = axislabel.assay_variable) %>% 
  format_logscale_y()

ggsave(plot_as(title_name, '-inhibition'), width = 5, height = 5)
