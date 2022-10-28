# Combine > 1 processed data and plot in different colours for comparison
# Run the pre-requisites from analysis.R initially

# User inputs ----

flnms <- c('q38_S055_8-orgs expression_22-10-22', # raw file name, without the "-processed" keyword
           'q038b_S055_10xdil_8-orgs-expression_26-10-22')

# Load processed data ----

proc_data <- map_dfr(flnms,
                     get_processed_datasets)

plt.comparison <- plot_facetted_assay(.data = proc_data, .yvar_plot = 40-CT, .colourvar_plot = run_ID)

ggsave(plot_as(base_title_name), width = 12, height = 4)


# additional analysis ----

assayvar_ordering <- 
  
  filter(proc_data, 
         run_ID == 'q038b',
         Target_name == 'U64') %>% # choose the latest file
  
  arrange(desc(Copies.per.ul.template)) %>% 
  mutate(across(assay_var.label, fct_inorder)) %>% 
  pull(assay_var.label) %>% 
  levels()

ordered_proc.data <- mutate(proc_data, 
                            across(assay_var.label, ~ fct_relevel(.x, assayvar_ordering)))

plt.comparison <- plot_facetted_assay(.data = ordered_proc.data, .yvar_plot = 40-CT, .colourvar_plot = run_ID)

ggsave(plot_as(base_title_name, '-compared'), width = 12, height = 4)
