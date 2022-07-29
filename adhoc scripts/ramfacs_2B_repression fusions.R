# Load multiple files for comparative analysis

# fig 2B will have only one plasmid set (51) ; 
# Using another figure with all fusions for presentations

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q028b_77 79-repressed only_18-5-22',
           'q28_77 79-fusion2_13-4-22',
           'q25_S037_RAM repression_14-2-22')

title_name <- 'Fusion : U64 RAM and fluorescence'

# name translators ----

assay_var_translation <- c('76' = 'Ribozyme',    # 'current name' = 'new name' format
                           # '48' = 'mScarlet',  # regex based translation to change x-axis labels
                           
                           '51' = 'Ribozyme-mScarlet',
                           '77' = 'mScarlet-Ribozyme v2',
                           '79' = 'mScarlet-Ribozyme',
                           
                           '67 \\+ ' = '', # remove extra information
                           ' \\+ 67' = '', # the special symbol '+' needs to be double escaped
                           'NTC' = 'ntc')

sample_category_cleaner <- c('test' = 'Repressed', # clean up redundant names
                             'repressed' = 'Repressed',
                             'positive' = 'Maximal')


# Input the data ----

# reading in files and merge : along with the run ID
.df <- get_processed_datasets(flnms)


# Processing ----

forplotting_cq.dat <- filter(.df, !str_detect(Sample_category, 'lysate')) # remove samples that are not relevant

# create a subset of data for presenting the most relevant data : 51, 77, 79 ; 48, 76 with repression and without
presentation.dat <- filter(forplotting_cq.dat, 
                           !str_detect(Sample_category, 'conjug'), # remove samples that are unimportant for presentation
                           !str_detect(assay_variable, 'Sort|sort|328|^67$|48'),
                           str_detect(Target_name, 'U64')) %>% 
  
  mutate('assay_var.label' = str_replace_all(assay_variable, assay_var_translation)) %>%  # assign labels for readability
  mutate(across(Sample_category, ~ str_replace_all(.x, sample_category_cleaner))) %>%  # clean up redundant names
  
  # arrange for plotting
  mutate(across(Sample_category, ~ fct_relevel(.x, 'Repressed', 'negative', 'Maximal'))) %>%  # set the order of the colours r, b, g
  mutate(across(assay_var.label, ~ fct_reorder(.x, Copies.per.ul.template))) %>%  # descending order of copies
  mutate(across(assay_var.label, ~ fct_relevel(.x, 'ntc'))) # Take ntc to the end -- first factor is being plotted closes to origin..
  

# Plot - presentation ----

plt.copy_ppt <-
  
  {plot_facetted_assay(.data = presentation.dat, .yvar_plot = Copies.per.ul.template, .xvar_plot = assay_var.label, .facetvar_plot = NULL) + 
      
      coord_flip() + 
      # facet_grid(rows = vars(Target_name), scales = 'free_y', space = 'free_y') +
      theme(legend.position = 'top', legend.title= element_blank()) + # stylistic editing of legend 
      xlab('')} %>% 
  format_logscale_y()

# save plot
ggsave(plot_as('all fusions, U64-red'), height = 3, width = 5)
  
# Plotting - all data ----

# plot 40 - Cq
plt.cq_straight <- plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_variable)

# Horizontal orientation : for label readability
horz.cq <- {plt.cq_straight + 
    coord_flip() + 
    facet_grid(rows = vars(Target_name), scales = 'free_y', space = 'free_y') +
    theme(legend.position = 'top') + 
    xlab('') + ylab('40 - Cq')} %>% 
  print()

# save plot ----

ggsave(plot_as('all fusions and conjug, U64-red'), plot = horz.cq, width = 7, height = 5)


# Interactive plot ----

inter_horz.cq <- horz.cq + aes(text = run_ID, label = biological_replicates)
ggplotly(inter_horz.cq)

