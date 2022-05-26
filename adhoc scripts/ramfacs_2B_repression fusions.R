# Load multiple files for comparative analysis

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q028b_77 79-repressed only_18-5-22',
           'q28_77 79-fusion2_13-4-22',
           'q25_S037_RAM repression_14-2-22')

title_name <- 'Fusion : U64 RAM and fluorescence'


# Input the data ----

# reading in files and merge
.df <- map_dfr(flnms, 
               ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv')) %>%  # read excel file exported by Quantstudio
                 mutate(run_ID = str_extract(.x, 'q[:alnum:]*')) # add the run_ID from the file name
)


# Processing ----

forplotting_cq.dat <- filter(.df, !str_detect(Sample_category, 'lysate'))
  
# Plotting ----

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

ggsave(plot_as('all fusions U64-red'), plot = horz.cq, width = 7, height = 5)


# Interactive plot ----

inter_horz.cq <- horz.cq + aes(text = run_ID, label = biological_replicates)
ggplotly(inter_horz.cq)
