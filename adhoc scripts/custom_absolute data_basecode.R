# Read in the processed data for specific plotting (papers etc.)
# Author: Prashant Kalvapalle;  Sep 4 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnm <- 'q10_RAM_So-pp-ps_19-7-21'  # raw filename, without the "-processed" keywork
title_name <-'q10_RAM So Pp Ps'

# options

# Labelling translators ----

# Subtitle labeller (for y axis variables)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template')

organism_translation <- c('So' = 'Shewanella', # regex based translation to make full name
                          'Pp' = 'putida',
                          'Ps' = 'stutzeri') 

assay_var_translation <- c('328' = 'Ribo', # regex based translation to change x-axis labels
                             '314' = '(-)')

# Input the data ----

# reading in file and polishing
.df <- read_csv(str_c('excel files/processed_data/',flnm, '-processed.csv')) # read excel file exported by Quantstudio


# Mess with the data ----

# any custom processing goes here
forplot_reduced.data <- filter(.df, str_detect(Sample_category, 'test')) %>% # filter necessary data  
  
  
  separate(assay_var.label,  # separate organism names into a new column
           into = c('organism', 'assay_var.label'),
           sep = ' ',
           fill = 'left') %>%  # for NTC, fill the left one with NA
  
  # translate to fullish name for organism
  mutate(across(organism, ~ str_replace_all(.x, organism_translation)) ) %>%  
  replace_na(replace = list('organism' = '')) %>%  # replace the NA of organism with empty string
  
  # translate to meaningful x axis labels
  mutate(across(assay_var.label, ~ str_replace_all(.x, assay_var_translation)) ) 

# Plotting ----

plt.copies_w.mean <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                                          .yvar_plot = Copies.per.ul.template,
                                          .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                          points_plt.style = 'jitter') + 
    
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + # plot mean
    ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


# Save plot

ggsave(str_c('qPCR analysis/Archive/', title_name, '.png'),
       plt.copies_w.mean,
       width = 8,
       height = 4)
