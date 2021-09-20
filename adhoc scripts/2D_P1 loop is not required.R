# Author: Prashant Kalvapalle;  Sep 20 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this


# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnm <- 'q12_328-330-RNAsimple_19-8-21'  # raw filename, without the "-processed" keywork
title_name <-'P1 loop not required_q12'

# Input the data ----

# reading in file and polishing
.df <- read_csv(str_c('excel files/processed_data/',flnm, '-processed.csv')) # read excel file exported by Quantstudio


# Mess with the data ----

# any custom processing goes here
forplot_reduced.data <- filter(.df, str_detect(Sample_category, 'GMO kit'),
                               !str_detect(assay_variable, 's328|NTC')) %>% # filter necessary data  
  
  # separate(assay_var.label, into = c(assay_var.label, NA), sep = '\n') %>%  # display x axis for publication
  mutate(across(assay_var.identifier, ~ str_replace(., 'empty', '(-)')))

# Plotting ----

plt.copies_w.mean <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                                          .yvar_plot = Copies.per.ul.template,
                                          .xvar_plot = assay_var.identifier,
                                          
                                          .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                          points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) # plot mean
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


# Save plot ----

ggsave(str_c('qPCR analysis/Archive/', title_name, '.png'),
       plt.copies_w.mean,
       width = 5,
       height = 4)

# Save data ----

# write data for repository in google drive
write_csv(forplot_reduced.data,
          str_c('excel files/paper_data/', title_name, '.csv', sep = ''),
          na = '')
