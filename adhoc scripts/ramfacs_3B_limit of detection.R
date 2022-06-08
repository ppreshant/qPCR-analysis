# custom analysis qPCR data: for fig 3B/ramfacs

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnm <- 'q32_S046 E coli dils_6-6-22'  
title_name <-'3B Limit of detection of splicing-qPCR'


# Sample name modifiers 

sample_name_translator <- c('Base strain' = 'Inf', # changes the LHS into the RHS
                            'ntc' = 'Inf',
                            '51' = '0.5',
                            '1/|,' = '') 

# Input the data ----

.df <- get_processed_datasets(flnm)

# Processing ----

# remove samples that are not relevant
forplotting_cq.dat <- filter(.df, !str_detect(Sample_category, 'conjug')) %>% 
  mutate('Plasmid containing cell fraction' = 
           (str_replace_all(assay_variable, sample_name_translator) %>% 
              as.numeric %>% {0.5/.} # convert to numbers and to fraction (1:1 -> 0.5)
            ),
         .before = assay_var.horz_label)


# Plotting ----

plt_copies <- {plot_facetted_assay(.yvar_plot = Copies.per.ul.template, .xvar_plot = `Plasmid containing cell fraction`,
                                  .xaxis.label.custom = NULL) +
    theme(legend.position = 'top')} %>% 
  format_logscale_x() %>% format_logscale_y()


plt_cq <- plt_copies <- {plot_facetted_assay(.yvar_plot = 40 - CT, .xvar_plot = `Plasmid containing cell fraction`,
                                             .xaxis.label.custom = NULL) +
    theme(legend.position = 'top')} %>% 
  format_logscale_x()

# interactive


# Save plot

ggsave(plot_as(title_name), plt_copies, width = 4.5, height = 4)
