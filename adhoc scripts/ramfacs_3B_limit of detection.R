# custom analysis qPCR data: for fig 3B/ramfacs

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnm <- 'q32_S046 E coli dils_6-6-22'  
title_name <-'3B Limit of detection of splicing-qPCR'


# Sample name modifiers 

sample_name_translator <- c('Base strain' = 'Inf', # changes the LHS into the RHS
                            'ntc' = 'Inf',
                            '51' = '0',
                            '1/|,' = '') # remove commas and convert the 1/x into x 

# Input the data ----

.df <- get_processed_datasets(flnm)

# Processing ----

# remove samples that are not relevant
forplotting_cq.dat <- filter(.df, !str_detect(Sample_category, 'conjug')) %>% 
  mutate('fraction of RAM cells' = 
           (str_replace_all(assay_variable, sample_name_translator) %>% 
              as.numeric %>% {1/(1+.)} # convert to numbers and to fraction (1:1 -> 0.5, 1:10 -> 1/11 ..)
            ),
         .before = assay_var.horz_label)


# Plotting ----

plt_copies <- {plot_facetted_assay(.yvar_plot = Copies.per.ul.template, .xvar_plot = `fraction of RAM cells`,
                                  .xaxis.label.custom = NULL) +
    theme(legend.position = 'top')} %>% 
  format_logscale_x() %>% format_logscale_y()


plt_cq <- {plot_facetted_assay(.yvar_plot = 40 - CT, .xvar_plot = `fraction of RAM cells`,
                                             .xaxis.label.custom = NULL) +
    theme(legend.position = 'top')} %>% 
  format_logscale_x()

# Single panel plot
plt_monopanel_copies_w_mean <- {plot_facetted_assay(.yvar_plot = Copies.per.ul.template, .xvar_plot = `fraction of RAM cells`,
                                             .facetvar_plot = NULL, .colourvar_plot = Target_name,
                                   .xaxis.label.custom = NULL) +
    
    # geom_point(aes(y = mean_Copies.per.ul.template), shape = '-', size = 5, show.legend = FALSE) + # show the means
    geom_smooth(aes(y = mean_Copies.per.ul.template), 
                data = filter(forplotting_cq.dat, `fraction of RAM cells` > 1e-5), # constrain data
                method = 'lm', 
                alpha = .01, show.legend = FALSE) +
    
    theme(legend.position = 'top')} %>% 
  format_logscale_x() %>% format_logscale_y()


# Save plot ----

ggsave(plot_as('q32_Cq'), plt_copies, width = 4.5, height = 4)
ggsave(plot_as(title_name), plt_copies, width = 4.5, height = 4)
ggsave(plot_as('q32_copies', '-w line'), plt_monopanel_copies_w_mean, width = 3.5, height = 4)


# Save data ----

# Save formatted data (goes with publication)

# arrange important columns first
select(forplotting_cq.dat,
       'fraction of RAM cells', Target_name, Copies.per.ul.template, mean_Copies.per.ul.template, 
       Sample_category,
       everything()) %>%
  
write_csv(str_c('excel files/paper_data/', 'ramfacs-', title_name, '.csv', sep = ''),
          na = '')
