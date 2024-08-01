# dead ribozyme?: re-plotting absolute copies of q18 data (pAOS054)

# params ----

# read the proper file 
flnm <- 'q18_dead U64 and isolated_4-10-21' # mention the file name without the "-linreg" or "-processed" suffixes

# read processed data ----
absolute_dat <- get_processed_datasets(flnm)
# replaced the old processed filename with raw/Cq values with this..

# subset data ----
absolute_subset <- 
  
  # filter data
  filter(absolute_dat, 
         str_detect(assay_var.label, 'Ribo|empty|NTC')) %>% 
  
  # make an order for plotting with factors
  mutate(across(assay_var.horz_label, factor, levels = rev(c('Ribo 295', 'dead-Ribo 54', 'empty 103', 'NTC')) ))



# plotting ----

horz.copies_w.mean <- {plot_facetted_assay(.data = absolute_subset, 
                                           .xvar_plot = assay_var.horz_label, 
                                           .xaxis.label.custom = axislabel.assay_variable,
                                           .yvar_plot = Copies.per.ul.template, 
                                           points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)} %>% 
  
  format_logscale_y()

# save plot
ggsave(plot_as(flnm, '-dead-ribo'), horz.copies_w.mean, width = 6, height = 6)


# fold change show ----

# show mean values
U64_mean_values <-
filter(absolute_subset, 
       assay_var.horz_label %in% c('Ribo 295', 'dead-Ribo 54'),
       Target_name == 'U64'
       ) %>% 
  select(assay_var.horz_label, mean_Copies.per.ul.template) %>% 
  unique() %>% 
  print()

# calculate fold change : Ribo 295 / dead-Ribo 54
U64_mean_values %>% 
  pull(mean_Copies.per.ul.template) %>%
  {.[2] / .[1]}



# summary fig ----

# name translators ----
# INSPIRED FROM S10c_q25_repression_summary.R

# used only for presentation stuff (extra analysis)
target_translator <- c('16s' = '16S rRNA',
                       'gfpbarcode' = 'cat-RNA',
                       'U64' = 'barcoded rRNA')

# Processing data ------------------------------

# create a subset of data for presenting : ONLY Ribozyme and dead-Ribozyme
processed_data <- filter(absolute_subset,
                         str_detect(assay_var.label, 'Ribo'), # Select only Ribozyme and dead-Ribozyme
                         !str_detect(Target_name, '16s') # Remove 16s (not plotting)
                         ) %>%  
  
 
  # rename target for ease
  mutate(across(Target_name, ~ str_replace_all(.x, target_translator)))
  


label_data <- 
  select(processed_data, 
         Target_name, assay_var.identifier, mean_Copies.per.ul.template) %>%
  unique() %>%  # select unique values
  
  # add a column for labels to cat-RNA (3,4)
  mutate(label_categories = c('Dead', 'Functional', NA, NA) )


# plot ------------------------------

{ggplot(processed_data, 
        aes(y = Target_name, x = Copies.per.ul.template, 
            shape = assay_var.identifier)) + 
    
    geom_jitter(width = 0, height = 0.2, alpha = 0.5) + 
    
    # change shapes
    scale_shape_manual('Expression', values = c('dead-Ribo' = 1, 'Ribo' = 19), ) +
    # remove legend manually?
    
    # Indicate the mean values
    geom_point(aes(x = mean_Copies.per.ul.template), 
               shape = '|', size = 5) + 
    
    # label values
    geom_text(aes(x = mean_Copies.per.ul.template,
                  # label = scales::label_scientific(digits = 1)(mean_Copies.per.ul.template),
                  label = scientific_10(mean_Copies.per.ul.template), # for power notation
                  
                  vjust = 2),  # adjust position of label
              
              parse = TRUE,  # parse the label
              data = label_data) +
    
    # direct label instead of legend
    geom_text(aes(x = mean_Copies.per.ul.template,
                  label = label_categories,
                  vjust = -2),  # adjust position of label
              
              data = label_data) +
    
    # theme and remove the legend
    theme_classic() + 
    
    theme(legend.position = 'none') +
    
    # rearrange order intuitively (reverse the alphabetical order)
    scale_y_discrete(limits = rev(c('barcoded rRNA', 'cat-RNA'))) +
    
    # adjust labels
    labs(title = NULL, x = 'Copies per ul template', y = NULL)
  
} %>%  # remove legend
  
  format_logscale_x() # make logscale on quantity on x-axis


# Save plot ------------------------------

ggsave(plot_as('q18_dead-Ribo_summary'), width = 6, height = 3) # save plot as png  

ggsave(filename = str_c('qPCR analysis/', 'q18_dead-Ribo_summary', '.pdf'), 
       dpi = 300, width = 4, height = 2.5)


# Save data ------------------------------

write.csv(processed_data, str_c('excel files/paper_data/', title_name, '.csv')) # save data as csv

# save clean data for paper
clean_data <- 
  select(processed_data, Target_name, Sample_category, 
         Copies.per.ul.template, mean_Copies.per.ul.template) %>%
  
  # change names to match figure
  mutate(across(Sample_category, ~ str_replace_all(.x, 'Maximal', 'Constitutive')))


output_path <- '../../Writing/RAM paper outputs/Archive'

write.csv(clean_data, str_c(output_path, '/', title_name, '.csv')) # save data as csv


# Statistics ------------------------------

# t-test for repressed vs constitutive
t_test_results <-
  
  filter(processed_data, Target_name == 'cat-RNA') %>%
  
  t.test(Copies.per.ul.template ~ Sample_category,
         data = .,
         alternative = 'greater') %>% # one-tailed test
  print()

t_test_results$p.value %>% fancy_scientific() # format p-value