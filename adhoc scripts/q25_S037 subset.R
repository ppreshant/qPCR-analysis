# adhoc code for q25 data : plotting a subset as a horizontal graph (similar to fluorescence)

# Run analysis.R on the desired file till line 269

# Cq data ----
forplotting_cq.dat %>% View()

repression.dat <- 
  forplotting_cq.dat %>% 
  filter(str_detect(assay_var.label, '51|76|67|328'),
         Target_name != '16s') %>%  # grab desired subset of samples
  mutate(across(assay_var.label, ~ str_replace(., 'Ribo\n', ''))) # Remove the helper from the variable 


repression.cq <- plot_facetted_assay(.data = repression.dat,  
                                     .yvar_plot = 40-CT)

rcq_horizontal <- 
  repression.cq + 
  coord_flip() + 
  facet_grid(rows = vars(Target_name)) +
  theme(legend.position = 'top')

ggsave(str_c('qPCR analysis/Archive/', flnm, '-concise.png'), 
       rcq_horizontal,
       width = 5,
       height = 3)

# additional horizontal plot for all data
plt.cq + 
  coord_flip() + 
  facet_grid(rows = vars(Target_name)) +
  theme(legend.position = 'top')


# linregpcr data ----

repression.linreg <- 
  linreg.selected %>% 
  filter(str_detect(assay_variable, '51|76|67|328'),
         Target_name != '16s') # grab desired subset of samples


present.linreg <- 
  repression.linreg %>% 
  select(all_of(metadata_unique_columns), mean_N0) %>% # select only the mean data
  
  arrange(Target_name) %>% 
  mutate(across(mean_N0, ~ formatC(.x, format = 'e', digits = 2))) %>% 
  unique()

write.csv(present.linreg, 'excel files/processed_data/q25_repressed_linreg.csv', na = '')

# linregpcr plots ----

# run linregPCR on spyder
# run qc_linregpcr till line 66 (linreg.selected)

repression_linreg.cq <- plot_facetted_assay(.data = repression.linreg,
                                            .xvar_plot = assay_variable,
                                            .yvar_plot = N0)

rlrcq_horizontal <- 
  {repression_linreg.cq + 
  geom_point(aes(y = mean_N0), shape = '|', size = 5, show.legend = FALSE) + 
  coord_flip() + 
  facet_grid(rows = vars(Target_name)) +
  theme(legend.position = 'top')} %>% 
  
  format_logscale_y() %>% 
  print()

ggsave(plot_as('q25_repression_linreg'), width = 5, height = 3)
