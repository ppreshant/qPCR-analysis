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


# linregpcr ----

# run linregPCR on spyder