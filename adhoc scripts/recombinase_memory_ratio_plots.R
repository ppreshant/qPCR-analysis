# run analysis.R first
# this makes the plots relevant for memory

# regular plot 
pltd <- plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = Sample_category, 
                            .colourvar_plot = assay_variable, flipped_plot = FALSE) + 
  geom_line(aes(group = assay_variable))

ggsave(plot_as(title_name), height = 4, width = 6)


# flipped - join replicates
plt_flipped_repl <- 
  {plot_facetted_assay(.data = ratio_data, 
                       .yvar_plot = flipped, .xvar_plot = Sample_category, 
                       .colourvar_plot = assay_variable,
                       flipped_plot = FALSE, .facetvar_plot = NULL) + 
      
      geom_line(aes(group = interaction(assay_variable, biological_replicates), alpha = 0.2 ))} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, 'flipped-repl'), width = 4, height = 4)


# pivot and make ratios ----
ratio_data <- select(forplotting_cq.dat, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone) %>% 
  
  group_by(assay_variable, Sample_category) %>% 
  mutate(median_flipped_fraction = median(flipped_fraction, na.rm = T),
         median_copy_number = median(plasmid_copy_number, na.rm = T))

# fraction flipped ----
plt_fraction <- 
  {plot_facetted_assay(.data = ratio_data, 
                       .yvar_plot = flipped_fraction, .xvar_plot = Sample_category, 
                       .colourvar_plot = assay_variable,
                       flipped_plot = FALSE, .facetvar_plot = NULL) + 
      geom_line(aes(y = median_flipped_fraction, group = assay_variable))} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, 'flipped-fraction'), width = 4, height = 4)


# fraction flipped / replicate ----

plt_fraction_repl <- 
  {plot_facetted_assay(.data = ratio_data, 
                       .yvar_plot = flipped_fraction, .xvar_plot = Sample_category, 
                       .colourvar_plot = assay_variable,
                       flipped_plot = FALSE, .facetvar_plot = NULL) + 
      
      geom_line(aes(group = interaction(assay_variable, biological_replicates), alpha = 0.2 ))} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, 'flipped-fraction-repl'), width = 4, height = 4)

# plasmid copy number ----
plt_copies <- 
  {plot_facetted_assay(.data = ratio_data, 
                       .yvar_plot = plasmid_copy_number, .xvar_plot = Sample_category, 
                       .colourvar_plot = assay_variable,
                       flipped_plot = FALSE, .facetvar_plot = NULL) + 
      geom_line(aes(y = median_copy_number, group = assay_variable))} %>% 
  
  format_logscale_y()

ggsave(plot_as(title_name, 'copy_number'), width = 4, height = 3)
