# S091_combined_analysis.R

# Prelims ----
source('./0-general_functions_main.R') # Source the general_functions file before running this


# User inputs
flnms <- c('S091-marine-withsort-5-5-24',
           'S091b-marine-higher conc') %>% 
  
  str_c('-ddPCR') # attach ddPCR to the filenames
  

title_name <- 'S091-marine-ddPCR'


# load data ----

combined_data <- get_processed_datasets(flnms) %>% 
  relocate(run_ID) %>% # bring run_ID to the front
  
  # get 'b' in S0xxb ; else below will be 'a'
  mutate(batch = if_else(str_detect(run_ID, 'b'), 'b', 'a'),
         .after = run_ID)


# Processing ----

# (pcn): Ratio: plasmid/chromosome Swetha's marine ddPCR data
ratio_data <- combined_data %>% 
  
  # remove flagged wells
  filter(!flag_wells == 'Manual',  # remove manually flagged wells
         
         # remove wells below LOD for chromosome
         !(flag_wells == 'Below LOD' &  Target_name == 'Chromosome') 
  ) %>% # remove negatives [junk ratios]
  
  # calculate copy number : ratio of plasmid to chromosome
  select(Sample_category, assay_variable, Well, Target_name, Concentration, batch) %>% 
  pivot_wider(names_from = Target_name, values_from = Concentration) %>% 
  mutate(copy_number = Plasmid / Chromosome)


# adhoc ----
filter(combined_data,
       flag_wells == 'Manual') %>%
       # (flag_wells == 'Below LOD' &  Target_name == 'Chromosome')) %>% 
  relocate(PositiveDroplets, LOD_max_pos_droplets) %>%
  view


# Plotting ----

## copies by target ----

# plot copies per target
copies_all_targets <- 
  plot_facetted_assay(.data = 
                        filter(combined_data, !flag_wells == 'Manual'), # removed flagged wells
                      
                      .yvar_plot = CopiesPer20uLWell, .xvar_plot = assay_variable, 
                      .label.var = Well,
                      .facetvar_plot = Target_name) + 
  
  # add LOD line
  geom_hline(data = LOD, aes(yintercept = LOD_max_copieswell), linetype = 'dashed', color = 'red') +
  geom_hline(data = LOD, aes(yintercept = LOD_mean_copieswell), linetype = 'dashed', color = 'gray')

# geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

# save the plot
ggsave(plot_as(title_name, '-ddPCR_raw'), plot = copies_all_targets,
       width = 4.2, height = 4)

# interactive plot
ggplotly(copies_all_targets)


## plasmid copy number (ratio) ----

copy_num <- 
  plot_facetted_assay(.data = ratio_data, 
                      .yvar_plot = copy_number, .xvar_plot = assay_variable, 
                      .label.var = Well,
                      .facetvar_plot = NULL) 

# copy_num + 
#   geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

# save the plot
ggsave(plot_as(title_name, '-ddPCR_plasmid-copies'), plot = copy_num, 
       width = 4.2, height = 3)

# interactive plot
ggplotly(copy_num)





# adhoc ----
if(0)
{
  # look at negatives ----
  filter(compiled_w_metadata, Sample_category == 'Negative', 
         !flag_wells == 'Manual') %>% 
    
    relocate(Well, Target_name, PositiveDroplets, flag_wells) %>%
    
    arrange(assay_variable, Target_name, desc(PositiveDroplets)) %>% 
    
    view
}
