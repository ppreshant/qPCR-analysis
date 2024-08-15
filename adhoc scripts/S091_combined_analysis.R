# S091_combined_analysis.R

# Prelims ----
source('./0-general_functions_main.R') # Source the general_functions file before running this


# User inputs
flnms <- c('S091-marine-withsort-5-5-24',
           'S091b-marine-higher conc') %>% 
  
  str_c('-ddPCR') # attach ddPCR to the filenames
  

title_name <- 'S091-marine-ddPCR'

## ordering -----
assay_var_order <- c('J23100', 
                     'high', 'med', 'low',
                     'high-growth', 'med-growth', 'low-growth',
                     'WT', 
                     'blank', 'NTC') # order the assay variables


# load data ----

combined_data <- get_processed_datasets(flnms) %>% 
  relocate(run_ID) %>% # bring run_ID to the front
  
  # fix NAs in flag_wells
  replace_na(list(flag_wells = '')) %>% 
  
  # get 'b' in S0xxb ; else below will be 'a'
  mutate(batch = if_else(str_detect(run_ID, 'b'), 'b', 'a'),
         .after = run_ID) %>% 
  
  ## custom processing ordering data ----
  mutate(across(assay_variable, ~ fct_relevel(.x, assay_var_order))) # order the assay variables


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


# TODO: calculate the mean in the same dataset or different?

# remove flagged wells (manual)
plotting_data <- combined_data %>% 
  filter(!flag_wells == 'Manual') # remove manually flagged wells


## get LOD ----

LOD <- select(combined_data, 
              Target_name, run_ID,
              LOD_mean_copieswell, LOD_max_copieswell, LOD_max_pos_droplets) %>% 
  unique() %>%
  
  # select the maximum ones between different runs
  reframe(across(matches('LOD'), max), .by = Target_name)


## adhoc/checking ----
if(0)
{
  
  # view stuff below LOD / == LOD_max_pos_droplets
  filter(combined_data,
         # flag_wells == 'Manual') %>%
    (flag_wells == 'Below LOD' &  Target_name == 'Chromosome'), # below LOD for chromosome
     Sample_category != 'Negative' # except negatives, which are obvious
    
        ) %>%
    
    relocate(PositiveDroplets, LOD_max_pos_droplets) %>%
    view
  
  # filter view
  filter(combined_data, flag_wells == 'Manual') %>% view
  
  # (obsolete) filter with NAs 
  # Source: https://stackoverflow.com/a/74756923/9049673
  filter_na <- function(tbl, expr){
    tbl %>% filter({{expr}} %>% replace_na(T))
  }
 
  
}


# save data ----

write_csv(ratio_data, 
          str_c('excel files/processed_data/', title_name, '-ratio-combined.csv'), # save the ratio data
          na = '') # remove NAs for clarity

# Plotting/paper ----

## plotting wrappers ----

#' plots grey points for replicates ; mean and 1 stdev range in black
#' mimics the figure2-b,c style of https://pubs.acs.org/doi/10.1021/ac501118v
#' @param plt_object ggplot object with data and aesthetics
#' @param ... Can pass additional arguments to the ggplot objects of 
#' `geom_point()` and `stat_summary()`
#' @return ggplot object
#' @example: ggplot(...) %>% geom_mean_w_replicates(shape = '+', aes(y = a2))
show_mean_w_replicates <- function(plt_object, ...)
{
  plt_object +
    
    # plot points in grey ; shift them to the left for clarity
    geom_point(colour = 'grey', position = position_nudge(x = -0.2), ...) + 
    
    # show means and stdev
    # source: HMisc package/ https://stackoverflow.com/a/41848876/9049673
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "pointrange", ...)
  
}


# Plot the individual copies of plasmid, chromosome; 2nd plot with plasmid copy number
#' @param .data_subset data subset
#' @param x_axis_label x axis label
#' @return ggplot object
#' @example: ggplot(...) %>% show_mean_w_replicates()
# TODO: need to show legend of the shapes (custom make?) in the raw copies/top plot
# TODO: provide conditionally outputs of intermediate plots?

plot_raw_copies_and_PCN <- function(.data_subset, x_axis_label)
{
  
  copy_number_plt <-
    
    ggplot(.data_subset, 
           aes(x = assay_variable, y = copy_number)) %>% 
    
    # plot points in grey (left shifted), mean and stdev in black
    show_mean_w_replicates() + 
    
    # facet by organism
    facet_grid(cols = vars(Sample_category), scales = 'free_x', space = 'free_x') +
    
    # label axes
    xlab(x_axis_label) +
    
    ylab(NULL) +
    labs(subtitle = 'Plasmid copy number')
  
  
  # plot individual targets
  
  # convert data into long format for plasmid and chromosome
  .data_long <- 
    select(.data_subset, -copy_number) %>% # remove copy number
    
    pivot_longer(cols = c(Plasmid, Chromosome), names_to = 'Target', values_to = 'copies')
  
  
  # plotting
  raw_copies_plt <- 
    
    ggplot(.data_long, 
           aes(x = assay_variable, y = copies, shape = Target)) %>% 
    
    # plot points in grey (left shifted), mean and stdev in black
    show_mean_w_replicates() + 
    
    # # Same for chromosome copies
    # show_mean_w_replicates(shape = '+', mapping = aes(y = Chromosome)) + 
    
    # shape by target
    scale_shape_manual(values = c(19, 3)) + 
    
    # facet by organism
    facet_grid(cols = vars(Sample_category), scales = 'free_x', space = 'free_x') +
    
    # thematic changes
    theme(legend.position = 'top') + # legend on top ;; c(0.1, .9) for top left
    
    # label axes
    xlab(x_axis_label) +
    ylab(NULL) +
    labs(subtitle = 'Copies of DNA/uL')
  
  
  # merge both plots
  library(patchwork)
  
  raw_copies_plt + copy_number_plt + 
    
    # layout with single column ; collect x axis labels
    plot_layout(ncol = 1, axis_titles = "collect")
  
}



## sorted data ----
# TODO : remove Rd (outlier) to zoom in

sorted_data <- filter(ratio_data, str_detect(assay_variable, 'low$|med$|high$'))

# make a plot of copy number by assay variable ; facet by sample_category
sorted_plot <- 
  plot_raw_copies_and_PCN(sorted_data, 'Bins of fluorescence, post sorting') %>% 
  print()

# save the plot
ggsave(plot_as('S091', '-sorted'),
       width = 7, height = 5)

## WT vs plasmid -----

lysate_data <- filter(ratio_data, str_detect(assay_variable, 'WT|J23100'))

# make a plot of copy number by assay variable ; facet by sample_category
cultures_plot <- 
  plot_raw_copies_and_PCN(lysate_data, 'Plasmid') %>% 
  print()

# save the plot
ggsave(plot_as('S091', '-lysate'),
       width = 7, height = 5)


## Sorted, post growth ----

sorted_growth <- 
  filter(ratio_data, str_detect(assay_variable, 'growth$')) %>% 
  
  # remove growth from assay_variable
  mutate(across(assay_variable, ~ str_remove(.x, '-growth')))
  

# make a plot of copy number by assay variable ; facet by sample_category

sort_growth_plot <- 
  plot_raw_copies_and_PCN(sorted_growth, 'Bins of fluorescence, post sorting and growth') %>% 
  print()


# save the plot
ggsave(plot_as('S091', '-sorted-after growth'),
       width = 7, height = 5)


# old code: plot copies per target faceted by organism

# copy_number_lysate <- 
#   plot_facetted_assay(.data = lysate_data, 
#                       .yvar_plot = copy_number, .xvar_plot = assay_variable,
#                       .colourvar_plot = NULL,
#                       .label.var = Well, 
#                       flipped_plot = FALSE, facet_scale_constraint = 'free',
#                       .facetvar_plot = Sample_category)

   
# Old plotting ---------------------------------

## copies by target ----

# plot copies per target
copies_all_targets <- 
  {plot_facetted_assay(.data = 
                         filter(combined_data, !flag_wells == 'Manual'), # removed flagged wells
                      
                      .yvar_plot = CopiesPer20uLWell, .xvar_plot = assay_variable, 
                      .label.var = Well,
                      .facetvar_plot = Target_name) + 
  
  # add LOD line
  geom_hline(data = LOD, aes(yintercept = LOD_max_copieswell), linetype = 'dashed', color = 'red') +
  geom_hline(data = LOD, aes(yintercept = LOD_mean_copieswell), linetype = 'dashed', color = 'gray')} %>% 
  
  print

# geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

# zoom in 
copies_all_targets + ylim (c(0, 100))

# save the plot
ggsave(plot_as(title_name, '-ddPCR_zoom'),
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
ggsave(plot_as(title_name, '-ddPCR_plasmid-copies'),
       width = 4.2, height = 3)

# interactive plot
ggplotly(copy_num)




# adhoc ----
if(0)
{
  ## look at negatives ----
  filter(compiled_w_metadata, Sample_category == 'Negative', 
         !flag_wells == 'Manual') %>% 
    
    relocate(Well, Target_name, PositiveDroplets, flag_wells) %>%
    
    arrange(assay_variable, Target_name, desc(PositiveDroplets)) %>% 
    
    view
}
