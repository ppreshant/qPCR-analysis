# adhoc testing of rdmlpython's linregpcr output data for sanity check
# Author: Prashant Kalvapalle;  March 1 2022

# run prelims, fl and plate layout from analysis.R
  # Not repeating here to avoid redundancy
# since this code might be some day integrated with analysis.R with a user switch

# run test_sliwin.R 15-45 after loading excel file


flnm <- 'q25_S037_RAM repression_14-2-22-linreg'


# Data input ----

# Read the rdmlpython output file (tab separated tsv file)
linreg.results <- read_tsv(str_c('excel files/', flnm, '.csv')) %>%
 
  select(-sample) %>% # remove irrelevant columns
  
  left_join(plate_template, by = c('well' = 'Well Position')) %>%  
  # can rename to Well Position if needed later
  
  
  # clean up targets : Giving primacy to the data file, removing the plate_template targets
  mutate(Target_name = target, .keep = 'unused') %>% 
  
  # remove un-named samples
  filter(!is.na(Sample_category)) # remove all samples that don't have a Sample_category - implies plate layout was empty
# This is intended to remove samples whose labels have been removed to prevent analysis


# select useful data ----

# define all the metadata columns coming from the plate_layout
metadata_unique_columns <- c('Sample_category', 'assay_variable', 'Target_name')

# change the data with the windows as line equation y - y2 = m(x - x2)
linreg.selected <- linreg.results %>% 
    
  # select only relevant columns
  select(any_of(metadata_unique_columns), # metadata unique
         well, biological_replicates, # addition metadata for replicates
         `n in log phase`, `last log cycle`, `n included`,  # window range
         
         # linregpcr actual outputs
         `mean PCR eff`, `N0 (mean eff)`, `Cq (mean eff)`,
         
         # window anchors
         `baseline`, `log lin cycle`, `log lin fluorescence`) %>% 
  
  # rename columns for ease of reference while plotting and for equations
  rename(n = 'n in log phase',
         x0 = 'log lin cycle', # mid point of log linear phase?
         y0 = 'log lin fluorescence', # Already baseline subtracted value
         xend = 'last log cycle',
         n_sml = 'n included',
         eff = 'mean PCR eff',
         N0 = 'N0 (mean eff)',
         Cq = 'Cq (mean eff)') %>%
  
  # make a col for mean of N0 -- for plotting as a line/boxplot
  group_by(across(all_of(metadata_unique_columns))) %>% 
  mutate(mean_N0 = mean(N0, na.rm = TRUE))



# Overlay window on curve ----
# will add lines to overview_raw_plot/log to show the window of linearity detected with the E value


# Make data to plot the window of linearity
linreg_window.data <- 
  
  linreg.selected %>% 
  
  # Make data for the window
  mutate(window = 
           pmap(list(n, xend, eff, x0, y0),
                ~ tibble(x_cycle = (..2 - ..1 + 1):xend, # Make a series of points window start to end
                         # y_fl_log = ..5 + log10(..3) * (x_cycle - ..4)
                         y_fl_bs = (..5) * (..3)^(x_cycle - ..4)
                         )
                )
         ) %>% 
  
  unnest(cols = window)
    

# Join linreg_results to raw data
qpcr_polished.data <- qpcr_amplification.data %>% 
  
  # join data with linreg results
  left_join(linreg.selected) %>% 
  
  # subtract background
  mutate(Rn_bg = Rn - baseline)


# Plotting ----

# plot only the linear window
plt_linear_window <- 
  {ggplot(linreg_window.data,
         aes(x_cycle, 
             y_fl_bs,
             colour = Target_name,
             group = interaction(Target_name, biological_replicates), # ensures each replicate is a different line
             label = well)) +
  geom_line() + 
  # make nested facets : category with children of assay variables
  ggh4x::facet_nested_wrap(. ~ Sample_category + assay_variable)} %>% 
  
  format_logscale_y() %>% 
  print()

ggsave(plot_as('q25_window overlaid'), width = 5, height = 6)


# Overlay window on top of raw data

overview_bs_plt <- {ggplot(qpcr_polished.data,
                            aes(x = Cycle,
                                y = Rn_bg,
                                colour = Target_name,
                                group = interaction(Target_name, biological_replicates), # ensures each replicate is a different line
                                label = `Well Position`)) +
    geom_line() +

    # make nested facets : category with children of assay variables
    ggh4x::facet_nested_wrap(. ~ Sample_category + assay_variable) +


    # add lines for window of linearity -- DOES not match
    geom_line(data = linreg_window.data,
              mapping = aes(x_cycle,
                  y_fl_bs,
                  # colour = Target_name,
                  group = interaction(Target_name, biological_replicates), # ensures each replicate is a different line
                  label = well),
              colour = 'black'
              )
    
  } %>%

  format_logscale_y() %>%

  print()


# could get this to work
# plot_facetted_assay()

# adhoc plot
# N0_linreg_plt <- 
#   {linreg.results %>% 
#   ggplot(aes(assay_variable, `N0 (mean eff)`, colour = Target_name)) + 
#   geom_point() +
#   
#   facet_wrap(~ Sample_category, scales = 'free_x')} %>% 
#   
#   format_logscale_y() %>% 
#   print()

# Plot N0 along with mean on logscale
N0_linreg_plt <- 
  {plot_facetted_assay(.data = linreg.selected, 
                      .yvar_plot = N0,
                      .colourvar_plot = Target_name,
                      .facetvar_plot = Sample_category,
                      .xvar_plot = assay_variable,
                      points_plt.style = 'jitter') + 
      geom_boxplot(aes(y = mean_N0), show.legend = FALSE)
      
      } %>% 
  
  format_logscale_y() %>% 
  print()

# saving plot
ggsave(plot_as('S025_linregpcr_N0 log'), width = 6, height = 4)

# plot 40 - Cq to compare
title_name <- 'q25_linregpcr cq'
cq_linregp <- plot_facetted_assay(.data = linreg_window.data, 
                    .yvar_plot = 40 - Cq,
                    .colourvar_plot = Target_name,
                    .facetvar_plot = Sample_category,
                    .xvar_plot = assay_variable)
# 
ggsave(plot_as('q25_linregpcr cq-2'), width = 6, height = 4)

