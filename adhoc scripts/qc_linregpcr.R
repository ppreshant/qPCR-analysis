# adhoc testing of rdmlpython's linregpcr output data for sanity check
# Author: Prashant Kalvapalle;  March 1 2022

# run prelims, fl and plate layout from analysis.R
  # Not repeating here to avoid redundancy
# since this code might be some day integrated with analysis.R with a user switch

# run test_sliwin.R 15-45 after loading excel file
# run linregpcr_post-process.R



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


