# Read in the qPCR file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  June 13 2020


# User inputs ----
# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) 
# Sample naming guide for plate template: in google sheet 'enter cell address of experiment name' (check 'plate_template_raw' for link)
# BCoV-608_A.2; ignore-facet_(x axis variable).biological replicate (BCoV is ignored, only for pipetting reference, actual target is taken from the qPCR results)
# If code fails, first thing: check the number of lines to skip before the data begins and tally with the code (including the headings)

flnm <- 'WW23_622_BRSV_TR2'  # set the filename

title_name <- flnm
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('BRSV_N', 'BCoV_M', 'N1_CoV2', 'N2_CoV2', 'N1_multiplex',  'N2_multiplex'),
  slope =  c(-3.61, -3.47, -2.98, -3.12, -2.95, -2.94),
  intercept = c(38, 39, 39, 40, 32, 34) # values for various targets
)
template_volume <- 4 # ul template volume in qPCR reaction

# Other parameters ----

# Additional parameters (changed rarely)
plot_assay_variable <- 'Sample' # printed on the x axis of the graph
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); Determines which variable is chosen for plotting in different colours
errorbar_width = 0.1; # width of errorbars - emperically change
plot_normalized_to_recovery <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to recovery % 

# Obsolete variables (will be phased out soon)
# can be implemented easily with generalized plotting function
plot_mode <-  'absolute_quantification'  # Options : ('absolute_quantification' or ''); absolute_quantification will calculate copy #'s based on intercept and slope from standard curve - manually entered below ; else, Cq values are plotted
plot_mean_and_sd <- 'yes' # Options: ('yes' or 'no'); plots mean and errorbars instead of each replicate as a point: Only in absolute_quantification mode
experiment_mode <- 'assay' # options ('small_scale' ; 'assay') ; future implementation: 'custom'. Explanation below : can be plotting tube labels or WWTP names
# 'assay' =  Plots for Assays (facetted by Sample category = control vs experiment ; naming: 'Sample Name'_variable primer pair)
# 'small_scale' = plots for troubleshooting expts : faceted by primer pair and sample name = template

# inclusion exclusion variables
plot_select_facet <- '' # Options ('' or 'something') ; filters a particular template name to plot 
plot_exclude_facet <- '^none' # Regex pattern: 'Controls2', '^MHT*', '^none; exclude categories for plotting; ex: Controls etc.: filters based on `Sample Name`: works only in assay mode
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; exclude assay_variables for plotting; ex: no template control etc.: filters based on assay_variable: works only in assay mode
