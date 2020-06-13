# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./general_functions.R') # Source the general_functions file before running this

# User inputs ----
# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) 
# Sample naming guide for plate template: in google sheet 'enter cell address of experiment name' (check 'plate_template_raw' for link)
  # BCoV-608_A.2; ignore-facet_(x axis variable).biological replicate (BCoV is ignored, only for pipetting reference, actual target is taken from the qPCR results)

flnm <- 'WW9_BCoV_troubleshooting'  # set the filename

title_name <- flnm
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('BRSV_N', 'BCoV_M', 'N1_CoV2', 'N2_CoV2'),
  slope =  c(-3.62, -3.49, -3, -3.12),
  intercept = c(39, 39, 39, 40) # values for various targets
)

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

# Data input ----

# Preperation steps
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
plate_template_raw <- read_sheet('https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118', sheet = 'Plate import setup', range = 'G1:S9')

# Read in qPCR data and labels from plate template
fl <- readqpcr(flpath) # read excel file exported by Quantstudio
plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') # convert plate template (sample names) into a single vector, columnwise

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order

# Load desired qPCR result sheet and columns
bring_results <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`) %>% rename(Target = `Target Name`) %>%  .[sample_order,] %>%  # select only the results used for plotting, calculations etc. and arrange them according to sample order
 select(-`Sample Name`) %>% right_join(plate_template, by = 'Well Position') # Incorporate samples names from the google sheet by matching well position
 
# Remove unneccesary data
rm(fl, plate_template_raw)  # remove old data for sparsity

# Data polishing ----


# Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)

# isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
polished_results <- bring_results %>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% separate(`Sample Name`,c('Sample Name','Tube ID'),'_') %>% 
  mutate(`Tube ID` = if_else(`Sample Name` == 'NTC', '0', `Tube ID`)) %>% 
  separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F) %>%  # Separate out biological replicates 
  mutate_if(is.character,as_factor) %>% arrange(`Well Position`) # Factorise the sample name and rearrange in column order of appearance on the plate (for plotting)

# select samples to plot (or to exclude write a similar command)
results_relevant <- polished_results %>% filter(str_detect(`Sample Name`, paste('^', plot_select_facet, sep = ''))) %>%  # Include only desired facets : str_detect will find for regular expression; ^x => starting with x
  filter(!str_detect(`Sample Name`, plot_exclude_facet)) %>%  # exclude unwanted facets (sample_name) 
  filter(!str_detect(assay_variable, plot_exclude_assay_variable)) # excluding unwanted x axis variables from assay_variable

# Computation ----


# Computing copy number from standard curve linear fit information
results_abs <- results_relevant_grouped %>% group_by(Target) %>% do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target

# Finding mean and standard deviation within replicates (both technical and biological)
if(plot_mean_and_sd == 'yes') {
  y_variable = quo(mean)
  concise_results_abs <- results_abs %>%  group_by(`Sample Name`, Target, assay_variable) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd)) # find mean and SD of individual copy #s for each replicate
  results_abs$`Copy #` %<>% replace_na(0) # make unamplified values 0 for plotting
  data_to_plot <- concise_results_abs
  
  } else {y_variable = quo(`Copy #`); data_to_plot <- results_abs}

plt <- data_to_plot %>% ggplot(aes(x = `assay_variable`, y = !!y_variable, color = !!plot_colour_by)) + ylab('Copies/ul RNA extract')    # Specify the plotting variables 

if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = errorbar_width)) + # plot errorbars if mean and SD are desired
  geom_jitter(data = results_abs, aes(x = `assay_variable`, y = `Copy #`, alpha = map_chr(`Copy #`, ~. == 0), size = map_chr(`Copy #`, ~. == 0)), colour = 'black', width = .2, show.legend = F) + # plot raw data
  scale_alpha_manual(values = c(.3, 1)) + scale_size_manual(values = c(1, 2)) } # manual scale for emphasizing unamplified samples

# Formatting plot
plt <- plt + geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting
plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting

print(plt.formatted)

# Data output ----
# this is usually commented out to prevent overwriting existing data; turn on only when needed for a single run

# write_sheet(results_abs,'https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=0', sheet = title_name) # save results to a google sheet
# ggsave('qPCR analysis/WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 5, height = 4)