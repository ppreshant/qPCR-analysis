# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./general_functions.R') # Source the general_functions file before running this

# User inputs ----
# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) 

flnm <- 'WW1_Baylor-bovine_test'  # set the filename
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
plate_template_raw <- read_sheet('https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118', sheet = 'Plate import setup', range = 'G1:S9')

title_name <-'Baylor RNA extracts_pilot run'
experiment_mode <- 'assay' # options ('small_scale' ; 'assay') ; future implementation: 'custom'. Explanation below
  # 'assay' =  Plots for Assays (facetted by Sample category = control vs experiment ; naming: 'Sample Name'_variable primer pair)
  # 'small_scale' = plots for troubleshooting expts : faceted by primer pair and sample name = template
plot_select_template <- '' # Options ('' or 'something') ; filters a particular template name to plot 
errorbar_width = 0.1; # width of errorbars - emperically change

# Assay mode features (choose if you want absolute quantification)
plot_assay_variable <- 'Sample' # printed on the x axis of the graph
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); Determines which variable is chosen for plotting in different colours
plot_mode <-  'absolute_quantification'  # Options : ('absolute_quantification' or ''); absolute_quantification will calculate copy #'s based on intercept and slope from standard curve - manually entered below ; else, Cq values are plotted
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('BRSV_N', 'BCoV_M', 'N1_CoV2', 'N2_CoV2'),
  slope =  c(-3.62, -3.49, -3, -3.12),
  intercept = c(39, 39, 39, 40) # values for various targets
)
plot_normalized_backbone <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to backbone 
plot_mean_and_sd <- 'yes' # Options: ('yes' or 'no'); plots mean and errorbars instead of each replicate as a point: Only in absolute_quantification mode
plot_exclude_category <- '^none' # Regex pattern: 'Controls2', '^MHT*', '^none; exclude categories for plotting; ex: Controls etc.: filters based on `Sample Name`: works only in assay mode
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; exclude assay_variables for plotting; ex: no template control etc.: filters based on assay_variable: works only in assay mode

# Input the data ----

# reading in file and polishing
fl <- readqpcr(flpath) # read excel file exported by Quantstudio

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
results_relevant <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`) %>% rename(Target = `Target Name`) %>%  .[sample_order,] # select only the results used for plotting, calculations etc. and arrange them according to sample order

plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') # convert plate template (sample names) into a single vector, columnwise
results_relevant %<>% mutate(`Sample Name` = plate_template$`Sample Name`) # Incorporate samples names from the google doc 
results_relevant$Target %<>% str_replace('BSRV', 'BRSV') # correcting mis-spelled name of BRSV target

rm(fl, plate_template_raw)  # remove old data for sparsity

# Plots for Assays ----
# (facetted by Sample category; naming: 'Sample Name'_variable primer pair)

# Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)

# isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
results_relevant %<>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% separate(`Sample Name`,c('Sample Name','Tube ID'),'_') %>% mutate(`Tube ID` = if_else(`Sample Name` == 'NTC', '0', `Tube ID`))

results_relevant %<>% separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F)

# Factorise the sample name in the order for plotting
results_relevant %<>% mutate_if(is.character,as_factor) 

# re-arrange the results in same order as the above factors (columnwise order of the plate)
results_relevant %<>% arrange(`Well Position`) 

# select samples to plot (or to exclude write a similar command)
results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x

results_relevant %<>% filter(!str_detect(`Sample Name`, plot_exclude_category)) # exclude unwanted samples categories (sample_name) 
results_relevant %<>% filter(!str_detect(assay_variable, plot_exclude_assay_variable)) # excluding unwanted samples from assay_variable


# Computing copy number from standard curve linear fit information
results_relevant_grouped <- results_relevant %>% group_by(Target) 
results_abs <- results_relevant_grouped %>% do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target

if(plot_mean_and_sd == 'yes') {
  y_variable = quo(mean)
  concise_results_abs <- results_abs %>%  group_by(`Sample Name`, Target, assay_variable) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd)) # find mean and SD of individual copy #s for each replicate
  data_to_plot <- concise_results_abs
  
  } else {y_variable = quo(`Copy #`); data_to_plot <- results_abs}

plt <- data_to_plot %>% ggplot(aes(x = `assay_variable`, y = !!y_variable, color = !!plot_colour_by)) + ylab('Copy #')    # Specify the plotting variables 

if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = errorbar_width)) + # plot errorbars if mean and SD are desired
  geom_jitter(data = results_abs, aes(x = `assay_variable`, y = `Copy #`), colour = 'black', size = 1, alpha = .2, width = .2)} # plot raw data


# Formatting plot
plt <- plt + geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting
plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting

print(plt.formatted)

# ggsave('qPCR analysis/WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 5, height = 4)