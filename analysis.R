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

# small_scale mode features


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

# Plots for small scale assays: Meant for troublshooting data (facetted by primer names; naming: 'Sample Name' primer-pair)----

if (experiment_mode == 'small_scale')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  results_relevant %<>%  separate(.,`Sample Name`,c('Sample Name','Primer pair'),' ')
  
  # Factorise the sample name in the order for plotting
  results_relevant %<>% mutate_if(is.character,as_factor) 
   
  # select samples to plot (or to exclude write a similar command)
  results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x
  
  # plot the Tm ; Graph will now show
  plttm <- plotalltms(results_relevant) # plots tms of multiple peaks in melting curve
  
  # plot the CT along with replicates
  plt <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 2, show.legend = T) +
    geom_boxplot(aes(x = `Sample Name`, y = `CT`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + ylab(expression(C[q])) + facet_grid(~Target, scales = 'free_x', space = 'free_x') 
  
  print(plt)
}



# Plots for Assays ----
# (facetted by Sample category; naming: 'Sample Name'_variable primer pair)

if (experiment_mode == 'assay')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  
  # isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
  results_relevant %<>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% separate(`Sample Name`,c('Sample Name','assay_variable'),'_') %>% mutate(assay_variable = if_else(`Sample Name` == 'NTC', 'NTC', assay_variable))
  
  results_relevant %<>% separate(assay_variable, c('assay_variable', 'biological_replicates'))
  
  # Factorise the sample name in the order for plotting
  results_relevant %<>% mutate_if(is.character,as_factor) 
  
  # re-arrange the results in same order as the above factors (columnwise order of the plate)
  results_relevant %<>% arrange(`Well Position`) 
  
  # select samples to plot (or to exclude write a similar command)
  results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x
  
  results_relevant %<>% filter(!str_detect(`Sample Name`, plot_exclude_category)) # exclude unwanted samples categories (sample_name) 
  results_relevant %<>% filter(!str_detect(assay_variable, plot_exclude_assay_variable)) # excluding unwanted samples from assay_variable
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    results_relevant_grouped <- results_relevant %>% group_by(Target) 
    results_abs <- results_relevant_grouped %>% do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target
    
    if(plot_mean_and_sd == 'yes') {
      y_variable = quo(mean)
      results_abs %<>% group_by(`Sample Name`, Target, assay_variable) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd)) # find mean and SD of individual copy #s for each replicate
      } 
    else {y_variable = quo(`Copy #`)}
    
    plt <- results_abs %>% ggplot(aes(x = `assay_variable`, y = !!y_variable, color = !!plot_colour_by)) + ylab('Copy #')    # Specify the plotting variables 

    if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = errorbar_width))} # plot errorbars if mean and SD are desired
    
  } 
  
  else plt <- results_relevant %>% ggplot(aes(x = `assay_variable`, y = CT, color = !!plot_colour_by))+ ylab(expression(C[q])) # plot CT values if absolute quantification is not needed
    
  # Formatting plot
  plt <- plt + geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting
  plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
  
  print(plt.formatted)

  # normalizing copy #s to backbone ----  
  if(plot_mode == 'absolute_quantification' & plot_normalized_backbone == 'yes')
  { # computing ratio of copy #s of targets : flipped and unflipped to the backbone
    
    sel <- results_abs %>% select(`Sample Name`,assay_variable,Target,`Copy #`) # select relevant columns (other numeric columns will throw errors)
    
    sel_b <- sel %>% filter(Target == 'Backbone') # filter out each Target
    sel_f <- sel %>% filter(Target == 'Flipped'); sel_u <- sel %>% filter(Target == 'Unflipped');
    
    sel_f %<>% mutate("Normalized copy #" = sel_f$`Copy #`/sel_b$`Copy #`); # make ratios to the backbone 
    sel_u %<>% mutate("Normalized copy #" = sel_u$`Copy #`/sel_b$`Copy #`);
    
    results_ratio <- bind_rows(sel_f, sel_u) # bind results into 1 tibble (for easy plotting)
    
    # plotting the normalized copy #'s
    plt_norm <- results_ratio %>% ggplot(aes(x = `assay_variable`, y = `Normalized copy #`, color = Target)) +   # plotting
    geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting
    
    plt_norm.formatted <- plt_norm %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
    
    print(plt_norm)
  }
}

# ggsave('qPCR analysis/WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 5, height = 4)
