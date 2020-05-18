# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./general_functions.R') # Source the general_functions file before running this


# User inputs -------------------------------------------------------------


# User inputs: choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject)

flnm <- 'Int_Assay1 MHT' 
flpath <- str_c('excel files/', flnm, '.xls') 
title_name <-'Arabinose induced flipping'
experiment_mode <- 'assay' # options ('small_scale' ; 'assay') ; future implementation: 'custom'. Explanation below
  # 'assay' =  Plots for Assays (facetted by Sample category = control vs experiment ; naming: 'Sample Name'_variable primer pair)
  # 'small_scale' = plots for troubleshooting expts : faceted by primer pair and sample name = template
plot_select_template <- '' # Options ('' or 'something') ; filters a particular template name to plot 

# small_scale mode features


# Assay mode features (choose if you want absolute quantification)
plot_assay_variable <- '[Arabinose] (mM)' # printed on the x axis of the graph
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); Determines which variable is chosen for plotting in different colours
plot_mode <-  'absolute_quantification'  # Options : ('absolute_quantification' or ''); absolute_quantification will calculate copy #'s based on intercept and slope from standard curve - manually entered below ; else, Cq values are plotted
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('Flipped', 'Unflipped', 'Backbone'),
  slope =  c(-3.36, -3.21, -3.55),
  intercept = c(42, 38, 42)
)
plot_normalized_backbone <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to backbone 
plot_mean_and_sd <- 'no' # Options: ('yes' or 'no'); plots mean and errorbars instead of each replicate as a point: Only in absolute_quantification mode
plot_exclude <- '^MHT' # Regex pattern: 'Controls2', '^MHT*' or 'none'; exclude categories for plotting; ex: Controls etc.: filters based on `Sample Name`: works only in assay mode


# Data loading and manipulation ------------------------------------------------------------


# Loading the data ----

# reading in file and polishing
fl <- readqpcr(flpath) # read excel file exported by Quantstudio

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
results_raw <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, `Ct Mean`, starts_with('Tm'),`Target Name`,`Target`) %>%  .[sample_order,] # select only the results used for plotting, calculations etc. and arrange them according to sample order

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
  
  # plot the CT mean along with replicates
  plt <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 2, show.legend = T) +
    geom_boxplot(aes(x = `Sample Name`, y = `Ct Mean`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + ylab(expression(C[q])) + facet_wrap(~Target, scales = 'free_x') 
  
  print(plt)
}


# Plots for Assays (facetted by Sample category; naming: 'Sample Name'_variable primer pair) ----
if (experiment_mode == 'assay')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  
  # isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
  results_relevant <- results_raw %>% separate(.,`Sample Name`,c('Sample Name','Primer pair'),' ') %>% separate(.,`Sample Name`,c('Sample Name','assay_variable'),'_')
  
  # Factorise the sample name in the order for plotting       re-arrange the results in same order as the above factors (columnwise order of the plate)
  results_relevant %<>% mutate_if(is.character,as_factor) %>% arrange(`Well Position`) 
  
  # select samples to plot (or to exclude write a similar command)
  results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x
  
  # plot the Tm of multiple peaks in melting curve ; Graph will now show
  plttm_assay <- plotalltms_assay(results_relevant) # plots tms of multiple peaks in melting curve
  
  # more data processing
  results_relevant %<>% filter(!str_detect(`Sample Name`, plot_exclude)) # exclude unwanted samples categories (sample_name) 

  results_relevant %<>% mutate(sample_type  = `Sample Name`, `Sample Name` = str_replace(`Sample Name`, 'GFP|MHT.*', 'Reporter')) # change GFP and MHT to reporter - for plot facetting
  results_relevant %<>% mutate(`Sample Name` = if_else(str_detect(assay_variable, 'Glu'), 'Controls', as.character(`Sample Name`)), `Sample Name` = fct_inorder(`Sample Name`), `Sample Name` = fct_relevel(`Sample Name`,'Controls') ) # change glucose data point into controls and plot controls first
  
  
  # Plot absolute quantification copy # : inferred from standard curve parameters (input in the start)
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    results_relevant_grouped <- results_relevant %>% group_by(Target) 
    results_abs <- results_relevant_grouped %>% do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target
    
    # changing text from assay variables into numbers (to be plotted on logscale) (change them back in inkscape)
    init_var <- c('Wa.*', 'rM.*', 'Gluc.*', 'fM.*','N'); fin_var <- c('.01','.03','.1','.3', '.01'); names(fin_var) = init_var; # N (NTC) merged with Water
    results_abs %<>% mutate( assay_variable = str_replace_all(assay_variable, fin_var), assay_variable = as.numeric(assay_variable)) 
    
    
    if(plot_mean_and_sd == 'yes') {
      y_variable = quo(mean)
      results_abs %<>% group_by(`Sample Name`, Target, assay_variable, sample_type) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd))
      } 
    else {y_variable = quo(`Copy #`)}
    
    # calculate hill function fit
    reporter_set <- results_abs %>% filter(str_detect(`Sample Name`,'Rep')) # filter only reporter values for hill function fitting
    hill_param <- reporter_set %>% rename(L = assay_variable, y = !!y_variable) %>%  hill_fitting_fn() # call hill fitting function on only reporter samples
    reporter_set  %<>% mutate(hill_fit = predict(hill_param)) # take the fit curve for plotting
    
    plt <- results_abs %>% ggplot(aes(x = `assay_variable`, y = !!y_variable, colour = sample_type)) + ylab('Copy #')    # plotting only mean
      
    if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = .1))} 
    
  } else plt <- results_relevant %>% ggplot(aes(x = `assay_variable`, y = CT, color = !!plot_colour_by))+ ylab(expression(C[q]))   
    
  # plot the CT mean along with replicates
  plt <- plt + geom_point(size = 2, show.legend = T) + geom_line(data = reporter_set, aes(x = assay_variable, y = hill_fit), linetype = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') # plot points and facetting
    
  plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale_x() # formatting plot, axes labels, title and logcale x plotting
  
  print(plt.formatted %>% format_logscale()) # print the formatted plot

  # normalizing copy #s to backbone ----  
  if(plot_mode == 'absolute_quantification' & plot_normalized_backbone == 'yes')
  { # computing ratio of copy #s of targets : flipped and unflipped to the backbone
    
    sel <- results_abs %>% select(`Sample Name`,assay_variable,`Primer pair`,Target,`Copy #`) # select relevant columns (other numeric columns will throw errors)
    
    sel_b <- sel %>% filter(Target == 'Backbone') # filter out each Target
    sel_f <- sel %>% filter(Target == 'Flipped'); sel_u <- sel %>% filter(Target == 'Unflipped');
    
    sel_f %<>% mutate("Normalized copy #" = sel_f$`Copy #`/sel_b$`Copy #`); # make ratios to the backbone 
    sel_u %<>% mutate("Normalized copy #" = sel_u$`Copy #`/sel_b$`Copy #`);
    
    results_ratio <- bind_rows(sel_f, sel_u) # bind results into 1 tibble (for easy plotting)
    
    # plotting the normalized copy #'s
    plt_norm <- results_ratio %>% ggplot(aes(x = `assay_variable`, y = `Normalized copy #`, color = Target)) +   # plotting
      scale_y_log10(  # logscale for y axis with tick marks
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
    
    plt_norm <- plt_norm + geom_point(size = 2) +
      theme_classic() + scale_color_brewer(palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
      ggtitle(title_name) + facet_wrap(~`Sample Name`, scales = 'free_x')
    
    print(plt_norm)
  }
}

# plate reader data ----
plateflnm <- '../Plate reader/plate reader data/S03_Ara memory GFP_26-9-18.xlsx' # file name for plate reader data

plate_data <- read_xlsx(plateflnm, sheet = 'in water', range = 'G68:J76')
# plate_data <- read_tsv(clipboard(), col_names = T) # read table from clipboard - really lazy
plate_data$`[Arabinose]`[1:2] <- c(10,0) # glucose will be 10 and 0 will be 0
plate_data %<>% mutate(category = c('Control', rep('Reporter',7)))

plate_subset <- plate_data %>% filter(`[Arabinose]` <= 1)
hill_plate <- plate_subset %>% rename(L = `[Arabinose]`, y = GFP ) %>% hill_fitting_fn() # call hill fitting function on only reporter samples
plate_subset  %<>% mutate(GFP.fit = predict(hill_plate)) # take the fit curve for plotting

plt.plate <- ggplot(plate_data, aes(`[Arabinose]`, GFP)) + geom_point(size = 2) + geom_line(data = plate_subset, aes(x = `[Arabinose]`, y = GFP.fit), linetype = 2) + ylab('GFP/OD (a.u.)') +
  facet_grid(~ category, scales = "free_x", space = "free_x")
plt.plate_formatted <- plt.plate %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() %>% format_logscale_x() # formatting plot, axes labels, title and logcale plotting

print(plt.plate_formatted) # print the formatted plot


# Custom plots (Transformed Assay data; Plots copy #; 'Sample Name'_variable 'primer pair') ----

# Save plots manually - copy paste this command to console
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)

# assay mode : plotting different colour and facet variables
# plt <- results_abs %>% ggplot(aes(x = `assay_variable`, y = `Copy #`, color = `Sample Name`)) + ylab('Copy #') +   # plotting
#   scale_y_log10(  # logscale for y axis with tick marks
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x) ))
# 
# plt + geom_point(size = 2, show.legend = T) +
#   theme_classic() + scale_color_brewer(palette="Set1") +
#   theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) +
#   ggtitle(title_name) + xlab(plot_assay_variable) + facet_wrap(~Target, scales = 'free_x')
