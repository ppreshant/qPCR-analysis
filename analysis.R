# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./general_functions.R') # Source the general_functions file before running this

# User inputs: choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) ----

flnm <- 'excel files/S06d_no rox.xls'  
title_name <-'qPCR on cell lysates'
experiment_mode <- 'assay' # options ('small_scale' ; 'assay') ; future implementation: 'custom'. Explanation below
  # 'assay' =  Plots for Assays (facetted by Sample category = control vs experiment ; naming: 'Sample Name'_variable primer pair)
  # 'small_scale' = plots for troubleshooting expts : faceted by primer pair and sample name = template
plot_select_template <- '' # Options ('' or 'something') ; filters a particular template name to plot 

# small_scale mode features


# Assay mode features (choose if you want absolute quantification)
plot_assay_variable <- 'Sample' # printed on the x axis of the graph
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); Determines which variable is chosen for plotting in different colours
plot_mode <-  'absolute_quantification'  # Options : ('absolute_quantification' or ''); absolute_quantification will calculate copy #'s based on intercept and slope from standard curve - manually entered below ; else, Cq values are plotted
std_par <- tibble(                       # Input the slope and intercept from standard curve of various primer pairs/targets here - Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('Flipped', 'Unflipped', 'Backbone'),
  slope =  c(-3.36, -3.23, -3.55),
  intercept = c(42, 38, 42) # values for primer pairs: Flipped:q4-5. Unflipped:q9-10, Backbone:q12-13
)
plot_normalized_backbone <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to backbone 
plot_mean_and_sd <- 'yes' # Options: ('yes' or 'no'); plots mean and errorbars instead of each replicate as a point: Only in absolute_quantification mode
plot_exclude_category <- '^MHT*' # Regex pattern: 'Controls2', '^MHT*', '^none; exclude categories for plotting; ex: Controls etc.: filters based on `Sample Name`: works only in assay mode
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; exclude assay_variables for plotting; ex: no template control etc.: filters based on assay_variable: works only in assay mode

# plotting functions for Melting temperature ----

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample Name`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample Name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}


# Input the data ----

# reading in file and polishing
fl <- readqpcr(flnm) # read excel file exported by Quantstudio

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
results_relevant <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, `Ct Mean`, starts_with('Tm'),`Target Name`,`Target`) %>%  .[sample_order,] # select only the results used for plotting, calculations etc. and arrange them according to sample order

# Plots for small scale assays: Meant for troublshooting data (facetted by primer names; naming: 'Sample Name' primer-pair)----

if (experiment_mode == 'small_scale')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  results_relevant %<>%  separate(.,`Sample Name`,c('Sample Name','Primer pair'),' ')
  
  # Factorise the sample name in the order for plotting
  results_relevant %<>% mutate_if(is.character,as_factor) 
  # results_relevant$`Well Position` %<>% factor(levels = unique(.[sample_order]))
  # results_relevant$`Sample Name` %<>% factor(levels = unique(.[sample_order]))
  # results_relevant$`Primer pair` %<>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  # results_relevant$Target %<>% factor(levels = unique(.[sample_order]))
  
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
  results_relevant %<>% separate(.,`Sample Name`,c('Sample Name','Primer pair'),' ') %>% separate(.,`Sample Name`,c('Sample Name','assay_variable'),'_')
  
  # Factorise the sample name in the order for plotting
  results_relevant %<>% mutate_if(is.character,as_factor) 
  # results_relevant$`Well Position` %<>% factor(levels = unique(.[sample_order]))
  # results_relevant$`Sample Name` %<>% factor(levels = unique(.[sample_order]))
  # results_relevant$`Primer pair` %<>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  # results_relevant$`assay_variable` %<>% factor(levels = unique(.[sample_order])) # assay_variable
  # results_relevant$Target %<>% factor(levels = unique(.[sample_order]))
  # results_relevant %<>% mutate_if(is.character, as_factor(levels = arrange(sample_order))) or factor(levels = unique(.[sample_order])) # fancy way of vectorizing - doesn't work
  
  # re-arrange the results in same order as the above factors (columnwise order of the plate)
  results_relevant %<>% arrange(`Well Position`) 
  
  # select samples to plot (or to exclude write a similar command)
  results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x
  
  # plot the Tm of multiple peaks in melting curve ; Graph will now show
  
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample Name`, `assay_variable`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`, -`assay_variable`)
  
  # plot the Tm ; Graph will now show
  plttm <- tmfl %>% ggplot(.) + aes(x = `assay_variable`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting')) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
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

    if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = .25))} # plot errorbars if mean and SD are desired
    
  } 
  
  else plt <- results_relevant %>% ggplot(aes(x = `assay_variable`, y = CT, color = !!plot_colour_by))+ ylab(expression(C[q])) # plot CT values if absolute quantification is not needed
    
  # plot the CT mean and formatting plots
  plt <- plt + geom_point(size = 2) + facet_wrap(~`Sample Name`, scales = 'free_x') # plot points and facetting
  plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
  
  print(plt.formatted)

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
    geom_point(size = 2) + facet_wrap(~`Sample Name`, scales = 'free_x') # plot points and facetting
    
    plt_norm.formatted <- plt_norm %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
    
    print(plt_norm)
  }
}

print(plttm)

# ggsave('qPCR analysis/S017.png')
# write.xlsx(results_abs, 'excel files/Test.xls', sheetName = 'analysis', append = TRUE, borders = 'surrounding') # saving data table of inferred copy #s to an excel sheet

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
