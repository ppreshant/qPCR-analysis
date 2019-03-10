# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# Source the general_functions file before running this
# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

# User inputs: choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject) ----

flnm <- 'excel files/S4_1 intrinsic flip.xls'  
title_name <-'Intrinsic flipping with time assay'
experiment_mode <- 'small_scale' # options ('small_scale' ; 'assay') ; future implementation: 'custom'. Explanation below
  # 'assay' =  Plots for Assays (facetted by Sample category = control vs experiment ; naming: 'Sample Name'_variable primer pair)
  # 'small_scale' = plots for troubleshooting expts : faceted by primer pair and sample name = template

# small_scale mode features
plot_select_template <- 'rMHT' # Options ('' or 'something') ; filters out a particular template name to plot 

# Assay mode features (choose if you want absolute quantification)
plot_mode <-  'absolute_quantification'  # Options : ('absolute_quantification' or ''); absolute_quantification will calculate copy #'s based on intercept and slope from standard curve - manually entered below ; else, Cq values are plotted
f_slope <- -3.36; f_intercept <- 42 
r_slope <- -3.23; r_intercept <- 38
plot_exclude <- '' # quo('Controls2') or ''; exclude categories for plotting; ex: Controls etc.: filters based on `Sample Name`: works only in assay mode

# plotting functions for Melting temperature ----

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(fl)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant%>% select(`Sample Name`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample Name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(fl)
{ 
  plttm <- results_relevant%>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}


# Input the data ----

# reading in file and polishing
fl <- readqpcr(flnm) # read excel file exported by Quantstudio

sample_order = order_columnwise(fl) # this orders the samples columnwise in the PCR plate or strip (data is shown row-wise) => This command will enable plotting column wise order
results_relevant <- fl$Results %>% select(`Sample Name`, CT, `Ct Mean`, starts_with('Tm')) # select only the results used for plotting, calculations etc.

# Plots for small scale assays: Meant for troublshooting data (facetted by primer names; naming: 'Sample Name' primer-pair)----

if (experiment_mode == 'small_scale')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  results_relevant <- separate(results_relevant,`Sample Name`,c('Sample Name','Primer pair'),' ')
  
  # Factorise the sample name in the order for plotting
  results_relevant$`Sample Name` %<>% factor(levels = unique(.[sample_order]))
  results_relevant$`Primer pair` %<>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  
  # select samples to plot or to exclude
  if(plot_select_template != '')  {results_relevant %<>% filter(str_detect(`Sample Name`, paste('^', plot_select_template, sep = '')))} # str_detect will find for regular expression; ^x => starting with x
  
  # plot the Tm ; Graph will now show
  plttm <- plotalltms(fl) # plots tms of multiple peaks in melting curve
  
  # plot the CT mean along with replicates
  plt <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 1, show.legend = T) +
    geom_boxplot(aes(x = `Sample Name`, y = `Ct Mean`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + ylab(expression(C[q])) + facet_wrap(~`Primer pair`, scales = 'free_x') 
  
  print(plt)
}


# Plots for Assays (facetted by Sample category; naming: 'Sample Name'_variable primer pair) ----
if (experiment_mode == 'assay')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  
  # isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
  results_relevant %<>% separate(.,`Sample Name`,c('Sample Name','Primer pair'),' ') %>% separate(.,`Sample Name`,c('Sample Name','assay_variable'),'_')
  
  # Factorise the sample name in the order for plotting
  results_relevant$`Sample Name` %<>% factor(levels = unique(.[sample_order]))
  results_relevant$`Primer pair` %<>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  results_relevant$`assay_variable` %<>% factor(levels = unique(.[sample_order])) # assay_variable
  
  # plot the Tm of multiple peaks in melting curve ; Graph will now show
  
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample Name`, `assay_variable`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`, -`assay_variable`)
  
  # plot the Tm ; Graph will now show
  plttm <- tmfl %>% ggplot(.) + aes(x = `assay_variable`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting')) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
  if(plot_exclude != '')  {results_relevant %<>% filter(`Sample Name` != (!!plot_exclude))}
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    flf <- results_relevant %>% filter(`Target Name` == 'fMHT') %>% mutate(`Copy #` = 10^( (CT - f_intercept)/f_slope ) )  # not vectorized
    flr <- results_relevant %>% filter(`Target Name` == 'rMHT') %>% mutate(`Copy #` = 10^( (CT - r_intercept)/r_slope ) )  # not vectorized
    fljoin <- bind_rows(flf, flr)
    
    plt <- fljoin %>% ggplot(aes(x = `assay_variable`, y = `Copy #`, color = `Sample Name`)) +   # plotting
      scale_y_log10(  # logscale for y axis with tick marks
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
  } else plt <- results_relevant %>% ggplot(aes(x = `assay_variable`, y = CT, color = `Sample Name`))+ ylab(expression(C[q]))   
    
  # plot the CT mean along with replicates
  plt <- plt + geom_point(size = 1, show.legend = T) +
    # geom_line(aes(x = `assay_variable`, y = `Ct Mean`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
  print(plt)
  
}


# Custom plots (Transformed Assay data; Plots copy #; 'Sample Name'_variable 'primer pair') ----

# Save plots manually - copy paste this command to console
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)