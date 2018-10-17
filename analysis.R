# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# Source the general_functions file before running this
# The last part of the script contains important codes from the console that were used for specific situations : These will be commented


# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject)----
flnm <- 'excel files/Int_Assay1 MHT.xls'  
title_name <-'Integrase induction'
experiment_mode <- 'assay' # 'small_scale' ; 'assay' ; 'absolute_quantification'  ; 'custom'

# Assay mode features for absolute quantification
plot_mode <- 'absolute_quantification'  # CT
f_slope <- -3.36; f_intercept <- 42 
r_slope <- -3.23; r_intercept <- 38

# plotting functions ----

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(fl)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- fl$Results %>% select(`Sample Name`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample Name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(fl)
{ 
  plttm <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_wrap(~`Primer pair`, scales = 'free_x')
}


# Input the data ----

# reading in file and transformations
fl <- readqpcr(flnm) # read excel file exported by Quantstudio

sample_order = order_columnwise(fl) # this order is columnwise (data is shown row-wise) => useful for plotting column wise order

# Plots for small scale assays (facetted by primer names;typical troubleshooting; naming: 'Sample Name' primer-pair)----

if (experiment_mode == 'small_scale')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  fl$Results <- separate(fl$Results,`Sample Name`,c('Sample Name','Primer pair'),' ')
  
  # Factorise the sample name in the order for plotting
  fl$Results$`Sample Name` <- fl$Results$`Sample Name` %>% factor(levels = unique(.[sample_order]))
  fl$Results$`Primer pair` <- fl$Results$`Primer pair` %>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  
  # plot the Tm ; Graph will now show
  plttm2 <- plotalltms(fl) # plots tms of multiple peaks in melting curve
  
  # plot the CT mean along with replicates
  plt <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 1, show.legend = T) +
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
  
  # isolate the primer pair and [Arabinose] into 2 columns
  fl$Results <- separate(fl$Results,`Sample Name`,c('Sample Name','Primer pair'),' ')
  fl$Results <- separate(fl$Results,`Sample Name`,c('Sample Name','[Arabinose]'),'_') 
  
  # Factorise the sample name in the order for plotting
  fl$Results$`Sample Name` <- fl$Results$`Sample Name` %>% factor(levels = unique(.[sample_order]))
  fl$Results$`Primer pair` <- fl$Results$`Primer pair` %>% factor(levels = unique(.[sample_order])) # Factorise the primer pairs
  fl$Results$`[Arabinose]` <- fl$Results$`[Arabinose]` %>% factor(levels = unique(.[sample_order])) # [Arabinose]
  
  # plot the Tm of multiple peaks in melting curve ; Graph will now show
  
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- fl$Results %>% select(`Sample Name`, `[Arabinose]`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`, -`[Arabinose]`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `[Arabinose]`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting')) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    flf <- fl$Results %>% filter(`Target Name` == 'fMHT') %>% mutate(`Copy #` = 10^( (CT - f_intercept)/f_slope ) )  
    flr <- fl$Results %>% filter(`Target Name` == 'rMHT') %>% mutate(`Copy #` = 10^( (CT - r_intercept)/r_slope ) )  
    fljoin <- bind_rows(flf, flr)
    
    plt <- fljoin %>% ggplot(aes(x = `[Arabinose]`, y = `Copy #`, color = `Sample Name`)) +   # plotting
      scale_y_log10(  # logscale for y axis with tick marks
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
  } else plt <- fl$Results %>% ggplot(aes(x = `[Arabinose]`, y = CT, color = `Sample Name`))+ ylab(expression(C[q]))   
    
  # plot the CT mean along with replicates
  plt <- plt + geom_point(size = 1, show.legend = T) +
    # geom_line(aes(x = `[Arabinose]`, y = `Ct Mean`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
  print(plt)
  
}


# Custom plots (Transformed Assay data; Plots copy #; 'Sample Name'_variable 'primer pair') ----


# Save plots manually - copy paste this command to console
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)