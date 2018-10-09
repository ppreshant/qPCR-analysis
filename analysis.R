# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# Source the general_functions file before running this
# The last part of the script contains important codes from the console that were used for specific situations : These will be commented


# choose file name, title for plots and experiment mode (file name starts in the same directory as Rproject)----
flnm <- 'excel files/Int_Assay1 MHT.xls'  
title_name <-'Integrase induction'
experiment_mode <- 'assay' # 'small_scale'  # 'assay' ; 'custom'

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
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`)
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(fl)
{ 
  plttm <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`)
}


# Input the data ----

# reading in file and transformations
fl <- readqpcr(flnm) # read excel file exported by Quantstudio

sample_order = order_columnwise(fl) # this order is columnwise (data is shown row-wise) => useful for plotting column wise order

# Transform and Plots for typical troubleshooting or small scale assays (facetted by primer names)----

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
    ggtitle(title_name) + facet_grid(~`Primer pair`)
  
  print(plt)
}


# Plots for Assays ----
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
  
  # plot the CT mean along with replicates
  plt <- fl$Results %>% ggplot(aes(x = `[Arabinose]`, y = CT, color = `Sample Name`)) + geom_point(size = 1, show.legend = T) +
    geom_line(aes(x = `[Arabinose]`, y = `Ct Mean`), show.legend = T) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(title_name) + facet_wrap(~`Sample Name`, scales = 'free_x')
  
  print(plt)
  
}


# Custom codes ----

# Save plots manually - copy paste this command to console
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)