# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# Source the general_functions file before running this
# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

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


# Input the data + transformations ----
# choose file name, starts in the same directory as Rproject
flnm <- 'excel files/Tr12 MHT.xls'  
title_name <-'Water troubleshooting: Tr12' 

fl <- readqpcr(flnm) # read excel file exported by Quantstudio

# Separate primer pair from sample name and make factors in the right order (same order as the plate setup)
fl <- separate_primer_names(fl)


# Plots for typical troubleshooting or small scale assays (facetted by primer names)----

# plot the Tm ; Graph will now show
plttm2 <- plotalltms(fl) # plots tms of multiple peaks in melting curve

# plot the CT mean along with replicates
plt <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 1, show.legend = T) +
  geom_boxplot(aes(x = `Sample Name`, y = `Ct Mean`), show.legend = T) +
  theme_classic() + scale_color_brewer(palette="Set1") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
  ggtitle(title_name) + facet_grid(~`Primer pair`)

print(plt)

# Save plots manually - copy paste this command to console ----
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)

# Plots for Assays ----


# Custom codes ----
