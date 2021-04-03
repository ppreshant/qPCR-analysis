# plotting functions for melting temperature (number 11)

# All tms ----
# plots all the Tm's if samples have multiple peaks in the melting curve
# Works only for small scale
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  Tm_data <- results_relevant %>% select(`Sample Name`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`, -`Primer pair`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- Tm_data %>% ggplot(.) + aes(x = `Sample Name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

# First Tm ----
# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}
