# qPCR related functions (good to copy to the general qPCR project)

# qPCR related  ----

# lookup table of primer pairs and respective targets - not useful as same information is in `Target` already
primer_table <- c('q1-3' = 'Flipped', 'q4-5' = 'Flipped', 
                  'q5-11' = 'Unflipped', 'q1-2' = 'Unflipped', 'q9-10' = 'Unflipped', 'q4-2' = 'Unflipped', 'q2-4' = 'Unflipped', 'q5-11' = 'Unflipped', 'q6-7' = 'Unflipped',
                  'q12-13' = 'Backbone')



# function to back-calculate CT using standard curve parameters
absolute_backcalc <- function(df, std_par)
{
  target_current <- df$Target %>% unique()
  std_current <- std_par %>% filter(str_detect(target_current, Target))
  
  df %>% mutate(`Copy #` = 10^( (CT - std_current$y_intercept)/std_current$Slope) )
}



# Plot Standard curve
plotstdcurve <- function(results_qpcr, plttitle, xlabel)
{
  plt <- results_qpcr %>% 
    ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target`) + geom_point() +
    theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle(plttitle) + xlab(xlabel) + ylab(expression(C[q])) +
    stat_smooth(data = filter(results_qpcr, Task == 'STANDARD'), method ="lm", se = F) # plots linear regression line
}


# output the difference between consecutive CT values
optional1 <- function()
{
  tsumrev <- trev %>% group_by(`Sample_name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
  diff(tsumrev$CT) %>% round(2)
}


# Tm plots ----


# plotting functions for Melting temperature
# Obsolete: Were used in small_scale mode and old_assay mode

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample_name`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample_name`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample_name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample_name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

