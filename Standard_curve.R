# Plots standard curve for a given file (make sure it has quantities and targets and it is standard curve data)

# choose file name, in the same directory as Rproject
flnm <- 'excel files/Std5.xls'  

fl <- readqpcr(flnm) # read file
# fl$Results <- fl$Results %>% filter(`Ct Mean` < 27) # filtering only standard curve within the linear range

plt <- plotstdcurve(fl,'qPCR Standard curve 4', 'log(Copy #)') # plot standard curve

# # Extract the names of the targets in use
# targets_used <- fl$Results %>% filter(Task == 'STANDARD') %>% pull(`Target Name`) %>% unique(.)  

# Isolating standard curve variables (Quantity,CT) of the different targets into groups
standard_curve_vars <- fl$Results %>% filter(Task == 'STANDARD')  %>% select(Quantity, CT,`Target Name`) %>% group_by(`Target Name`) # select required columns and group

# Apply linear regression and find the model fitting results (equation and slope, R2 values) for each target
std_table <- standard_curve_vars %>% do(., equation = lm_eqn(.), params = lm_eqn(., trig = 'coeff'), dat = .[[1]] )

# function to add labels to plot - linear regression equation
function(std_table_entry){
  
}

# if both targets are present
if(nrow(trev) > 0)
  
{  # labelling the plot with the regression line equation and R2 value  
  plt2 <- plt + geom_label_repel(data = trev[3,], label = lm_eqn(trev), parse = TRUE, show.legend = F, nudge_x = -4, nudge_y = -1) +
    geom_label_repel(data = tfw[3,], label = lm_eqn(tfw), parse = TRUE, show.legend = F, nudge_x = -4, nudge_y = -6)
  
  print(plt2)
  
  # save plot manually
  
  # make a data table output
  crev <- lm_eqn(trev,'coeff'); cfw <- lm_eqn(tfw,'coeff')
  
  efficiency_table <- tibble(Slope = c(crev[1],cfw[1]), Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = c(crev[2],cfw[2])) %>% round(2)
  rownames(efficiency_table) <- c('rMHT','fMHT')
  View(efficiency_table)
  
} else {  # if rMHT target is absent in the data
  
  # labelling the plot with the regression line equation and R2 value  
  plt2 <- plt + geom_label_repel(data = tfw[3,], label = lm_eqn(tfw), parse = TRUE, show.legend = F, nudge_x = -4, nudge_y = -6)
  
  print(plt2)
  
  # save plot manually
  
  # make a data table output
  cfw <- lm_eqn(tfw,'coeff')
  
  efficiency_table <- tibble(Slope = cfw[1], Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = cfw[2]) %>% round(2)
  rownames(efficiency_table) <- c('fMHT')
  View(efficiency_table)
}