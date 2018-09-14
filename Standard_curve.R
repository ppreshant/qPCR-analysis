# Plots standard curve for a given file (make sure it has quantities and it is standard curve data)

# choose file name, in the same directory as Rproject
flnm <- 'excel files/2018-07-11_MHT Std3.xls'  

fl <- readqpcr(flnm) # read file
plt <- plotstdcurve(fl,'qPCR Standard curve : Run 2') # plot standard curve

# Isolating standard curves of the two different targets
tfw <- fl$Results %>% filter(Task == 'STANDARD', `Target Name` == 'fMHT')
trev <- fl$Results %>% filter(Task == 'STANDARD', `Target Name` == 'rMHT') 

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