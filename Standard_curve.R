# Plots standard curve for a given file (make sure it has quantities and targets and it is standard curve data)

# choose file name, in the same directory as Rproject
flnm <- 'excel files/Std5.xls'  

fl <- readqpcr(flnm) # read file
# fl$Results <- fl$Results %>% filter(`Ct Mean` < 27) # filtering only standard curve within the linear range

plt <- plotstdcurve(fl,'qPCR Standard curve 5', 'log(Copy #)') # plot standard curve

# # Extract the names of the targets in use
# targets_used <- fl$Results %>% filter(Task == 'STANDARD') %>% pull(`Target Name`) %>% unique(.)  

# Isolating standard curve variables (Quantity,CT) of the different targets into groups
standard_curve_vars <- fl$Results %>% filter(Task == 'STANDARD')  %>% select(Quantity, CT,`Target Name`) %>% group_by(`Target Name`) # select required columns and group

# Apply linear regression and find the model fitting results (equation and slope, R2 values) for each target
std_table <- standard_curve_vars %>% do(., equation = lm_eqn(.), params = lm_eqn(., trig = 'coeff'), dat = .[1,] )
std_table$params %<>% bind_rows() # make parameters and data into tibbles : do function makes lists
std_table$dat %<>% bind_rows()  

# Add labels to plot - linear regression equation
plt + geom_label_repel(data = std_table$dat, label = std_table$equation, parse = TRUE, show.legend = F)

# processing linear regression out
efficiency_table <- tibble(Slope = std_table$params %>% pull(slope), Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = std_table$params %>% pull(r_square) %>% round(2))
rownames(efficiency_table) <- std_table$dat$`Target Name`
View(efficiency_table)
