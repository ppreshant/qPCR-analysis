# Plots standard curve for a given file (make sure it has quantities and targets and task set to standard curves (S))

source('./general_functions.R') # Source the general_functions file before running this

# User inputs ----
# choose file name, in the same directory as Rproject
flnm <- 'Std5'  # set the filename
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
target_variable_set_in <- 'Import' # 'GUI' if targets were entered in quantstudio GUI vs 'Import' if entered by loading excel template into quantstudio

# Data input, processing ----
fl <- readqpcr(flpath) # read file

if(target_variable_set_in == 'GUI') fl$Results %<>% rename('Target' = `Target Name`) # All subsequent calls use 'Target' so change the name 

# optional filtering to remove low concentration points in standard curve
# fl$Results <- fl$Results %>% filter(`Quantity` > 1e4) # filtering only standard curve within the linear range

# Plotting ----
plt <- plotstdcurve(fl,'qPCR Standard curve 5', 'log(Copy #)') # plot standard curve

# Isolating standard curve variables (Quantity,CT) of the different targets into groups
standard_curve_vars <- fl$Results %>% filter(Task == 'STANDARD')  %>% select(Quantity, CT,`Target`) %>% group_by(`Target`) # select required columns and group

# linear regression ----
# Apply linear regression and find the model fitting results (equation and slope, R2 values) for each target
std_table <- standard_curve_vars %>% do(., equation = lm_eqn(.), params = lm_eqn(., trig = 'coeff'), dat = .[1,] ) # "do" applies functions to each group of the data
std_table$params %<>% bind_rows() # Convert parameters and data into tibbles : "do" function makes htem lists
std_table$dat %<>% bind_rows()  

# Add labels to plot - linear regression equation
plt + geom_label_repel(data = std_table$dat, label = std_table$equation, parse = TRUE, show.legend = F, nudge_x = -3, force = 10)

# Data output ----
# processing linear regression out
efficiency_table <- tibble(Slope = std_table$params %>% pull(slope), Efficiency = 10^(-1/Slope), '% Efficiency' = (Efficiency -1)*100 , 'R-square' = std_table$params %>% pull(r_square) %>% round(2))
rownames(efficiency_table) <- std_table$dat$`Target`
View(efficiency_table)
