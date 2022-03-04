# FUnctions to load qPCR data and manipulate it. The functions can be called from another R file

# read in excel file (.xls) of qPCR exported from Quantstudio 3 
  # Make sure to include raw data as well

# Library calling  ----
# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel); library(googlesheets4); library(rlang); library(plotly)
# library(lubridate); # for advanced date operations

# python integration
# reticulate::use_python('C:/ProgramData/Miniconda3/python.exe') # Change path to your python installation
# use_condaenv(condaenv = 'base', conda = 'C:/ProgramData/Miniconda3/Scripts/conda')

# R-ladies baltimore recommends this - https://www.youtube.com/watch?v=U3ByGh8RmSc
# Sys.setenv(RETICULATE_PYTHON = 'C:\\ProgramData\\Miniconda3\\python.exe')

# one-time configuration of 'reticulate' python interface
# Have your own miniconda/anaconda with environments you need to load
# Sys.setenv(RETICULATE_MINICONDA_PATH = "/path/to/r-miniconda")  # https://github.com/rstudio/reticulate/issues/745
  

# google sheets ----
sheeturls <- list(plate_layouts_PK = 'https://docs.google.com/spreadsheets/d/1RffyflHCQ_GzlRHbeH3bAkiYo4zNlnFWx4FXo7xkUt8/edit#gid=0')

# Convenience function for plotting directory and .png suffix for adhoc plots
plot_as <- function(plt_name)  str_c('qPCR analysis/Archive/', plt_name, '.png')

# dummy data  ---- 
# or test data for testing simple functions 

# dummy test tibble
a <- tibble(a1 = 1:6, a2 = 6:1, a3 = rep(c('a', 'b'),3), a4 = a2 ^2)

# expression to test on plotting
y_namr_test <- list( 'a2' = expression(paste('this is a ', mu, 'L')),
                     'a4' = expression(paste('super large ', sigma, 'L')))

# test ggplot
a_plt <- ggplot(a, aes(a1, a2, colour = a3)) + 
  geom_point() + 
  geom_line() + 
  ylab(y_namr_test[['a4']])



# calling more funs ----

list_of_general_functions <- c("10-process_qpcr.R",
                               "1-reading_files_funs.R",
                               "2-tibble_columns_funs.R",
                               "3-obsolete_arcane_funs.R",
                               "4-qPCR_specific_funs.R",
                               "5-mathematical_fitting_funs.R",
                               "6-formatting_plot_funs.R",
                               "7-COVID specific_writing_funs.R",
                               "8-plot_mean_sd_jitter.R",
                               "9-plot_scatter.R",
                               '12-plot_facetted_assay.R',
                               'g.13-std_curve_processing.R') # "0-old_general_functions.R"

# Source all the functions listed above
map(str_c('./scripts_general_fns/', list_of_general_functions),
    source)






