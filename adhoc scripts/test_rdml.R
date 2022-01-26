# playing with rdml files
# Prashant K ; 17/jan/2022

# load modules
library(RDML)
source('./0-general_functions_main.R') # Source the general_functions file before running this



# User inputs
file_rdml <- 'S019 25-11-19.rdml'

dirpath <- 'excel files/'

flpath <- str_c(dirpath, file_rdml)

# Load file
fl <- RDML$new(flpath)
fl.data <- fl$GetFData()

fl$target

