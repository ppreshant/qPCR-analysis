# playing with rdml files
# Prashant K ; 17/jan/2022

# load modules ----
library(RDML)
source('./0-general_functions_main.R') # Source the general_functions file before running this



# User inputs ----
file_rdml <- 'S019_25-11-19_v152_dyefixed'

dirpath <- 'RDML_files/'

flpath <- str_c(dirpath, file_rdml, '.rdml')

# Load file ----
fl <- RDML$new(flpath)

# get fluorescence data (amplification step)
fl.data <- fl$GetFData()


# validation ----

# check if dye is recognized : fixed v152 files
fl$dye$SYBR$id

# check that dyeID is attached to one target
fl$target$flipped$dyeId


# write RDML ----
fl$AsXML(str_c(dirpath, file_rdml, "_AsXML.rdml")) # This file can't be read by rdmlpython, why?

# try reading back in : WORKS
fl2 <- RDML$new(str_c(dirpath, file_rdml, "_AsXML.rdml"))

fl2$target$flipped
