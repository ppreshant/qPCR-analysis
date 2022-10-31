# Merge RDMLs - for combined analysis when negative samples are run separately
# Change one of the files using RDML ninja from Run001 to Run002

# load modules ----
library(RDML) # RDML library
library(tidyverse) # general data structures library


# User inputs ----
base_rdml.file <- 'q16c_Uxx-positives_15-10-21fixpython'
base_dirpath <- 'RDML_files/results/'

toadd_rdml.file <- 'q016b_Uxx-negatives_15-10-21_editformerge'
toadd_dirpath <- 'RDML_files/Archive/'

merged_file_name <- 'q16_Uxx_15-10-21-merged'

# Load file ----
basefl <- RDML$new(str_c(base_dirpath, base_rdml.file, '.rdml'))
toaddfl <- RDML$new(str_c(toadd_dirpath, toadd_rdml.file, '.rdml'))


# Merge ----
mergedfl <- MergeRDMLs(list(basefl, toaddfl))

# write RDML ----
mergedfl$AsXML(str_c(base_dirpath, merged_file_name, "_AsXML.rdml")) # write merged file, can be read in online RDML-RunView
# Online RDML-Linregpcr : Corrupted RDML file was repaired (no rdml_data.xml) - please save fixed RDML file!
# The online tool still analyzes only 1 run at a time. Such as waste of time ... :(

# This file can't be read by rdmlpython, why?
# error: No rdml_data.xml in compressed RDML file found.




# optional steps
# Linregpcr on merged RDML ----


flnm <- str_c(merged_file_name, '_AsXML') # mention the file name without the "-linreg" or "-processed" suffixes

base_title_name <- 'q16c_Uxx-merged' # This will be the title name of the plots and the name of the html file

source('scripts_general_fns/17-run_linregpcr.R') 

# run linregpcr
run_linregpcr(flnm)
