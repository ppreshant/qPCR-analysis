# Testing reticulate to call python module from R : g14_linregpcr_analyze.py

# Load package ----
library(reticulate)

# set python ----

# try to add python to path
Sys.setenv(RETICULATE_PYTHON = 'C:\\ProgramData\\Miniconda3\\python.exe')
# R-ladies baltimore recommends this - https://www.youtube.com/watch?v=U3ByGh8RmSc

# check if python was discovered
py_config() # query the current python in use + list other discovered versions
# Error in Sys.setenv(PATH = new_path) : wrong length for argument


# Troubleshooting python connection ----

# more test : https://github.com/rstudio/reticulate/issues/1155
python <- "C:\\ProgramData\\Miniconda3/python.exe"
conda_info <- reticulate:::get_python_conda_info(python)
new_path <- reticulate:::conda_run(
  "python",
  c("-c", shQuote("import os; print(os.environ['PATH'])")),
  conda = conda_info$conda,
  envname = conda_info$root,
  stdout = TRUE
)
new_path




# Find the python that is in the system path
Sys.which('python')

# set Miniconda environment variable 

# -- otherwise use_miniconda() will say miniconda is not installed
Sys.setenv('RETICULATE_MINICONDA_PATH' = 'C:/ProgramData/Miniconda3/')
# one-time configuration of 'reticulate' python interface
# Have your own miniconda/anaconda with environments you need to load
# https://github.com/rstudio/reticulate/issues/745

# set conda environment -- does not work
use_miniconda('C:/ProgramData/Miniconda3/python.exe')

conda_list() # check if conda environments are detected
# use_condaenv(condaenv = 'base', conda = 'C:/ProgramData/Miniconda3/Scripts/conda')

use_condaenv('base')




use_python("C:/ProgramData/Miniconda3/python.exe") # Change path to your python installation
# Error in Sys.setenv(PATH = new_path) : wrong length for argument



# Testing ----

# Load a python library
os <- import("os")

# import as a module
linregpcr <- import('scripts_general_fns.g14_linregpcr_analyze.py')
# Error in Sys.setenv(PATH = new_path) : wrong length for argument

# source python script
source_python('scripts_general_fns/g14_linregpcr_analyze.py')



