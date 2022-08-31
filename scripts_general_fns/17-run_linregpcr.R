# 17-run_linregpcr.R

#' R-reticulate interface to run linregPCR through python
#' 
#' Ensure the correct filename is chosen in 0.5-user-inputs.R file

run_linregpcr <- function(rdml_file_name)
{
  # Load package
  library(reticulate)
  
  # set python
  
  # try to add python to path
  Sys.setenv(RETICULATE_PYTHON = 'C:\\ProgramData\\Miniconda3\\python.exe')
  
  # import linregpcr python interfacing function as a module
  linregpcr <- import('scripts_general_fns.g14_linregpcr_analyze')
  
  # Run the function with the rdml file name
  linregpcr$linregpcr_analyze(rdml_file_name)
  
  
}