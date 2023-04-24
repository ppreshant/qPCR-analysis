# 0.5-user_inputs.R

# Provide the file name to read and other variables that change for each qPCR run.

# File name ----

flnm <- 'q48a_S071_22-4-23' # mention the file name without the "-linreg" or "-processed" suffixes

# Options ----

# get templates from google sheet or from data excel file
template_source <- 'googlesheet' # googlesheet/excel = parse through a list of templates in the respective formats and
# get the template with the matching qxx ID. 'excel' looks for the file 'excel files/Plate layouts.xlsx'

# Linregpcr options
run_linregpcr <- TRUE

# Traditional Cq analysis - options
plot_mode <-  'raw_quantification' # Options : ('absolute_quantification' or 'raw_quantification'); 
# absolute_quantification = Will calculate copy #'s based on y_intercept and slope from standard curve.. 
# ..calculated or gathered from old std curves 
# raw_quantification = Cq values are plotted

# Standard curve options
skip.std.curves_already.exist <- TRUE # If TRUE, will retrieve std curve data from the google sheet
# This will pick the standard curve within the filename '.._Stdx_..', or use default if not found
# This pro-active user setting prevents duplicates being processed into the sheet when rerunning script
# FALSE implies stds will be processed from file
dilutions_to_truncate <- 0 # indicate the last n dilutions to be trimmed away before making std curve, default 0

# Default standards : When using same standards for multiple plates (not considered best practice..)
default_std.to.retrieve <-  'Std30' # if the file name doesn't hold any std curve, it will default to this
force.use_default.std <- TRUE # forces the use of default std from google sheet, instead of from current file


# error check ----
if(stringr::str_detect(flnm, '-linreg|-processed')) 
{stop('Invalid filename (flnm) in 0.5-user-inputs.R,
  flnm contains suffixes -linreg or -processed')}

# Calculation ----

date_regex <- '[:digit:]*-[:digit:]*-[:digit:]*' # Date regex

# title_name : appears on the html file name and header of selected plots, change as required
base_title_name <- stringr::str_remove(flnm, stringr::str_c("_", date_regex)) 

