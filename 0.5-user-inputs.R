# 0.5-user_inputs.R

# Provide the file name to read and other variables that change for each qPCR run.

# File name ----

flnm <- 'S091-marine-withsort-5-5-24' # mention the file name without the "-linreg" or "-processed" suffixes

# Default standards : When using same standards for multiple plates (not considered best practice..but I do it anyways)
default_std.to.retrieve <-  'Std30' # if the file name doesn't hold any std curve, it will default to this (Std30 for RAM)
skip.std.curves_already.exist <- TRUE # If TRUE, will retrieve std curve data from the google sheet, instead of file. 

# Options ----

# get templates from google sheet or from data excel file
template_source <- 'googlesheet' # googlesheet/excel = parse through a list of templates in the respective formats and
# get the template with the matching qxx ID. 'excel' looks for the file 'excel files/Plate layouts.xlsx'


## qPCR options ----

# Linregpcr options
run_linregpcr <- TRUE

# Traditional Cq analysis - options
plot_mode <-  'absolute_quantification' # Options : ('absolute_quantification' or 'raw_quantification'); 
# absolute_quantification = Will calculate copy #'s based on y_intercept and slope from standard curve.. 
# ..calculated or gathered from old std curves 
# raw_quantification = Cq values are plotted

# Standard curve options
# This will pick the standard curve within the filename '.._Stdx_..', or use default if not found
# This proactive user setting prevents duplicates being processed into the sheet when rerunning script
# FALSE implies stds will be processed from file - if they exist
dilutions_to_truncate <- 3 # indicate the last n dilutions to be trimmed away before making std curve, default 0

# Options for default standards
force.use_default.std <- TRUE # forces the use of default std from google sheet, if skip.std.curves_already.exist is TRUE


## ddPCR options ----


# flag wells to remove from analysis
# wells_to_remove <- 'none' # 'none' or a regex of well names to remove

# Enable the following line to manually flag wells to remove from plots and analysis (use A1, not A01)
wells_to_remove <-
  c(
    # S091
    'A1', # few total droplets
    'F2' # Smear in positive

    # S091b
    # 'A7', 'B7', 'C7', 'F8', 'G8', 'H8', # different master mix, exclude NTC (C07)
    # 'F12' # wonky plasmid amplification, exclude NTC
  )

# older wells..
# 'A1', 'B1', 'A2', 'B2', 'H2', 'D4', 'H4', # wells with no chromosome amplification
# 'B6', 'G8', 'F9', 'F10', # wells with no chromosome amplification (outliers ~ Rd)


# flag categories to remove from analysis
categories_to_remove <- '^Ec' # 'none' or a regex of categories to remove
  

# TODO: add a LOD calculation (positive droplets) and remove stuff below it


# custom options 

# highly customized for S091b-marine-higher conc
assay_var_order <- c('J23100', 
                 'high', 'med', 'low',
                 'high-growth', 'med-growth', 'low-growth',
                 'WT', 
                 'blank', 'NTC') # order the assay variables


# error check ----
if(stringr::str_detect(flnm, '-linreg|-processed')) 
{stop('Invalid filename (flnm) in 0.5-user-inputs.R,
  flnm contains suffixes -linreg or -processed')}

# Calculation ----

date_regex <- '[:digit:]*-[:digit:]*-[:digit:]*' # Date regex

# title_name : appears on the html file name and header of selected plots, change as required
base_title_name <- stringr::str_remove(flnm, stringr::str_c("_", date_regex)) 

