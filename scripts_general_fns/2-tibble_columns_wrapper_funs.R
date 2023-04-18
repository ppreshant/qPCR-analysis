# Function wrappers for general operations with tibbles/data frames (manipulating columns)


# takes in a qPCR file (read in columnwise) and arranges elements rowise
columnwise_index <- function(fl)
{ # this orders the samples columnwise in the PCR plate or strip : order will give indices for A1,B1,C1 ... A2,B2,C2... (wherever samples exist)
  
  fl$Results$`Well Position` %>%  str_sub(2) %>% as.integer() %>% order() # Take the well position A1 etc., extract the number, read as integer (instead of char), order it (1,1,1.. 2,2,... 12,12...) and regurn indices 
}

# Convert the 96 well into a single column, alongside the Well position
read_plate_to_column <- function(data_tibble, val_name)
{ # transforms a plate reader table into a column (named after the top left cell, unless mentioned)
  # eliminates plate row,column numbering ; Select 1 row above the plate (even if it doesn't contain a label)
  
  val_name <- enquo(val_name)
  # colnames(data_tibble) <- data_tibble[1,] # set column names as the first row
  data_tibble[,] %>% gather(key = 'col_num', value = !!val_name, -`<>`) %>% rename(row_num = `<>`) %>% unite('Well Position', c('row_num', 'col_num'), sep = '') %>% drop_na()
}


# mutates a subset of data and returns a new array (works for multiple conditions)
mutate_when <- function(data, ...) 
{ # Source: Stackoverflow - https://stackoverflow.com/a/34170176/9049673
  
  dots <- eval(substitute(alist(...)))
  for (i in seq(1, length(dots), by = 2)) 
  {
    condition <- eval(dots[[i]], envir = data)
    mutations <- eval(dots[[i + 1]], envir = data[condition, , drop = FALSE])
    data[condition, names(mutations)] <- mutations
  }
  data
}


# mutates a subset of data and returns a new array (does multiple mutations on same condition)
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) 
{ # Source: Stackoverflow -  https://stackoverflow.com/a/34096575/9049673
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


# paste without NAs
paste_without_NAs <- function(.string1, .string2, .sep = "-")
{
  # if string2 is NA, it leaves string1 unchanged, else joins them as string1 *separater* string2
  if_else(is.na(.string2), as.character(.string1), paste(.string1, .string2, sep = .sep))
}

# Calculate the percentage of the B117 variant in the data after pivoting
# where is the function?




#' function to remove last n quantities in std curve data 

remove_last_n_dilutions <- function(.data, dilutions_to_truncate = 0)
{
  # Group data by each target_name 
  group_by(.data, Target_name) %>% 
    
    filter(!str_detect(Sample_category, 'Std') | # either the sample is not Std or satisfies below critereon
             
             Quantity %in% # check for Quantity present in the top (total - n) unique values for each target
             {filter(., str_detect(Sample_category, 'Std')) %>%  # only Stds
                 select(Target_name, Quantity) %>%  # select only relevant columns
                 unique() %>% 
                 slice_max(order_by = Quantity, n = -dilutions_to_truncate) %>% # top (total - n) unique values
                 pull(Quantity)
             }) 
  # Warning : imperfect matching of all quantity values without context of Target_name 
  # -- fails in rare case when multiple targets have identical concentration values 
  # that are shifted so one must be included and other excluded
  
}


#' Wrapper to change "Water" samples to control and retain their numbering for reference to cross-contamination
#' Samples named as 'T_Water_3' -- this way preserves the numbering as biological replicate, converts Water into control
#' and makes assay variable as 'water'. The final plots need to be annotated selectively to show water (S065_q43_plotting.R)

clean_up_water_wells <- function(.df, new_sample_category_for_Water = 'control')
{
  mutate_cond(.df,
              str_detect(Sample_category, 'Water'),  # clean up names for 'Water' samples only
              
              # apply below transformations to only the "Water" rows
              across(assay_variable, ~ str_replace(., 'spillH4', '11')), # change the name to a number -- expt specific
              
              biological_replicates = as.numeric(assay_variable), 
              across(contains('assay_var'), ~ 'water'),
              across(Sample_category, ~ str_replace(., 'Water', new_sample_category_for_Water))) # change "Water" to "control"
  
}
