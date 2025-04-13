# general functions to read files
# reading files and manipulating columns

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  # bring the first sheet to count the number of rows to skip
  setup_sheet <- flnm %>%  
    readxl::read_excel(., sheet = 'Sample Setup', col_names = FALSE) %>% # read any sheet (Sample Setup/Results)
    select(1:3) # retain only the first 3 columns
    
  to.skip  <- setup_sheet %>% 
    pull(1) %>% # select first column only
    {which(. == 'Well') - 1} # check where "Well" appears, and subtract 1
  
  fl <- flnm %>%  
    readxl::excel_sheets() %>% 
    set_names(.,.) %>% 
    map(readxl::read_excel, path = flnm, skip = to.skip) # change skip # if channels are not calibrated
  
  # future: put an error here if column names are not read properly (maybe look for well position and one other column..)
  
  # Check if probes are being used - look for TAQMAN chemistry => multiplexing
  fl['chemistry_type'] <- setup_sheet %>% filter(str_detect(`...1`, 'Chemistry')) %>% pull(2)

  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}


# Gets the 96 well layout with template names matching the experiment ID from filename in a google sheet
get_template_for <- function(bait, sheet_url = sheeturls$plate_layouts_PK)
{ # Looking for WWx or Stdx - example WW21 or Std7 within the filename; 
  # Assumes plate spans from row B to N (1 row below the matching ID)
  
  # get template from google sheets or excel file
  
  # Finding the plate to be read
  plate_names_row <- if(template_source == 'googlesheet') # googlesheet vs excel options
  
    {googlesheets4::read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = 'C:C', col_types = 'c') } else {
      
      readxl::read_excel(path = 'excel files/Plate layouts.xlsx', range = readxl::cell_cols('C:C'), col_types = 'text')
    } 
  
  # get the plate ID from the full filename / bait
  plate_id <- str_extract(bait, plate_id_regex) # get the first part containing the plate ID
  
  # find the row number that matches the plate_id
  m_row <- plate_names_row %>% unlist() %>% as.character() %>% 
    {str_extract(., plate_id_regex) == plate_id} %>% # select the q0xya digits/letters part of the filename
    which() + 1 # extract the row with an exact match to the plate_id and add 1 to get the plate contents
  
  # Eror message and terminate if plate ID is not unique
  if(length(m_row) > 1) {stop( str_c('Plate ID :', plate_id, 
                                     'of filename :', bait, 
                                     'repeats in', paste0(m_row, collapse = ' & '), 
                                     'row numbers. Please fix and re-run the script', sep = ' ')) 
    
    # or if no matching plate is found
    } else if(!length(m_row)) stop( str_c('Plate ID of :', bait, 'does not match anything on the plate layout. 
    Please fix and re-run the script', sep = ' '))
  
  # make the full range for the specific plate
  range_to_get <- str_c('B', m_row + 1, ':N', m_row + 9) 
  
  # read the template corresponding to the file name -- from the range determined above
  plate_template_raw <- 
    if(template_source == 'googlesheet') # googlesheet vs excel options
    
      {googlesheets4::read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = range_to_get)} else {
      
      readxl::read_excel(path = 'excel files/Plate layouts.xlsx', range = range_to_get)
    }
  
  # Convert the 96 well into a single column, alongside the Well position (column wise order)
  plate_template <- read_plate_to_column(plate_template_raw, 'Sample_name_bulk')
  
}



#' gets plate layout from 96 well format and parse it into individual columns
#' @param flnm filename to get the plate layout from
#' @param read_target_name logical, whether to read the target name column
#' @return a tibble with the plate layout in long format (each well in a row)
get_and_parse_plate_layout <- function(flnm, read_target_name = TRUE)
{
  
  plate_template <- get_template_for(flnm, sheeturls$plate_layouts_PK) %>% # read sample names from googlesheets
    
    # Parsing sample names from the google sheet table  
    separate_wider_delim(
      cols = `Sample_name_bulk`, # Split the components of the sample name bulk by delimiter ('_')
      delim = '_', 
      names = c(if(read_target_name) 'Target_name', # if target name is to be read, else NULL
                'Sample_category',
                'assay_variable',
                'biological_replicates'),
      too_few = 'align_start') %>% # fill missing values starting from the right (replicates)
    
    
    mutate(across('assay_variable', as.character)) %>% # useful when plasmid numbers are provided, will convert to text
    
    # fix biological replicates
    {if(pull(., biological_replicates) %>% is.na %>% all) { # when biological replicates are not provided..
      
      group_by(., across(!any_of(c('biological_replicates', 'Well Position')))) %>% 
        mutate(biological_replicates = row_number()) %>% # ..make them up within each group (in column wise order)
        ungroup()
      
      # for few samples missing the replicates, puts a blank string ('')
      } else mutate(., across('biological_replicates', ~str_replace_na(., '')) ) 
    } %>% 
    
    # account for water controls for lysis
    clean_up_water_wells() # changes `Water` to `control` sample_category, assay_var 'water' and numbers-> replicates
    
  
  
}


#' Read a bunch of processed data sheets and join them, including the qxx as run_ID
#' @param .flnms vector of strings with name of file, without the '-processed.csv' 
get_processed_datasets <- function(.flnms)
{
  .df <- map_dfr(.flnms, # for each file, 
                 
                   # read excel file exported by Quantstudio
                 ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv')) %>%  
                   mutate(run_ID = str_extract(.x, plate_id_regex)) # add the run_ID from the filename
  )
  
  
}