# general functions to read files
# reading files and manipulating columns

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  # bring the first sheet to count the number of rows to skip
  setup_sheet <- flnm %>%  
    read_excel(., sheet = 'Sample Setup', col_names = FALSE) %>% # read any sheet (Sample Setup/Results)
    select(1:3) # retain only the first 3 columns
    
  to.skip  <- setup_sheet %>% 
    pull(1) %>% # select first column only
    {which(. == 'Well') - 1} # check where "Well" appears, and subtract 1
  
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = to.skip) # change skip # if channels are not calibrated
  
  # future: put an error here if column names are not read properly (maybe look for well position and one other column..)
  
  # Check if probes are being used - look for TAQMAN chemistry => multiplexing
  fl['chemistry_type'] <- setup_sheet %>% filter(str_detect(`...1`, 'Chemistry')) %>% pull(2)

  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}


# Gets the 96 well layout with template names matching the experiment ID from filename in a google sheet
get_template_for <- function(bait, sheet_url = sheeturls$plate_layouts_PK)
{ # Looking for WWx or Stdx - example WW21 or Std7 within the filename; Assumes plate spans from row B to N (1 row below the matching ID)
  
  # get template from google sheets or excel file
  
  # Finding the plate to be read
  plate_names_row <- if(template_source == 'googlesheet') # googlesheet vs excel options
  
    {read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = 'C:C', col_types = 'c') } else {
      
      readxl::read_excel(path = 'excel files/Plate layouts.xlsx', range = cell_cols('C:C'), col_types = 'text')
    } 
  
  
  m_row <- plate_names_row %>% unlist() %>% as.character() %>% 
    # find the row with standard beginnings matching the filename
    str_detect(., str_c('^', bait %>% str_match('^(q)[:digit:]*') %>% .[1]) ) %>% # select the q0xy digits part of the file
    # str_detect(., bait) %>% 
    which() + 1
  range_to_get <- str_c('B', m_row + 1, ':N', m_row + 9)
  
  # Eror message and terminate if plate ID is not unique
  if(length(m_row) > 1) {stop( str_c('Plate ID of :', bait, 'repeats in', paste0(m_row, collapse = ' & '), 
                                     'row numbers. Please fix and re-run the script', sep = ' ')) 
    # or if no matching plate is found
    } else if(!length(m_row)) stop( str_c('Plate ID of :', bait, 'does not match anything on the plate layout. 
    Please fix and re-run the script', sep = ' '))
  
  # read the template corresponding to the file name -- from the range determined above
  plate_template_raw <- 
    if(template_source == 'googlesheet') # googlesheet vs excel options
    
      {read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = range_to_get)} else {
      
      readxl::read_excel(path = 'excel files/Plate layouts.xlsx', range = range_to_get)
    }
  
  # Convert the 96 well into a single column, alongside the Well position
  plate_template <- read_plate_to_column(plate_template_raw, 'Sample_name_bulk') # convert plate template (Sample_names) into a single vector, columnwise
  
}



# gets plate layout from 96 well format and parse it into individual columns
get_and_parse_plate_layout <- function(flnm)
{
  
  plate_template <- get_template_for(flnm, sheeturls$plate_layouts_PK) %>% # read samplenames from googlesheets
    
    # Parsing sample names from the google sheet table  
    separate(`Sample_name_bulk`, # Split the components of the sample name bulk by delimiters ('-', '_', '.')
             c('Target_name', 
               'Sample_category',
               'assay_variable',
               'biological_replicates'),
             sep = '-|_|\\.') %>% 
    
    
    mutate(across('assay_variable', as.character)) %>% # useful when plasmid numbers are provided, will convert to text
    mutate(across('biological_replicates', ~str_replace_na(., '')) ) # if no replicates are provided, puts a blank string ('')
  
  
}


#' Read a bunch of processed data sheets and join them, including the qxx as run_ID
#' @param .flnms vector of strings with name of file, without the '-processed.csv' 
get_processed_datasets <- function(.flnms)
{
  .df <- map_dfr(.flnms, 
                 ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv')) %>%  # read excel file exported by Quantstudio
                   mutate(run_ID = str_extract(.x, 'q[:alnum:]*')) # add the run_ID from the filename
  )
  
  
}