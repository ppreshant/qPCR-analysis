# general functions to read files
# reading files and manipulating columns

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = 41) # change skip # if channels are not calibrated
  
  # put an error here if column names are not read properly (maybe look for well position and one other column..)
  
  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}


# Gets the 96 well layout with template names matching the experiment ID from filename in a google sheet
get_template_for <- function(bait, sheet_url = sheeturls$plate_layouts_PK)
{ # Looking for WWx or Stdx - example WW21 or Std7 within the filename; Assumes plate spans from row B to N (1 row below the matching ID)
  
  # Finding the plate to be read
  plate_names_row <- read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = 'C:C', col_types = 'c')
  m_row <- plate_names_row %>% unlist() %>% as.character() %>% 
    # find the row with standard beginnings matching the filename
    str_detect(., str_c('^', bait %>% str_match('^(q)[:digit:]*') %>% .[1]) ) %>% 
    # str_detect(., bait) %>% 
    which() + 1
  range_to_get <- str_c('B', m_row + 1, ':N', m_row + 9)
  
  # Eror message and terminate if plate ID is not unique
  if(length(m_row) > 1) stop( str_c('Plate ID of :', bait, 'repeats in', paste0(m_row, collapse = ' & '), 'row numbers. Please fix and re-run the script', sep = ' '))
  
  # read the template corresponding to the file name
  plate_template_raw <- read_sheet(sheet_url, sheet = 'qPCR plate layouts', range = range_to_get)
  
  # Convert the 96 well into a single column, alongside the Well position
  plate_template <- read_plate_to_column(plate_template_raw, 'Sample_name_bulk') # convert plate template (Sample_names) into a single vector, columnwise
  
}