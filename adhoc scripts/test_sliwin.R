# adhoc inserting sliwin/pcrfit (linregpcr like) workflow into analysis.R

library(qpcR)

# Load a dataset into file by running analysis.R until line 75

plate_template <- get_template_for(flnm, sheeturls$plate_layouts_PK) %>% # read samplenames from googlesheets
  
  # Parsing sample names from the googlesheet table
  separate(`Sample_name_bulk`, # Split the components of the sample name bulk by delimiters ('-', '_', '.')
           c('Target_name', 
             'Sample_category',
             'assay_variable',
             'biological_replicates'),
           sep = '-|_|\\.') %>% 
  
  
  mutate(across('assay_variable', as.character)) %>% # useful when plasmid numbers are provided
  mutate(across('biological_replicates', ~str_replace_na(., '')) ) # if no replicates are provided, puts a blank string ('')
  

# get the amplification data (Rn) from the qPCR excel file
qpcr_amplification.data <- fl$`Amplification Data` %>% 
  left_join(plate_template) # join the samples names from the spreadsheet

# quick visual check
ggplot(qpcr_amplification.data,
       aes(x = Cycle,
           y = Rn,
           colour = Target_name)) +
  geom_line() + 
  facet_grid(~ assay_variable)
  

# trying sliwin

# data should be in wide format before doing this fitting with each sample's fluor in one column? 
qpcr_fit <- pcrfit(qpcr_amplification.data,
                         cyc = 3,
                         fluo = 5)


qpcr_sliwin <- sliwin(qpcr_fit)
# .Error in GRID[i, 1]:GRID[i, 4] : NA/NaN argument -- Is this expecting 4 channels of fluorescence?