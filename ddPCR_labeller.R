# Read in the ddPCR data file and attach the names from plate template
# Author: Prashant Kalvapalle;  June 30 2020

# User inputs ----


flnm <- 'dd.WW4_629_N1-N2'  # set the filename
template_volume <- 10 /22 * 20 # ul template volume per well of the ddPCR reaction

# Biobot ID sheet to draw sample names from
bb_sheets <- c('Week 12 (6/29)')

# Loading pre-reqisites ----


# Loading libraries, functions and user inputs
source('./general_functions.R') # Source the general_functions file

# Data input ----

templates_sheet <- 'https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118'

# Read the ddPCR master file

# Read in qPCR data and labels from plate template
fl <- read_sheet('https://docs.google.com/spreadsheets/d/1jdO_P9SZGezSTLiIARtSmA7qaUuX3wA-jCe7YiQ1sCI/edit#gid=0', sheet = flnm) # read excel file exported by Quantstudio
plate_template <- get_template_for(flnm, templates_sheet)  # Get the plate template matching file name, convert to 1 column 

# Polishing ----


# Load desired qPCR result sheet and columns
bring_results <- fl %>% # select only the results used for plotting, calculations etc. and arrange them according to sample order
  select(-Sample) %>% 
  mutate_at('Well', ~ str_replace(., '0', '')) %>% rename('Well Position' = Well) %>% 
  right_join(plate_template, by = 'Well Position') %>%  # Incorporate samples names from the google sheet by matching well position
  mutate_at('Target', ~str_replace(., 'BSRV', 'BRSV')) %>% 
  filter(!is.na(Target)) %>% 
  mutate('Copies/ul RNA' = CopiesPer20uLWell/ template_volume) %>% 
  select(`Sample Name`, `Copies/ul RNA`, everything())

# Bring WWTP names from google sheet: "Biobot Sample IDs"
biobot_lookup <- map_df(bb_sheets , ~read_sheet('https://docs.google.com/spreadsheets/d/1ghb_GjTS4yMFbzb65NskAlm-2Gb5M4SNYi4FHE4YVyI/edit#gid=677034958', sheet = .x, range = 'H:J')) %>% rename('Biobot ID' = matches('Biobot|Comments', ignore.case = T), 'WWTP' = contains('SYMBOL', ignore.case = T)) %>% mutate('Biobot ID' = str_remove(`Biobot ID`,'\\.'), WWTP = as.character(WWTP))

# polishing qPCR data - Make Biobot ID column clean
# isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
polished_results <- bring_results %>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% separate(`Sample Name`,c('Sample Name','Tube ID'),'_') %>% 
  mutate(`Tube ID` = if_else(`Sample Name` == 'NTC', '0', `Tube ID`)) %>% 
  separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F) %>%  # Separate out biological replicates 
  unite('Tube ID', c(assay_variable, biological_replicates), sep = '.', remove = F, na.rm = T) %>% # remaking Tube ID - removes spaces after 'dot'
  unite('Biobot ID', c(`Sample Name`, assay_variable), sep = '', remove = F) %>%

  mutate_at('assay_variable', as.character) %>% 
  mutate_at('biological_replicates', ~str_replace_na(., ''))
  


# regex based biobot ID cleaning
# qpcr_polished <- bring_results %>% mutate('Biobot ID' = str_match(`Sample Name`, '-.*_.*') %>% str_remove_all(., '-|_|\\..*')) %>%
#   
#   mutate_at('assay_variable', as.character) %>% 
#   mutate_at('biological_replicates', ~str_replace_na(., ''))

# Join WWTP names to qPCR dataset
bb_qpcr <- left_join(polished_results, biobot_lookup, by = 'Biobot ID') %>% mutate_at(c('WWTP', 'FACILITY NAME'), ~if_else(str_detect(., '^X')|is.na(.), assay_variable, .)) %>% 
  mutate_at('Tube ID', ~str_remove(., "\\.")) %>%  unite('Label_tube', c('Sample Name', 'Tube ID'), sep = "", remove = F) %>%  # make a unique column for matching volumes
  rename('Facility' = `FACILITY NAME`) %>% 
  mutate_at('Sample Name', ~as.character(.)) %>% 
  mutate_at('Facility', ~if_else(. == assay_variable, str_c(`Sample Name`, '/', assay_variable), .)) %>% 

  arrange(Facility, biological_replicates) %>% 
  unite('Facility', c(Facility, biological_replicates), sep = "-", na.rm = T) %>% 
  select( -Label_tube, - assay_variable, -`Tube ID`) %>%  
  select(Facility, WWTP, Target, `Copies/ul RNA`, everything())


# Data output ----
# this is usually commented out to prevent overwriting existing data; turn on only when needed for a single run

write_sheet(bb_qpcr,'https://docs.google.com/spreadsheets/d/1I8E8hyF7rFIARTG3ns5khpywafKlkJ6W0qPBxpiUVug/edit#gid=0', sheet = flnm) # save results to a google sheet
# ggsave('qPCR analysis/', WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 8, height = 4)