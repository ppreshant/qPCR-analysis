# S057j_ddPCR_analysis.R

# Prelims ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Load data ----

flnm <- 'S057j_Marine lysate dils touchup_31Mar23'

get_data <- read_csv(str_c('excel files/', flnm, '.csv'))

# Process data ----

processed_data <- get_data %>% 
  separate(Sample, into = c('organism', 'type'), sep = '_') # new columns


ratio_data <- processed_data %>% 
  select(organism, type, Well, Target, Concentration) %>% 
  pivot_wider(names_from = Target, values_from = Concentration) %>% 
  mutate(copy_number = Plasmid / Chromosome)

# temp save
# ratio_data %>% filter(organism == 'Rd', str_detect(type, 'gDNA')) %>% write_csv('S057j_tst.csv')