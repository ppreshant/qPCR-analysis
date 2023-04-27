# ddPCR_analysis.R
# S057j ; S072 ..

# Prelims ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

flnm <- 'S072dd_marine 16S tst_26-4-23' # 'S057j_Marine lysate dils touchup_31Mar23'

# create a titlename : appears on the html file name and header of selected plots, change as required
date_regex <- '[:digit:]*-[:digit:]*-[:digit:]*' # Date regex
title_name <- stringr::str_remove(flnm, stringr::str_c("_", date_regex))

# Load data + meta ----

get_data <- read_csv(str_c('excel files/', flnm, '.csv'))

# Read the sample names and metadata from google sheet
template_source <- 'googlesheet' ; 
plate_template <- get_and_parse_plate_layout(flnm) %>% select(-Target_name) %>% rename(Well = 'Well Position')


# Process data ----

processed_data <- get_data %>% 
  left_join(plate_template) %>%  # join template
  rename(Target_name = Target) %>% # rename target to qPCR format
  # separate(Sample, into = c('Sample_category', 'assay_variable'), sep = '_') # new columns for S057j

  mutate(across(Concentration, as.numeric)) # force conc to be numeric

ratio_data <- processed_data %>% 
  select(Sample_category, assay_variable, Well, Target_name, Concentration) %>% 
  pivot_wider(names_from = Target_name, values_from = Concentration) %>% 
  mutate(copy_number = Plasmid / Chromosome)

# temp save
# ratio_data %>% filter(Sample_category == 'Rd', str_detect(assay_variable, 'gDNA')) %>% write_csv('S057j_tst.csv')

# plot ----

# plot copy number
plot_facetted_assay(.data = ratio_data, 
                    .yvar_plot = copy_number, .xvar_plot = assay_variable, 
                    .facetvar_plot = NULL)

ggsave(plot_as(title_name, '-ddPCR'), width = 4.2, height = 3)

# plot copies per target
plot_facetted_assay(.data = filter(processed_data, !str_detect(assay_variable, '1e3')), 
                    .yvar_plot = CopiesPer20uLWell, .xvar_plot = assay_variable, 
                    .facetvar_plot = Target_name)

ggsave(plot_as(title_name, '-ddPCR_raw'), width = 4.2, height = 4)
