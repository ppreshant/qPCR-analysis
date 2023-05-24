# ddPCR_analysis.R
# S057j ; S072 ..

# Prelims ----

source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script

title_name <- base_title_name

# Load data + meta ----

get_data <- read_csv(str_c('excel files/', flnm, '.csv'))

# Read the sample names and metadata from google sheet
plate_template <- get_and_parse_plate_layout(flnm) %>% select(-Target_name) %>% rename(Well = 'Well Position')


# Process data ----

processed_data <- get_data %>% 
  left_join(plate_template) %>%  # join template
  rename(Target_name = Target) %>% # rename target to qPCR format
  # separate(Sample, into = c('Sample_category', 'assay_variable'), sep = '_') # new columns for S057j
  mutate(across(Concentration, as.numeric)) %>%  # force conc to be numeric
  
  mutate(across(Sample_category, ~ fct_relevel(.x, 'C-', 'ntc', 'control', after = Inf))) # bring controls to the end

ratio_data <- processed_data %>% 
  select(Sample_category, assay_variable, Well, Target_name, Concentration) %>% 
  pivot_wider(names_from = Target_name, values_from = Concentration) %>% 
  mutate(copy_number = Plasmid / Chromosome)

# temp save
# ratio_data %>% filter(Sample_category == 'Rd', str_detect(assay_variable, 'gDNA')) %>% write_csv('S057j_tst.csv')

# plot ----

# plot copy number
copy_num <- 
  plot_facetted_assay(.data = ratio_data, 
                    .yvar_plot = copy_number, .xvar_plot = assay_variable, 
                    .facetvar_plot = NULL) 

copy_num + 
  geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

ggsave(plot_as(title_name, '-ddPCR'), width = 4.2, height = 3)

ggplotly(copy_num)

# plot copies per target
copies_all_targets <- 
plot_facetted_assay(.data = filter(processed_data, !str_detect(assay_variable, '1e3')), 
                    .yvar_plot = CopiesPer20uLWell, .xvar_plot = assay_variable, 
                    .facetvar_plot = Target_name) + 
  geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

ggsave(plot_as(title_name, '-ddPCR_raw'), width = 4.2, height = 4)

ggplotly(copies_all_targets)

# custom plotting

copies_all_targets + ylim(c(0, 4300))
ggsave(plot_as(title_name, '-ddPCR_raw'), width = 4.2, height = 4)
