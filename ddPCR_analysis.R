# ddPCR_analysis.R
# S057j ; S072 ..

# Prelims ----

source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script

title_name <- base_title_name


# pre-processing function ----


#' Renames colnames from Quantasoft analysis Pro to match the Quantasoft run data (typical)
#' source: wastewater covid data analysis pipeline/g.10-sheet_columns_renaming_funs.R
#' @param .df dataframe to rename
#' @return dataframe with renamed columns
 
raw_ddpcr_renamer <- function(.df)
{
  .df %>% 
    
    rename(PositiveDroplets = Positives) %>%  # rename to be more meaningful
    
    # rename the column names - if exported from Quantasoft analysis Pro
    rename(CopiesPer20uLWell = matches('Copies/20.*LWell'),
           Concentration = matches('Conc\\(copies/.*L)'),
           AcceptedDroplets = any_of('Accepted Droplets'),
           Threshold = any_of('Threshold1'),
           MeanAmplitudeofPositives = any_of('MeanAmplitudeOfPositives'),  # 'Of' to 'of'
           MeanAmplitudeofNegatives = any_of('MeanAmplitudeOfNegatives')) %>%  
    # rename the column name - if exported from Quantasoft analysis Pro
    
    
    mutate(across(any_of('Concentration'), as.numeric)) %>%  # Remove the NO CALLS and make it numeric column  
    
    # convert ambiguous columns into numeric
    mutate(across(matches('Total|Poisson|Mean|Ch|Ratio|Abundance|Linkage|CNV|Copies|Det'), as.numeric)) %>% 
    mutate(across(where(is.list), as.character)) %>%  # convert any stray lists into character
  
    # remove leading 0 (zero) from A01..A09..H01..H09
    mutate(across('Well', ~ str_replace(., '(?<=[:alpha:])0(?=[:digit:])', '') )) 
  
}

# Load data + meta ----

get_data <- read_csv(str_c('excel files/', flnm, '.csv'))

# Read the sample names and metadata from google sheet
plate_template <- get_and_parse_plate_layout(flnm, read_target_name = FALSE) %>% 
  
  rename(Well = 'Well Position') # rename well position to match the data and to be consise


# Process data ----

processed_data <- get_data %>%
  
  raw_ddpcr_renamer %>% # rename columns (for quantasoft analysis pro exported data)
  
  # join the plate layout/template
  left_join(plate_template) %>%  
  relocate(colnames(plate_template)) %>% # move the columns to the front
  
  rename(Target_name = Target) %>% # rename target to qPCR format
  
  # flag wells to remove from analysis (keeps them in the data output for transparency)
  # remove before doing ratios!
  mutate(flag_wells_to_remove = str_detect(Well, wells_to_remove), .after = Well) %>% # flag wells to remove
  
  # control column types and order
  mutate(across(Concentration, as.numeric)) %>%  # force conc to be numeric
  
  # bring controls to the end
  mutate(across(Sample_category, ~ fct_relevel(.x, 'C-', 'ntc', 'Negative', after = Inf))) %>%
  
  ## custom processing ordering data ----
  mutate(across(assay_variable, ~ fct_relevel(.x, assay_var_order))) # order the assay variables
  

## custom/ratios ----

# custom analysis for Swetha's marine ddPCR data
ratio_data <- processed_data %>% 
  
  # remove flagged wells
  filter(!str_detect(Well, wells_to_remove), 
         !str_detect(Sample_category, 'Negative')) %>% # remove negatives with junk values
  
  # calculate copy number : ratio of plasmid to chromosome
  select(Sample_category, assay_variable, Well, Target_name, Concentration) %>% 
  pivot_wider(names_from = Target_name, values_from = Concentration) %>% 
  mutate(copy_number = Plasmid / Chromosome)


# Save data ----

# save processed data
write_csv(processed_data, 
          str_c('excel files/processed_data/', flnm, '-ddPCR-processed.csv'), 
          na = '')

# Save ratio data
write_csv(ratio_data, 
          str_c('excel files/processed_data/', flnm, '-ddPCR-ratio.csv'), 
          na = '')

# Plotting ----

## copies by target ----

# plot copies per target
copies_all_targets <- 
  plot_facetted_assay(.data = 
                      filter(processed_data, !flag_wells_to_remove), 
                    
                    .yvar_plot = CopiesPer20uLWell, .xvar_plot = assay_variable, 
                    .label.var = Well,
                    .facetvar_plot = Target_name) + 
  
  geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

# save the plot
ggsave(plot_as(title_name, '-ddPCR_raw'), plot = copies_all_targets,
       width = 4.2, height = 4)

# interactive plot
ggplotly(copies_all_targets)


## plasmid copy number (ratio) ----

copy_num <- 
  plot_facetted_assay(.data = ratio_data, 
                      .yvar_plot = copy_number, .xvar_plot = assay_variable, 
                      .label.var = Well,
                      .facetvar_plot = NULL) 

# copy_num + 
#   geom_line(aes(group = Sample_category), alpha = 0.2, show.legend = F) # connect points for easy visual

# save the plot
ggsave(plot_as(title_name, '-ddPCR_plasmid-copies'), plot = copy_num, 
       width = 4.2, height = 3)

# interactive plot
ggplotly(copy_num)


## custom ----

copies_all_targets #+ ylim(c(0, 4300))
ggsave(plot_as(title_name, '-ddPCR_raw'), width = 4.2, height = 4)


