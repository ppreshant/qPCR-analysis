# Script to analyze raw data and multicomponent data of qPCR files by hand
# Author : Prashant
# Updated on = 3 Feg 2022

# Initialization ----

# Load a dataset into file by running analysis.R until line 92


# Processing ----

# get the amplification data (Rn) from the qPCR excel file
raw.data <- fl$`Raw Data` %>% 
  left_join(plate_template) %>% # join the samples names from the spreadsheet
  
  # remove un-named samples
  filter(!is.na(Target_name)) %>%  # remove all samples that don't have a target name - implies plate layout was empty
  # This is intended to remove samples whose labels have been removed to prevent analysis
  
  # long format
  pivot_longer(cols = matches('x.-m.'),
               names_to = 'channel',
               values_to = 'signal')

# Plotting ----

# quick visual check
overview_raw_plt <- 
  raw.data %>% 
  
  
  ggplot(aes(x = Cycle,
             y = signal,
             colour = channel,
             group = biological_replicates, # ensures each replicate is a different line
             label = `Well Position`)) +
  geom_line() + 
  facet_grid(Sample_category ~ assay_variable)

ggplotly(overview_raw_plt)  # Interactive plot. Use it to figure out the wells containing outliers


