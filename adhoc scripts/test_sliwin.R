# adhoc inserting sliwin/pcrfit (linregpcr like) workflow into analysis.R

# loading libraries ----
library(qpcR) # causing problems with dplyr::select()
select <- dplyr::select # forcing the MASS::select() that masks when loading qpcR


# Processing ----

# Load a dataset into file by running analysis.R until line 92

# get the amplification data (Rn) from the qPCR excel file
qpcr_amplification.data <- fl$`Amplification Data` %>% 
  left_join(plate_template) %>% # join the samples names from the spreadsheet
  
  # Target name readjustment for probe assays (TAQMAN), for multiplexing
  
  # rename the target name from the quantstudio file
  rename('quantstudio_target_name' = 'Target Name') %>% # Target Name from the Quantstudio

  # for TAQMAN (assumed multiplexing)     # replace Target_name from layout with quantstudio
  {if(fl["chemistry_type"] == 'TAQMAN') mutate(., Target_name = quantstudio_target_name) else .} %>% 
  select(-quantstudio_target_name) %>% 
  
  
  # remove un-named samples
  filter(!is.na(Target_name)) # remove all samples that don't have a target name - implies plate layout was empty
  # This is intended to remove samples whose labels have been removed to prevent analysis

  
# quick visual check
overview_raw_plt <- ggplot(qpcr_amplification.data,
       aes(x = Cycle,
           y = Rn,
           colour = Target_name,
           group = interaction(Target_name, biological_replicates), # ensures each replicate is a different line
           label = `Well Position`)) +
  geom_line() + 
  facet_grid(Sample_category ~ assay_variable)

ggplotly(overview_raw_plt)  # Interactive plot. Use it to figure out the wells containing outliers


# qpcR fitting ----

# nested data for model fitting
nested_amplification.data <- qpcr_amplification.data %>% 
  
  # each unique sample and replicate in a separate group
  group_by(Sample_category, assay_variable, biological_replicates, Target_name, `Well Position`) %>% 
  nest() %>% 
  
  
  # flag samples with negative Rn (problems with pcrfit)
  mutate(negative_Rn_flag = map_lgl(data, ~ any(.$Rn < 0)))

# Inform of flagged samples and ask to proceed..
# print('Samples in Wells : ', )
  
# pcrfit to data
pcrfitted_amplification.data <- 
  
  nested_amplification.data %>% 
  filter(!negative_Rn_flag) %>% # remove flagged samples with negative Rn
  
  # do pcrfit to default parameters
  mutate(qpcr_fit = 
           map(data, 
               ~ pcrfit(.x,
                        cyc = 2,
                        fluo = 3)))

# Find flags for non amplifying curves..
# sliwin has errors in fitting flat data with no amplification

# trying sliwin
processed_amplification.data <- 
  
  pcrfitted_amplification.data %>% 
  
  # find window of linearity (fit sliding window)
  mutate(sliwin_fit = 
           map(qpcr_fit, 
               ~ sliwin(object = .x))
        )
         
  
  




qpcr_sliwin <- sliwin(qpcr_fit)
# .Error in GRID[i, 1]:GRID[i, 4] : NA/NaN argument -- Is this expecting 4 channels of fluorescence?