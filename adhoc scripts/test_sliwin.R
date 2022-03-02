# adhoc inserting sliwin/pcrfit (linregpcr like) workflow into analysis.R

# Load the general functions and the excel file into `fl` using analysis.R

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
  filter(!is.na(Sample_category)) # remove all samples that don't have a Sample_category - implies plate layout was empty
  # This is intended to remove samples whose labels have been removed to prevent analysis

  
# quick visual check
overview_raw_plt <- {ggplot(qpcr_amplification.data,
       aes(x = Cycle,
           y = Rn,
           colour = Target_name,
           group = interaction(Target_name, biological_replicates), # ensures each replicate is a different line
           label = `Well Position`)) +
  geom_line() + 
    # make nested facets : category with children of assay variables
  ggh4x::facet_nested_wrap(. ~ Sample_category + assay_variable)} %>% 
  print()

ggplotly(overview_raw_plt)  # Interactive plot. Use it to figure out the wells containing outliers


# Make plot in logscale    
overview_raw_plt %>% format_logscale_y() %>% print()
# ggsave(plot_as('q25_logy'), width = 5, height = 6)


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
  # 
  # # rename columns to do pcrfit + augment workflow
  # mutate(data = map(data,
  #            ~ rename(., Fluo = Rn, Cycles = Cycle))) %>% 
  
  # do pcrfit to default parameters
  
  # fitting equation (default l4 model): $$f(x) = c + \frac{d - c}{1 + exp(b(log(x) - log(e)))}$$)
  # params: b = hill coeff~ , c = min, d = max, e = Kd like. 
  # Note: I noticed that b is negative, hence the min and max order is c and d..
  
  mutate(qpcr_fit = 
           map(data, 
               ~ pcrfit(.x,
                        cyc = 2,
                        fluo = 3)),
         
         # extract pcrfit parameters to flag no amplifications
         params = map(qpcr_fit,
                      broom::tidy),
         
         # estimate of amplification: max signal/initial
         fold_amplification = map_dbl(params,
                                      ~ .$estimate[3]/.$estimate[2]), # d/c from params
         
         no_amplification_flag = fold_amplification < 3)  # empirically determined for probe qPCR: q25 data                              
 
# aug_amplification.data <- pcrfitted_amplification.data %>% 
#   
#   mutate(
#          # get the fitted values to plot
#          augmented = map(qpcr_fit,
#                          broom::augment) # augment gives error that data has 0 rows.. due to different names? 
#          )



# plot fold_amplification derived from parameters and compare it with amplification curves
pcrfitted_amplification.data %>% 
  ggplot(aes(assay_variable, fold_amplification, colour = Target_name)) +
  geom_point() +
  
  facet_wrap(~ Sample_category)

# Find flags for non amplifying curves..
# sliwin has errors in fitting flat data with no amplification

# trying sliwin
processed_amplification.data <- 
  
  pcrfitted_amplification.data %>% 
  
  # remove samples without amplification
  filter(!no_amplification_flag) %>% 
  
  # find window of linearity (fit sliding window)
  mutate(sliwin_fit = 
           map(qpcr_fit, 
               ~ sliwin(object = .x))
        )
# works now but takes a long time.. 3 mins~

data_to_save <- processed_amplification.data %>% 
  
  # remove unnecessary columns
  select(-data, -qpcr_fit, -params) %>% 
  
  # unpack parameters from sliwin
  mutate(efficiency = map_dbl(sliwin_fit, ~ .$eff),
         initial_conc = map_dbl(sliwin_fit, ~ .$init),
         baseline = map_dbl(sliwin_fit, ~ .$base),
         window_min = map_dbl(sliwin_fit, ~ .$window[1]),
         window_max = map_dbl(sliwin_fit, ~ .$window[2])) %>% 
  
  select(-sliwin_fit)

# save data ---- 
# save data to CSV file

write_csv(data_to_save, file = str_c('excel files/', flnm, '-sliwin.csv'), na = '')
