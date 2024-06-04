# S10_q25.R

# Prerequisites ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

# source('./0.5-user-inputs.R') # source the user inputs from a different script


# Load data ------------------------------

flnm <- 'q25_S037_RAM repression_14-2-22' # overwrites 
title_name <- 'q25_S10c_RAM repression' # make title_name for plot and saving

absolute_dat <- get_processed_datasets(flnm) %>% 
  
  

# name translators ----
# COPIED FROM ramfacs_2b_repression fusions.R

assay_var_translation <- c('76' = 'Ribozyme',    # 'current name' = 'new name' format
                           # '48' = 'mScarlet',  # regex based translation to change x-axis labels
                           
                           '51' = 'Ribozyme-mScarlet',
                           '77' = 'mScarlet-Ribozyme v2',
                           '79' = 'mScarlet-Ribozyme',
                           
                           '67 \\+ ' = '', # remove extra information
                           ' \\+ 67' = '', # the special symbol '+' needs to be double escaped
                           'NTC' = 'ntc')

sample_category_cleaner <- c('test' = 'Repressed', # clean up redundant names
                             'repressed' = 'Repressed',
                             'positive' = 'Maximal',
                             'negative' = 'Negative')


# used only for presentation stuff (extra analysis)
target_translator <- c('16s' = '16S',
                       'gfpbarcode' = 'Ribozyme expression',
                       'U64' = 'Spliced product')

# Processing data ------------------------------

# create a subset of data for presenting : ONLY 76
processed_data <- filter(absolute_dat,
                           str_detect(assay_variable, '76')) %>% # Select only 76
  
  mutate(across(assay_variable, ~ str_replace_all(.x, assay_var_translation))) %>%  # assign labels for readability
  mutate(across(Sample_category, ~ str_replace_all(.x, sample_category_cleaner))) %>%  # clean up redundant names
  
  # rename target for ease
  mutate(across(Target_name, ~ str_replace_all(.x, target_translator)))


label_data <- 
  select(processed_data, Target_name, Sample_category, mean_Copies.per.ul.template) %>%
  unique() # select unique values
  # 
  # # round off mean values
  # mutate(mean_Copies.per.ul.template = round(mean_Copies.per.ul.template, 1))
  

# plot ------------------------------

{ggplot(processed_data, 
        aes(y = Target_name, x = Copies.per.ul.template, 
            shape = Sample_category)) + 
    
    geom_jitter(width = 0.2, alpha = 0.5) + 
    
    # change shapes
    scale_shape_manual('Expression', values = c('Repressed' = 1, 'Maximal' = 19), ) +
    
    
    # Indicate the mean values
    geom_point(aes(x = mean_Copies.per.ul.template), 
               shape = '|', size = 3) + 
    
    # label values
    geom_text(aes(x = mean_Copies.per.ul.template,
                  label = scales::label_scientific(digits = 1)(mean_Copies.per.ul.template),
                  # label = scales::label_parse(mean_Copies.per.ul.template), # for power notation
                  
                  vjust = -2),  # adjust position of label
              
              data = label_data) +

    
    theme_minimal() + theme(legend.position = 'top') +
    
    # adjust labels
    labs(title = NULL, x = 'Copies per ul template', y = NULL)
  
      } %>%  # remove legend
  
  format_logscale_x() # make logscale on quantity on x-axis


# Save plot ------------------------------

ggsave(plot_as(title_name), width = 6, height = 3) # save plot as png  

ggsave(filename = str_c('qPCR analysis/', title_name, '.pdf'), 
       dpi = 300, width = 6, height = 3)
