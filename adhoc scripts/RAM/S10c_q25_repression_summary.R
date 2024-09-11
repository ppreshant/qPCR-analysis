# S10_q25.R

# Prerequisites ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

# source('./0.5-user-inputs.R') # source the user inputs from a different script


# Load data ------------------------------

flnm <- 'q25_S037_RAM repression_14-2-22' # overwrites 
title_name <- 'q25_S10c_RAM repression' # make title_name for plot and saving

absolute_dat <- get_processed_datasets(flnm)
  
  

# name translators ----
# COPIED FROM ramfacs_2b_repression fusions.R

sample_category_cleaner <- c('test' = 'Repressed', # clean up redundant names
                             'repressed' = 'Repressed',
                             'positive' = 'Maximal',
                             'negative' = 'Negative')


# used only for presentation stuff (extra analysis)
target_translator <- c('16s' = '16S rRNA',
                       'gfpbarcode' = 'cat-RNA',
                       'U64' = 'barcoded rRNA')

# Processing data ------------------------------

# create a subset of data for presenting : ONLY 76
processed_data <- filter(absolute_dat,
                           str_detect(assay_variable, '76'), # Select only 76
                         !str_detect(Target_name, 'U64')) %>%  # remove U64
  
  mutate(assay_variable = 'Ribozyme') %>%  # assign labels for readability
  mutate(across(Sample_category, ~ str_replace_all(.x, sample_category_cleaner))) %>%  # clean up redundant names
  
  # rename target for ease
  mutate(across(Target_name, ~ str_replace_all(.x, target_translator)))


label_data <- 
  select(processed_data, Target_name, Sample_category, mean_Copies.per.ul.template) %>%
  unique() %>%  # select unique values
  
  # add a column for labels to cat-RNA (3,4)
  mutate(label_categories = c(NA, NA, 'Repressed', 'Constitutive') )
  

# plot ------------------------------

{ggplot(processed_data, 
        aes(y = Target_name, x = Copies.per.ul.template, 
            shape = Sample_category)) + 
    
    geom_jitter(width = 0, height = 0.2, alpha = 0.5) + 
    
    # change shapes
    scale_shape_manual('Expression', values = c('Repressed' = 1, 'Maximal' = 19), ) +
    # remove legend manually?
    
    # Indicate the mean values
    geom_point(aes(x = mean_Copies.per.ul.template), 
               shape = '|', size = 5) + 
    
    # label values
    geom_text(
      aes(x = mean_Copies.per.ul.template,
                  label = scientific_10(mean_Copies.per.ul.template),
      ),
      parse = TRUE,
      # label = scales::label_scientific(digits = 1)(mean_Copies.per.ul.template),
      # label = scales::label_parse(mean_Copies.per.ul.template), # for power notation
                  
      vjust = 2,  # adjust position of label
      
      data = label_data) +
    
    # direct label instead of legend
    geom_text(aes(x = mean_Copies.per.ul.template,
                  label = label_categories,
                  vjust = -2),  # adjust position of label
              
              data = label_data) +
    
    # theme and remove the legend
    theme_classic() + 
    
    theme(legend.position = 'none') +
    
    # adjust labels
    labs(title = NULL, x = 'Copies per ul template', y = NULL)
  
      } %>%  # remove legend
  
  format_logscale_x() # make logscale on quantity on x-axis


# Save plot ------------------------------

ggsave(plot_as(title_name), width = 6, height = 3) # save plot as png  

ggsave(filename = str_c('qPCR analysis/', title_name, '.pdf'), 
       dpi = 300, width = 4, height = 2.5)


# Save data ------------------------------

write.csv(processed_data, str_c('excel files/paper_data/', title_name, '.csv')) # save data as csv

# save clean data for paper
clean_data <- 
  select(processed_data, Target_name, Sample_category, 
         Copies.per.ul.template, mean_Copies.per.ul.template) %>%
  
  # change names to match figure
  mutate(across(Sample_category, ~ str_replace_all(.x, 'Maximal', 'Constitutive')))


output_path <- '../../Writing/RAM paper outputs/Archive'

write.csv(clean_data, str_c(output_path, '/', title_name, '.csv')) # save data as csv


# Statistics ------------------------------

# t-test for repressed vs constitutive
t_test_results <-

filter(processed_data, Target_name == 'cat-RNA') %>%

  t.test(Copies.per.ul.template ~ Sample_category,
         data = .,
         alternative = 'greater') %>% # one-tailed test
  print()

t_test_results$p.value %>% fancy_scientific() # format p-value
