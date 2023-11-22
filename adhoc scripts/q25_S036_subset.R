# q25_S036_subset.R

# Prelims -----
source('./0-general_functions_main.R') # Source the general_functions file before running this
source('./0.5-user-inputs.R') # source the user inputs from a different script


# User inputs -----

# override title name from user inputs 
title_name <- 'q25_S036'

# load data ----
processed.data <- get_processed_datasets(flnm) %>% 
  
  # S036 data is lysate ; remove U64 target
  filter(Sample_category == 'lysate', Target_name != 'U64') # subset relevant data 
  

# normalizing ----  

normalized.data <- processed.data %>%  
  
  filter(assay_variable != 'ntc') %>% # remove ntcs -- which are NAs
  
  mutate(across(Copies.per.ul.template, # impute the values for missing standard curve Cqs (gfp:gfp)
                ~ if_else(is.na(.), 2^(40-CT), .))) %>% 
  
  pivot_wider(id_cols = c(assay_variable, biological_replicates), 
              names_from = Target_name, values_from = `Copies.per.ul.template`) %>% 
  
  # normalize to 16S
  mutate(signal_per_16S = `gfp:gfp` / `16s`,
         ribozyme_per_16S = gfpbarcode/`16s`,
         .keep = 'unused') %>% 
  
  mutate('signal/ribozyme' = signal_per_16S / ribozyme_per_16S) %>%

  # bring into long format
  pivot_longer(cols = c(signal_per_16S, ribozyme_per_16S), 
               names_pattern = '(.*)(?=_per_16S)',
               
               names_to = 'Target_name', values_to = 'Copies_per_16S') %>% 
  
  mutate(Sample_category = 'Ribogem lysate') # add for plotting

# plotting ----

# plotting Cq
horz.cq <- plot_facetted_assay(.data = processed.data, 
                               .yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                               .xaxis.label.custom = axislabel.assay_variable)

ggsave(str_c('qPCR analysis/', title_name, '.png'), width = 5, height = 3)


# plotting signal normalized to 16S

normalized.plt <- 
  plot_facetted_assay(.data = normalized.data, 
                      .yvar_plot = Copies_per_16S, .xvar_plot = assay_variable)

ggsave(str_c('qPCR analysis/', title_name, '_normalized.png'), width = 5, height = 3)
ggplotly(normalized.plt)
