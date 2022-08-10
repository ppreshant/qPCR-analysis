# Read in the processed data for specific plotting (papers etc.)
# Author: Prashant Kalvapalle;  Sep 27 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q017a_all organisms_16s_Std17', 
           'q17b_all organisms_duplx U64 gfp_Std8')  # enter the raw filename, without the "-processed" keyword
title_name <-'q17_all organisms'

# options

# x-axis label
plot_assay_variable <- 'Template name' # printed on the x axis of the graph


# Labelling translators ----

# Convert organism to long names
organism_translation <- c('So' = 'S. oneidensis', # regex based translation to make full name
                          'Pp' = 'P. putida',
                          'Ps' = 'P. stutzeri',
                          'Ec' = 'E. coli',
                          'Sm' = 'S. meliloti',
                          'Vn' = 'V. natriegens') 

assay_var_translation <- c('328' = 'CatRNA', # regex based translation to change x-axis labels
                           '314' = '(-)')

target_translation <- c('16s' = '16S rRNA', # regex to change the target names for publication
                        'gfpbarcode' = 'unspliced CatRNA',
                        'U64' = 'barcoded 16S rRNA')

# Input the data ----

# reading in file and polishing
.df <- map_dfr(flnms, 
               ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv'), # read processed csv files
                          col_types = 'cccncnnnncnlcnncc') # specifying column types -- avoids '16s' being read as numbers
)

# Mess with the data ----

# any custom processing goes here
forplot_reduced.data <- 
  
  .df %>% # optional filtering step here
  
  # recalculate mean -- was giving NA for missing replicates
  group_by(Target_name, assay_variable) %>%  # group for calculating mean of replicates
  mutate(mean_Copies.per.ul.template = mean(Copies.per.ul.template, na.rm = TRUE)) %>% # find mean of the replicates
  ungroup() %>% 
  
  separate(assay_var.label,  # separate organism names into a new column
           into = c('organism', 'assay_var.label'),
           sep = ' ',
           fill = 'left') %>%  # for NTC, fill the left one with NA
  
  # translate to fullish name for organism
  mutate(across(organism, ~ str_replace_all(.x, organism_translation)) ) %>%  
  replace_na(replace = list('organism' = '')) %>%  # replace the NA of organism with empty string
  
  # translate to meaningful x axis labels
  mutate(across(assay_var.label, ~ str_replace_all(.x, assay_var_translation)) ) %>% 
  
  # change the target names for publication
  mutate(across(Target_name, ~ str_replace_all(.x, target_translation)))


# Calculating spliced fraction of the Ribozyme/16s
wider_reduced.dat <- 
  
  forplot_reduced.data %>% 
  # mutate(across(Copies.per.ul.template, ~ coalesce(.x, Copies.per.ul.template))) %>%  # take inferred copies into Copies.per.ul.template
  
  select(Target_name, organism, Copies.per.ul.template, assay_variable, assay_var.label, biological_replicates) %>%  # select only required columns
  
  pivot_wider(names_from = Target_name, values_from = Copies.per.ul.template, names_prefix = 'Copies_') %>%  # each target to a col
  
  # calculate ratios of each targets per replicate
  mutate(Spliced_fraction_CatRNA = `Copies_barcoded 16S rRNA` / `Copies_unspliced CatRNA`,
         Spliced_fraction_16s = `Copies_barcoded 16S rRNA` / `Copies_16S rRNA`, 
         expression_ratio = `Copies_unspliced CatRNA` / `Copies_16S rRNA`) %>% 
  
  # find the mean across biological replicates for all numeric cols
  group_by(assay_variable) %>% 
  mutate(across(where(is.numeric), 
                mean, ignore.na = TRUE,
                .names = "mean_{.col}")) 
  

# Save data ----

# save raw data
write_csv(forplot_reduced.data,
          str_c('excel files/paper_data/2H_all organisms', '-raw', '.csv'),
          na = '')

# Replicate data of test/328 data only
individual_copies_data <- 
  
  wider_reduced.dat %>% 
  filter(str_detect(assay_variable, '328')) %>% # only retain the test data with ribozyme
  
  ungroup() %>% # can remove assay_variable column after this
  select(organism, matches('^Copies')) %>%  # select the minimal columns
  
  rename_with(~ str_remove(.x, 'Copies_')) %>%  # remove the Copies in each column name
  mutate(Sample = '+ CatRNA (64)')

write_csv(individual_copies_data,
          str_c('excel files/paper_data/2H_all organisms', '-copies', '.csv'),
          na = '')

# Summary (mean) of replicates
summary_data <- 
  individual_copies_data %>% 
  group_by(organism, Sample) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

write_csv(summary_data,
          str_c('excel files/paper_data/2H_all organisms', '-mean', '.csv'),
          na = '')


# Plotting ----

# plotting all organisms in same panel - with different colours for targets

plt.copies_singlepanel <- 
  forplot_reduced.data %>% 
  filter(str_detect(assay_variable, '328')) %>% # only retain the test data with ribozyme
  
  arrange(mean_Copies.per.ul.template) %>%  # arrange in ascending order
  mutate(across(organism, fct_inorder)) %>% # freeze the order for plotting
  
  {plot_facetted_assay(.data = .,
                      .yvar_plot = Copies.per.ul.template,
                      .xvar_plot = organism,
                      .colourvar_plot = Target_name, 
                      .facetvar_plot = NULL,
                      points_plt.style = 'jitter') 
  
  # show mean 
  # geom_point(aes(y = mean_Copies.per.ul.template), size = 7, shape = '_', show.legend = FALSE)
    } %>% 
  
  format_logscale_y()
  
plt.copies_withcontrols_singlepanel <- 
  forplot_reduced.data %>% 
  
  # rename negative samples to be plotted under organism.. 
  mutate(across(
    organism, 
    ~ if_else(str_detect(assay_variable, '314|NTC'), assay_var.label, .x)
  )) %>% 
  
  arrange(mean_Copies.per.ul.template) %>%  # arrange in ascending order
  mutate(across(organism, fct_inorder)) %>% # freeze the order for plotting
  mutate(across(organism, ~ fct_relevel(.x, c('NTC', '(-)')))) %>% # change order of (-) and NTC
  
  plot_facetted_assay(.data = .,
                       .yvar_plot = Copies.per.ul.template,
                       .xvar_plot = organism,
                       .colourvar_plot = Target_name, 
                       .facetvar_plot = NULL,
                       points_plt.style = 'jitter') %>% 
  
  format_logscale_y()  



# multi-panels for target and organism
plt.copies_w.mean <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                                          .yvar_plot = Copies.per.ul.template,
                                          .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                          points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + # plot mean
    ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


# plotting splicing fractions

plt.splicing_ratio_ribozyme <- {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                                          .yvar_plot = Spliced_fraction_CatRNA,
                                          .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                          .facetvar_plot = organism,
                                          points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Spliced_fraction_CatRNA), show.legend = FALSE)  # plot mean
    # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


plt.splicing_ratio_16s <- {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                                                    .yvar_plot = Spliced_fraction_16s,
                                                    .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                                    .facetvar_plot = organism,
                                                    points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Spliced_fraction_16s), show.legend = FALSE)  # plot mean
  # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()

# plot expression ratio of ribozyme

plt.expression_ratio_ribozyme <- {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                                                    .yvar_plot = expression_ratio,
                                                    .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                                    .facetvar_plot = organism,
                                                    points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_expression_ratio), show.legend = FALSE)  # plot mean
  # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


# Correlation between ribozyme expression and splicing
plt.expression_splicing_correlation <- 
{ggplot(filter(wider_reduced.dat, assay_variable != 'NTC'), 
       aes(x = `Copies_unspliced CatRNA`, y = `Copies_barcoded 16S rRNA`, 
           colour = organism, shape = assay_var.label)) + 
  geom_point() + 
    geom_smooth(method = 'lm', aes(colour = NULL, shape = NULL), show.legend = FALSE) + # draw a correlation line
    scale_shape_manual(values = c(1, 16)) } %>% # points hollow vs filled
  format_logscale_x() %>% 
  format_logscale_y()

# presentation: Correlation between ribozyme expression and splicing (without negatives)
plt.positivesonly_correlation <- 
  {ggplot(filter(wider_reduced.dat, !str_detect(assay_var.label, 'NTC|(-)')), 
          aes(x = `Copies_unspliced CatRNA`, y = `Copies_barcoded 16S rRNA`,
              label = organism)) + 
      geom_point() + 
      
      # ggforce::geom_mark_rect(aes(label = organism)) + labelling points does not work
      geom_smooth(method = 'lm', show.legend = FALSE, alpha = 0.2) # draw a correlation line
      } %>%
  format_logscale_x() %>%
  format_logscale_y()


# Save plot

ggsave(str_c('qPCR analysis/Archive/', title_name, '-all targets.png'),
       plt.copies_singlepanel,
       width = 7,
       height = 4)

ggsave(str_c('qPCR analysis/Archive/', title_name, '-all targets w controls.png'),
       plt.copies_withcontrols_singlepanel,
       width = 8,
       height = 4)


ggsave(str_c('qPCR analysis/Archive/', title_name, '.png'),
       plt.copies_w.mean,
       width = 8,
       height = 4)

# splicing fractions
ggsave(str_c('qPCR analysis/Archive/', 'q17_Splicing ratio ribozyme', '.png'),
       plt.splicing_ratio_ribozyme,
       width = 8,
       height = 4)

ggsave(str_c('qPCR analysis/Archive/', 'q17_Splicing ratio 16s', '.png'),
       plt.splicing_ratio_16s,
       width = 8,
       height = 4)


# expression ratios
ggsave(str_c('qPCR analysis/Archive/', 'q17_expression ratio ribozyme', '.png'),
       plt.expression_ratio_ribozyme,
       width = 8,
       height = 4)

# Correlation
ggsave(plot_as(title_name, '-expression_correlation'), plt.expression_splicing_correlation, 
       width = 5, height = 4)

# correlation - positives only
ggsave(str_c('qPCR analysis/Archive/', title_name, '-correlation_positives.pdf'), plt.positivesonly_correlation,
       width = 4, height = 4)
