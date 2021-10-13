# Read in the processed data for specific plotting (papers etc.)
# Author: Prashant Kalvapalle;  Sep 27 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q017a_all organisms_16s_Std17', 
           'q17b_all organisms_duplx U64 gfp_Std8')  # raw filename, without the "-processed" keywork
title_name <-'q17_all organisms'

# options

# x-axis label
plot_assay_variable <- 'Template name' # printed on the x axis of the graph


# Labelling translators ----

# Subtitle labeller (for y axis variables)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template',
                       'Splicing_fraction_ribozyme' = 'Fraction of the ribozyme mRNA spliced',
                       'Splicing_fraction_16s' = 'Fraction of the 16s rRNA spliced',
                       'expression_ratio' = 'Expression of ribozyme mRNA per 16s rRNA')

organism_translation <- c('So' = 'Shewanella', # regex based translation to make full name
                          'Pp' = 'putida',
                          'Ps' = 'stutzeri',
                          'Ec' = 'E coli',
                          'Sm' = 'S meli',
                          'Vn' = 'Vibrio') 

assay_var_translation <- c('328' = 'Ribo', # regex based translation to change x-axis labels
                           '314' = '(-)')

# Input the data ----

# reading in file and polishing
.df <- map_dfr(flnms, 
               ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv')) # read excel file exported by Quantstudio
)

# Mess with the data ----

# any custom processing goes here
forplot_reduced.data <- 
  
  .df %>% # optional filtering step here
  
  separate(assay_var.label,  # separate organism names into a new column
           into = c('organism', 'assay_var.label'),
           sep = ' ',
           fill = 'left') %>%  # for NTC, fill the left one with NA
  
  # translate to fullish name for organism
  mutate(across(organism, ~ str_replace_all(.x, organism_translation)) ) %>%  
  replace_na(replace = list('organism' = '')) %>%  # replace the NA of organism with empty string
  
  # translate to meaningful x axis labels
  mutate(across(assay_var.label, ~ str_replace_all(.x, assay_var_translation)) ) 


# Calculating spliced fraction of the Ribozyme/16s
wider_reduced.dat <- 
  
  forplot_reduced.data %>% 
  mutate(across(Copies_proportional, ~ coalesce(.x, Copies.per.ul.template))) %>%  # take inferred copies into copies_proportional
  
  select(1:4, Copies_proportional, assay_variable, biological_replicates) %>%  # select only required columns
  
  pivot_wider(names_from = Target_name, values_from = Copies_proportional, names_prefix = 'Copies_') %>%  # each target to a col
  
  # calculate ratios of each targets per replicate
  mutate(Splicing_fraction_ribozyme = Copies_U64 / Copies_gfpbarcode,
         Splicing_fraction_16s = Copies_U64 / Copies_16s, 
         expression_ratio = Copies_gfpbarcode / Copies_16s) %>% 
  
  # find the mean across biological replicates for all numeric cols
  group_by(assay_variable) %>% 
  mutate(across(where(is.numeric), 
                mean, ignore.na = TRUE,
                .names = "mean_{.col}")) 
  


# Plotting ----

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
                                          .yvar_plot = Splicing_fraction_ribozyme,
                                          .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                          .facetvar_plot = organism,
                                          points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Splicing_fraction_ribozyme), show.legend = FALSE)  # plot mean
    # ggh4x::facet_nested(cols = vars(Target_name, organism), scales = 'free_x', space = 'free_x')   # redo facets
} %>% 
  
  format_logscale_y() %>% # show logscale
  print()


plt.splicing_ratio_16s <- {plot_facetted_assay(.data = wider_reduced.dat,  # plotting function call with defaults
                                                    .yvar_plot = Splicing_fraction_16s,
                                                    .colourvar_plot = NULL, # remove colour: since only 1 Sample_category
                                                    .facetvar_plot = organism,
                                                    points_plt.style = 'jitter') + 
    
    geom_boxplot(aes(y = mean_Splicing_fraction_16s), show.legend = FALSE)  # plot mean
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



# Save plot

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
