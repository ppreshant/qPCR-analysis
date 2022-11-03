# 16-user_parameters (optional)

# Labelling translators ----


# Parameters for changing y axis variables from short names to long names for labelling on plot

# Subtitle labeller (replaces the y axis names to be readable in the plots)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template',
                       'Splicing_fraction_ribozyme' = 'Fraction of the ribozyme mRNA spliced',
                       'Splicing_fraction_16s' = 'Fraction of the 16s rRNA spliced',
                       'expression_ratio' = 'Expression of ribozyme mRNA per 16s rRNA')


# Label x axis (assay_variable) : attaches plasmid numbers with informative names for plotting
lst_assay.vars_translation <- list('gfp' = c('89', '315'), # informative_name -> c('assay_variables' ..)
                                   'Ribo' = c('328', '295', '297', '298', '299', '300', '186'),
                                   'Ribo-P1' = '330',
                                   'dead-Ribo' = '54',
                                   'empty' = c('314', '103') ) # use (-)? 

tbl_assay.vars_translation <- lst_assay.vars_translation %>% # convert the list into tibble
  map2_dfr(., names(.), ~ tibble('assay_variable' = .x, 'assay_var.identifier' = .y))
# Add another column in this tibble for specific conversion of 295, etc. into descriptive names?
# make a lookup named vector; use str_replace() to make this new column


# obsolete options ----
# might incorporate into analysis.R if needed later - Nov/22 PK

experiment_mode <- 'assay' # options ('old_assay' ; 'assay') ; future implementation: 'custom'. Explanation below
# 'assay' =  Plots for Assays (facetted by Target_name, colour by Sample_category = control vs experiment ; 
# naming: primerpairname-overall name_templatename.biologicalreplicatenumber)
# Feature : What will custom include?

# Assay mode features : Obsolete
plot_normalized_backbone <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to backbone for memory

# exclude Sample_category for plotting; ex: Controls etc.
plot_exclude_category <- '^MHT*' # Regex pattern: 'Controls2', '^MHT*', '^none; 

# exclude assay_variables for plotting; ex: no template control etc.
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; 

