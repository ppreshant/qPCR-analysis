# 16-user_parameters (optional)

# Labelling translators ----


# Parameters for changing y axis variables from short names to long names for labelling on plot

# Subtitle labeller (for y axis variables)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template',
                       'Splicing_fraction_ribozyme' = 'Fraction of the ribozyme mRNA spliced',
                       'Splicing_fraction_16s' = 'Fraction of the 16s rRNA spliced',
                       'expression_ratio' = 'Expression of ribozyme mRNA per 16s rRNA')


# Label x axis (assay_variable) in easy to interpret form 
lst_assay.vars_translation <- list('gfp' = c('89', '315'), # informative_name -> c('assay_variables' ..)
                                   'Ribo' = c('328', '295', '297', '298', '299', '300', '186'),
                                   'Ribo-P1' = '330',
                                   'dead-Ribo' = '54',
                                   'empty' = c('314', '103') ) # use (-)? 

tbl_assay.vars_translation <- lst_assay.vars_translation %>% # convert the list into tibble
  map2_dfr(., names(.), ~ tibble('assay_variable' = .x, 'assay_var.identifier' = .y))
# Add another column in this tibble for specific conversion of 295, etc. into descriptive names?
# make a lookup named vector; use str_replace() to make this new column
