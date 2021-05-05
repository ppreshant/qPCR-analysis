# Ad-hoc script for polishing names and sub-selecting samples for presentation plots
# Prashant
# Date: 5/5/21
# made relevant graph for Joff's army grant (high level overview) : q07 data

# Run analysis.R first to get data in ready format
# forplotting_cq.dat // or absolute_dat

# Tailoring data ----
# self named vector of the subset of names that will be retained in the final name
target.name_translation <- c('U64', '16s', 'barcode') %>% # put things matching multiple values last (like 16s)
  setNames(nm = str_c('.*', ., '.*') ) # name is used for matching as a regular expression

sel_absolute.dat <- absolute_dat %>% 
  filter(assay_variable == '295' & Sample_category == 'Test') %>% # select only U64 and remove NoRT controls
  mutate(across(Target_name, ~ str_replace_all(., target.name_translation)), # translate names to be more descriptive
         across(Target_name, ~ str_replace(., 'U64', 'Spliced')) ) # adhoc naming
  

# View(sel_absolute.dat)


# Plotting ----

tailor.plt.copies_w.mean <- {plot_facetted_assay(.data = sel_absolute.dat, .yvar_plot = Copies.per.ul.template,
                                                .colourvar_plot = NULL) + # make base plot with points
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + # add bars for mean
    
    geom_text(data = ~ select(., -matches('Well|Tm|CT|replicates|^Copies')) %>% unique(), # retain unique data for mean
              mapping = aes(y = mean_Copies.per.ul.template, # label in scientific format
                             label = format(mean_Copies.per.ul.template, digits = 1, scientific = TRUE)),
              parse = TRUE,
              nudge_y = -1) + # move label downwards
    
    # text
    xlab('') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # remove x axis ticks and text
    ggtitle('', subtitle = 'Copies/uL RNA') # remove title and change subtitle
  } %>% 
  format_logscale_y %>% # make y-axis logscale
  print()

ggsave('qPCR analysis/Archive/q07_U64.pdf', tailor.plt.copies_w.mean)
