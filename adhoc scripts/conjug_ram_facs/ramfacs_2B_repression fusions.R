# Load multiple files for comparative analysis

# fig 2B will have only one plasmid set (51) ; 
# Fig S1 (supplementary) : another figure with all fusions (for presentations too)

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q028b_77 79-repressed only_18-5-22',
           'q28_77 79-fusion2_13-4-22',
           'q25_S037_RAM repression_14-2-22')

title_name <- 'Fusion : U64 RAM and fluorescence'

# name translators ----

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
                       'gfpbarcode' = 'Ribozyme',
                       'U64' = 'Spliced')


# Input the data ----

# reading in files and merge : along with the run ID
.df <- get_processed_datasets(flnms)


# Processing ----

forplotting_cq.dat <- filter(.df, !str_detect(Sample_category, 'lysate')) # remove samples that are not relevant

# create a subset of data for presenting the most relevant data : 51, 77, 79 ; 48, 76 with repression and without
presentation.dat <- filter(forplotting_cq.dat, 
                           !str_detect(Sample_category, 'conjug'), # remove samples that are unimportant for presentation
                           !str_detect(assay_variable, 'Sort|sort|328|^67$|48')) %>% # remove sorted samples and other controls
  
  mutate('assay_var.label' = str_replace_all(assay_variable, assay_var_translation)) %>%  # assign labels for readability
  mutate(across(Sample_category, ~ str_replace_all(.x, sample_category_cleaner))) %>%  # clean up redundant names
  
  # arrange for plotting
  mutate(across(Sample_category, ~ fct_relevel(.x, 'Repressed', 'negative', 'Maximal'))) %>%  # set the order of the colours r, b, g
  mutate(across(assay_var.label, ~ fct_reorder(.x, Copies.per.ul.template, .na_rm = FALSE))) %>%  # descending order of copies
  mutate(across(assay_var.label, ~ fct_relevel(.x, 'ntc'))) # Take ntc to the end -- first factor is being plotted closes to origin..

only_u64_nontc <- presentation.dat %>% 
  filter(str_detect(Target_name, 'U64'))


# Plot - presentation ----

plt.copy_ppt <-
  
  {plot_facetted_assay(.data = only_u64_nontc, 
                       .yvar_plot = Copies.per.ul.template, .xvar_plot = assay_var.label, .facetvar_plot = NULL) + 
      
      coord_flip() + 
      # facet_grid(rows = vars(Target_name), scales = 'free_y', space = 'free_y') +
      theme(legend.position = 'top', legend.title= element_blank()) + # stylistic editing of legend 
      xlab('')} %>% 
  format_logscale_y()

# save plot
ggsave(plot_as('all fusions, U64-red'), height = 3, width = 5)
  
# Plotting - all data ----

# plot 40 - Cq
plt.cq_straight <- plot_facetted_assay(.data = presentation.dat, # FOR ALL DATA USE : forplotting_cq.dat
  .yvar_plot = 40-CT) # use ".xvar_plot = assay_variable" for 76 + 67 ; 51 etc on the y axis

# Horizontal orientation : for label readability ; ALREADY DONE WITH ABOVE CODE
# horz.cq <- {plt.cq_straight + 
#     coord_flip() + 
#     facet_grid(rows = vars(Target_name), scales = 'free_y', space = 'free_y') +
#     theme(legend.position = 'top') + 
#     xlab('') + ylab('40 - Cq')} %>% 
#   print()

# save plot ----

ggsave(plot_as('all fusions and conjug, U64-red'), plot = horz.cq, width = 7, height = 5)


# plot copies.per.ul.template for fusions
plt.fusions <- 
  plot_facetted_assay(.data = presentation.dat, 
                      .yvar_plot = Copies.per.ul.template) %>% 
  format_logscale_y() # make logscale on quantity on x-axis (says y since plot is coord_flipped)

ggsave(plot_as('all fusions, U64-red, all targets'), plot = horz.cq, width = 5, height = 6)


# Interactive plot ----

interactive_fusions <- plt.fusions + aes(text = run_ID, label = biological_replicates)
ggplotly(interactive_fusions)


# Calculate fold repression ----

fold_repression <- 
  presentation.dat %>% 
  filter(!str_detect(run_ID, 'q028b'), # remove repeated repressed samples
         !(str_detect(run_ID, 'q28') & str_detect(assay_variable, '51|76')), # remove lone Maximal samples
         !str_detect(assay_var.label, 'ntc'), # remove ntc
         !str_detect(Target_name, '16s')) %>% # remove 16S -- not expected biologically to see variation
  
  # select relevant columns and pivot to make Maximal and Repressed into separate columns
  select(Target_name, Sample_category, assay_var.label, mean_Copies.per.ul.template) %>% 
  unique() %>% # remove the replicates 
  pivot_wider(names_from = Sample_category, values_from = mean_Copies.per.ul.template) %>% 
  
  mutate(repression_fold = Maximal/Repressed, # calculate ratio of the averages of Max / Repressed
         Sample_category = 'all') # add column for plotting
# Note: doing ratio of averages since the replicates are not paired

# show mean data for quick reference
summarise(fold_repression, 
          avg_fold_repression = repression_fold, # take average
          .by = c(Target_name, assay_var.label)) %>%  # for each target and design
  pivot_wider(names_from = Target_name, values_from = avg_fold_repression)


# plot fold change due to repression ----
plt.fold_repression <-
  plot_facetted_assay(.data = fold_repression, 
                    .yvar_plot = repression_fold) + 
  
  aes(label = biological_replicates) # add label for interactive plot

# interactive plot
ggplotly(plt.fold_repression)


# Save data ----

write_csv(presentation.dat, 
          file = str_c('excel files/paper_data/conjug_ram_facs/', 'S1B', '-q52.csv', sep = ''),
          na = '')


# Extra analysis ----

# compare 328 to 76 (since 328 is removed from the all encompassing plot)
.df %>% filter(str_detect(assay_variable, '328|76'), str_detect(Target_name, 'U64')) %>% 
  select(1:3, contains('mean')) %>% unique


## 51 only -----

# plot only 51-maximal for presentation (lab meeting: 3/Apr/24)

sel_51 <- 
  filter(presentation.dat, assay_variable == '51') %>% 
  
  # rename target for ease
  mutate(across(Target_name, ~ str_replace_all(.x, target_translator)))

# print a summary for presentation
summ_sel_51 <- 
  summarise(sel_51, 
            mean_copies = mean(Copies.per.ul.template), 
            
            mean_label = scales::label_scientific(digits = 1)(mean_copies),
            y_position = if_else(mean_copies > 1e4, mean_copies/10, mean_copies * 10),
              
            .by = c(Target_name, assay_variable, Sample_category)) %>% 
  print()

plt.51 <- 
  {plot_facetted_assay(.data = sel_51, 
                      flipped_plot = F, .xvar_plot = assay_variable, 
                      .xaxis.label.custom = 'plasmid', 
                      .yvar_plot = Copies.per.ul.template) + 
      
      # label values
      geom_text(aes(y = y_position, label = mean_label), colour = 'black',
                alpha = 0.5,
                data = summ_sel_51) + 
      
      ggtitle('Reactants vs spliced product') + 
      theme(legend.position = 'none') # remove legend
      
      } %>% 
  format_logscale_y() # make logscale on quantity on x-axis (says y since plot is coord_flipped)

ggsave(plot_as('q25_51-maximal'), plot = plt.51, width = 3, height = 4)


## Ratios to 16S ----

# work in progress!

# Modified code from 2H_all_organisms.R script of RAM
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