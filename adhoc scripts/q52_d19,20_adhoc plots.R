# q52_d19,20_adhoc plots.R

# Prerequisites ----
source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Load data ----
absolute_dat <- get_processed_datasets(flnm)

# Subset data ----

# D19 data : 1280 series + ntc + MG1655
d19 <- absolute_dat %>% 
  filter(str_detect(assay_variable, '^12..|MG1655|ntc'), # select only the correct ones
         # !str_detect(Target_name, '16s') # remove 16S -- this zooms into more interesting data
         )

plt_d19 <- {plot_facetted_assay(.data = d19, 
                                .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = Copies.per.ul.template, 
                                points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + 
    
    ggtitle('q52_D19_new RAM fusions_repression')} %>% 
  
  format_logscale_y()


ggsave(plot_as('q52_D19_zoom'), plt_d19, width = 5, height = 5)


# D20 data : 15. series + ntc + MG1655
d20 <- absolute_dat %>% 
  filter(str_detect(assay_variable, '^15.|MG1655|ntc'), # select only the correct ones
         !str_detect(Target_name, 'gfpbarcode') # remove 16S
  )

plt_d20 <- {plot_facetted_assay(.data = d20, 
                                .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = Copies.per.ul.template, 
                                points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + 
    
    ggtitle('q52_D19_new RAM fusions_repression')} %>% 
  
  format_logscale_y()

ggsave(plot_as('q52_D20_zoom'), plt_d20, width = 5, height = 5)


# Target_ratio ----
source('scripts_general_fns/23-RAM_ratios_wrappers.R')

# calculate ratios of U64/16S, gfpbarcode/16S and U64/gfpbarcode per biological replicate
target_ratio <- calc_ram_ratios()

## plot target ratio ----

plt.target_ratio <- plot_facetted_assay(.data = target_ratio %>% filter(str_detect(Target_name, '_per_')),
                                        .yvar_plot = Copies.per.ul.template,
                                        .xaxis.label.custom = axislabel.assay_variable) %>% 
  # TODO : add a boxplot for the mean later..
  format_logscale_y()

ggplotly(plt.target_ratio) # interactive plot


# Repression_ratio ----

# calculate ratios of Maximal/Repressed per target or signal variable, per biological replicate
fold_repression <- calc_repression_ratios()

# show only 1280 and 1281 U64 data
fold_repression %>% 
  filter(str_detect(Target_name, '^U64$')) %>% 
  drop_na(repression_fold)

# pivot for quick view of 1280 and 1281 repression : save to lab notebook
fold_repression %>% 
  drop_na(repression_fold) %>% 
  select(Target_name, assay_var.label, repression_fold) %>% 
  # filter(str_detect(Target_name, '^U64$|^gfpbarcode$')) %>% 
  pivot_wider(names_from = Target_name, values_from = repression_fold) %>% 
  view
  



## plot fold change due to repression ----
plt.fold_repression <-
  plot_facetted_assay(.data = fold_repression %>% filter(str_detect(Target_name, 'U64')) %>% drop_na(repression_fold), 
                      .yvar_plot = repression_fold)


# Metaanalysis ---- 

# compare Cqs with q25, q28 data

# Shorten category names
shorten_assay_var <- c('mScarlet-Ribozyme v2' = 'Red-Ribo_v2_77',
                           'mScarlet-Ribozyme' = 'Red-Ribo_79',
                           'Ribozyme-mScarlet' = 'Ribo-Red',
                           'Ribozyme' = 'Ribo')    # 'current name' = 'new name' format

identifier_cols <- c('Target_name', 'Sample_category', 'assay_var.label') # use these columns for each unique replicate

# get q25,28, 28b data
old_fusions <- read_csv(str_c('excel files/paper_data/conjug_ram_facs/', 
                              'S1B-q52.csv')) %>% # read data
  
  # cleanup repeated replicate numbering : for positive control (Max, 51, 76) in q28 (1 to 4)
  mutate_cond((run_ID == 'q28' &
                 str_detect(Sample_category, 'Maximal') &
                 str_detect(assay_var.label, '^Ribo') &
                 biological_replicates == 1),
              biological_replicates = 4) %>% 
  
  # shorten assay variables : readable, less space
  mutate(across(assay_var.label, ~ str_replace_all(.x, shorten_assay_var)))

# combine with Shyam's fusions : q52_d19 data
combined_data <- bind_rows(old_fusions, d19) %>%  # combine old fusions with D19 data
  select(-c(assay_variable, assay_var.horz_label, assay_var.identifier)) %>%  # remove redundant label columns
  
  # clean replicate numbering for Negatives (ntc)
  mutate(biological_replicates = row_number(), .by = all_of(identifier_cols))


# check for clash in biological_replicates numbering
# tst <- select(combined_data, 1,2, assay_var.label, biological_replicates, run_ID)
# group_by(tst, across(c(1:4))) %>% mutate(n = n()) %>% filter(n > 1) %>% view


  
# Cq plots all fusions (S044 and D19) ----

# purpose : to compare head to head with S044_q28 data (pPK077, 79 fusion variants) : is leak significant?
plt_combined <- {plot_facetted_assay(.data = combined_data, 
                                .xvar_plot = assay_var.label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = 40 - CT, 
                                points_plt.style = 'jitter') +
    
    ggtitle('q25,28,52_old and new RAM fusions_repression')} 

ggsave(plot_as('q25,28,52_fusions_Cq'), plt_combined, width = 5, height = 6)

# copies plot
plt_comb_copies <- 
  {plot_facetted_assay(.data = filter(combined_data, Target_name == 'U64'), 
                       .xvar_plot = assay_var.label, 
                       .xaxis.label.custom = axislabel.assay_variable,
                       .yvar_plot = Copies.per.ul.template, 
                       points_plt.style = 'jitter') +
      
      geom_point(aes(y = mean_Copies.per.ul.template), shape = '|', show.legend = FALSE) + # add line for mean
      
      ggtitle('Leak: old and new RAM fusions_q25,28,52')} %>% 
  
  format_logscale_y()

ggsave(plot_as('q25,28,52_fusions_Copies_U64'), plt_comb_copies, width = 6, height = 3)


## processing ----

source('scripts_general_fns/23-RAM_ratios_wrappers.R') # ratios functions

##  ratios calculation ----
comb_target_ratio <- calc_ram_ratios(combined_data) # ratio between targets

comb_fold_repression <- calc_repression_ratios(comb_target_ratio) # ratio between max/repressed

### subset for reporting ----

# subset U64/16S ratio for printing
comb_U64_16s_ratio <- filter(comb_target_ratio, Target_name == 'U64_per_16S') # filter only U64/16S

# Make readable summary of leaky data for Shyam
filter(comb_fold_repression, Target_name == 'U64') %>% # select only U64 data
  
  # mutate(across(where(is.numeric), ~ num(.x, digits = 2))) %>% # round off to 2 decimals
  
  rename(plasmid = assay_var.label) %>% # rename to readable name
  select(-c(Target_name, Sample_category)) %>% 
  view



## plotting ----
plt_combined_U64_16s <- 
  {plot_facetted_assay(.data = comb_U64_16s_ratio, 
                       .xvar_plot = assay_var.label, 
                       .xaxis.label.custom = assay_var.label,
                       .yvar_plot = Copies.per.ul.template, 
                       points_plt.style = 'jitter') +
      
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + # add line for mean
      
      ggtitle('Leak with U64/16S: old and new RAM fusions')} %>% 
  
  format_logscale_y()


plt.fold_repression <-
  plot_facetted_assay(.data = fold_repression %>% filter(str_detect(Target_name, 'U64')) %>% drop_na(repression_fold), 
                      .yvar_plot = repression_fold)