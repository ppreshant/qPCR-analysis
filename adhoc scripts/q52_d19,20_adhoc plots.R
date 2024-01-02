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


# Cq plots for D19 only ----

# purpose : to compare head to head with S044_q28 data (pPK077, 79 fusion variants) : is leak significant?
plt_d19_cq <- {plot_facetted_assay(.data = d19, 
                                .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = 40 - CT, 
                                points_plt.style = 'jitter') +
    
    ggtitle('q52_D19_new RAM fusions_repression')} 

ggsave(plot_as('q52_D19_Cq'), plt_d19_cq, width = 5, height = 5)
