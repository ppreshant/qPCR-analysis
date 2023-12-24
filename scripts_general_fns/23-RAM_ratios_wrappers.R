# 23-RAM_ratios_wrappers.R

#' wrapper function to make ratios of U64/16S, gfpbarcode/16S and U64/gfpbarcode per biological replicate 
#' NOTE: all replicates are assumed biological ; need to think how to incorporate technical which I don't do..
#' @param .df data frame holding the absolute copy number data (or est. copy). Defaults to `absolute_dat`
#' @returns pivot data with ratios designated as U64_per_16S, gfpbarcode_per_16S, U64_per_gfpbarcode

calc_ram_ratios <- function(.df = absolute_dat)
{
  # change the copy number column name if default (Copies.per.ul.template) is not found
  copies_colname <- 'Copies.per.ul.template' %in% colnames(.df) %>% 
    if(.) 'Copies.per.ul.template' else 'Copies_proportional'
  
  ram_ratios <- 
    .df %>% 
    
    # pivot to bring each target into a separate column
    pivot_wider(id_cols = c(Sample_category, biological_replicates, matches('assay_var')),
                names_from = Target_name, values_from = any_of(copies_colname)) %>% 
    
    # calculate ratios
    mutate(U64_per_16S = U64/`16s`,
           gfpbarcode_per_16S = gfpbarcode/`16s`,
           U64_per_gfpbarcode = U64/gfpbarcode
           ) %>% 
    
    # summarize (for future summarizing)
    summarize(across(c(U64, gfpbarcode, `16s`, # choose each raw column
                       U64_per_16S, gfpbarcode_per_16S, U64_per_gfpbarcode), # and each ratio
                     
                     mean), # and get the mean across biological replicates
              .by = c(Sample_category, matches('assay_var'))) # for each category and sample
  
}

#' wrapper function to make ratios of Maximal/Repressed per target or signal variable, per biological replicate
#' @param .df data frame holding the absolute copy number data (or est. copy). Defaults to `absolute_dat`
#' @returns pivot data with ratios designated as repression_ratio

calc_repression_ratios <- function(.df)
{
  fold_repression <- 
    .df %>% 
    
    # cleanup samples that are not paired for repressed and maximal
    filter(!str_detect(assay_var.label, 'ntc'), # remove ntc
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
  
  
}