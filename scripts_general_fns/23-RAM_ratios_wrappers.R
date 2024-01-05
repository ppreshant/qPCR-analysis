# 23-RAM_ratios_wrappers.R

#' wrapper function to calculate ratios of U64/16S, gfpbarcode/16S and U64/gfpbarcode per biological replicate 
#' NOTE: all replicates are assumed biological ; need to think how to incorporate technical which I don't do..
#' @param .df data frame holding the absolute copy number data (or est. copy). Defaults to `absolute_dat`
#' @returns pivot data with ratios designated as U64_per_16S, gfpbarcode_per_16S, U64_per_gfpbarcode added within the Target_name column

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
  
    # bring all targets and ratios back into the Target_name column
    pivot_longer(cols = matches('U64|gfpbarcode|16s'), # select all targets and ratios ('_per_' for only ratios)
                 names_to = 'Target_name',
                 values_to = copies_colname) %>% 
    
    # take average of biological replicates
    mutate(across(any_of(copies_colname),
           # c(U64, gfpbarcode, `16s`, # choose each raw column
           #      U64_per_16S, gfpbarcode_per_16S, U64_per_gfpbarcode), # and each ratio
           
           ~ mean(.x, na.rm = TRUE), # and get the mean across biological replicates / remove NAs for robustness
           .names = 'mean_{.col}'), # name columns as mean_..
  
           .by = c(Sample_category, Target_name, matches('assay_var')))  # for each category and sample
  

  
}

#' wrapper function to calculate ratios of Maximal/Repressed per target or signal variable, per biological replicate
#' @param .df data frame holding the absolute copy number data (or est. copy). Defaults to `absolute_dat`
#' @returns pivot data with ratios designated as repression_ratio

calc_repression_ratios <- function(.df = target_ratio)
{
  fold_repression <- 
    .df %>% 
    
    # cleanup samples that are not paired for repressed and maximal
    filter(!str_detect(assay_var.label, 'ntc'), # remove ntc
           !str_detect(Target_name, '^16s')) %>% # remove 16S -- not expected biologically to see variation
    
    # select relevant columns and pivot to make Maximal and Repressed into separate columns
    select(Target_name, Sample_category, assay_var.label, mean_Copies.per.ul.template) %>% 
    unique() %>% # remove the replicates 
    pivot_wider(names_from = Sample_category, values_from = mean_Copies.per.ul.template) %>% 
    
    mutate(repression_fold = Maximal/Repressed, # calculate ratio of the averages of Max / Repressed
           Sample_category = 'all') # add column for plotting
  # Note: doing ratio of averages since the replicates are not paired

}