# 22-memory_wrappers_ratio.R


#' wrapper function to make ratio from processed data
#' @param .df data frame holding processed data
#' @returns tibble pivot_wider data with ratios of flipped/backbone and backbone/chromosome 

calculate_memory_ratios <- function(.df = forplotting_cq.dat)
{
  
  select(.df, -CT, 
         -any_of(c('mean_Copies.per.ul.template', 'std_curve equation'))) %>% # remove the non unique columns
    pivot_wider(names_from = Target_name, values_from = !!column_for_copies) %>% 
    
    mutate(plasmid_copy_number = !!backbone/!!chromosome,
           flipped_fraction = !!flipped/!!backbone) %>% 
    
    group_by(across(any_of(grouping_vars_for_ratio))) %>% 
    mutate(median_flipped_fraction = median(flipped_fraction, na.rm = T), # median to secure from outliers?
           median_copy_number = median(plasmid_copy_number, na.rm = T))
  
}



#' Set memory target names by version
#' @param target_version target version, suffix to target names : blank = '' or '-v0' 
set_memory_targets <- function(target_version = '')
{
  # how to do this as an expr?
  # future idea
}