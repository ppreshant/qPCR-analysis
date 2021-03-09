# writing data and other specific COVID wastewater related functions

# This function writes to the specified google sheet if the current sheet does
# not exist. If the sheet does exist it will ask the user before writing.
check_ok_and_write <- function(data, sheet_URL, title_name)
{
  write_ok <- TRUE
  sheet_dne <- FALSE
  
  # this block sets sheet_dne to true if the sheet does not exist
  sheet_dne <- tryCatch(
    expr = {
      read_sheet(sheet_URL, sheet = title_name)
      message("Sheet already exists")
      FALSE
    },
    error = function(e){
      message('Sheet does not exist')
      return(TRUE)
      print(e)
    }
  )
  
  # if the sheet exists (sheet_dne is false), then ask the user if
  # they want to overwrite. If the user selects cancel, then abort
  if (!sheet_dne) {
    write_ok <- menu(c('Yes', 'No'), title = paste("A sheet with the name", title_name, "already exists. Do you want to overwrite?", sep=" "))
    if (write_ok == 2){
      # stop("Cancel selected, script aborted.")
      print(str_c("Not overwriting: ", title_name))
    }
  }
  
  if (write_ok == 1) {
    write_sheet(data, sheet_URL, sheet=title_name)
  }
}


# Manhole/Bayou Naming Utilities ----
get_manhole_names <- function() {
  lookup_table <- read_sheet(sheeturls$biobot_id, 'All manhole')
  manhole_codes <- lookup_table %>% pull('Facility SYMBOL') %>% str_replace_all(" ", "")
  scrub_spaces <- manhole_codes %>% str_replace_all(" ", "") %>% paste(collapse = "|")
  return(scrub_spaces)
}


get_bayou_names <- function() {
  lookup_table <- read_sheet(sheeturls$biobot_id, 'All Bayou')
  manhole_codes <- lookup_table %>% pull('Facility SYMBOL') %>% str_replace_all(" ", "")
  scrub_spaces <- manhole_codes %>% str_replace_all(" ", "") %>% paste(collapse = "|")
  return(scrub_spaces)
}



# LOD Computation Utilities ----

#' This function runs append_LOD_info() as many time as there are unique targets
#' in the input dataframe under the Target column. Then sticks the data back
#' together in a dataframe. The result should have the same length as the input
#' 
#' fl: a dataframe from qPCR data dump
complete_LOD_table <- function(fl) {
  targs <- unique(fl$Target)
  totalsheet <- fl[0,]
  
  for(i in 1:length(targs)) {
    nextsheet <- append_LOD_info(fl, targs[i])
    totalsheet <- rbind(totalsheet, nextsheet)
  }
  
  return(totalsheet)
  
}



#' This function appends info related to the LOD/LOB/LOQ to the sheet. This
#' function can only do this for one PCR target at a time.
#' 
#' fl: a dataframe from qPCR data dump
#' 
#' targ: the name of a target from the Target column of the dataframe. Example
#' "N1_multiplex", "N2_multiplex", "BCoV" etc.
cappend_LOD_info <- function(fl, targ) {
  fl <- fl %>% filter(Target == targ)
  
  # Pull negative controls out
  negative_controls <- fl %>% filter(Sample_name == 'NTC'| assay_variable == 'DI'| assay_variable == 'BLANK')
  
  # Pull any rows with 3 droplets a.k.a the LOQ
  threes <- fl %>% filter(Positives == 3)
  # If no rows has 3 droplets then the concentration is hard coded to 6
  # Otherwise, take the mean of 3 droplet concentrations
  if(dim(threes)[1] == 0) {
    three_copies <- 0.7
  } else {
    three_copies <- mean(threes$`Copy #`)
  }
  
  # Calculate the LOB
  limit_blank <- mean(negative_controls$`Copy #`) + (1.6 * sd(negative_controls$`Copy #`))
  # Calculate the LOD
  LOD <- three_copies + limit_blank
  
  # Put everything into the table
  new_table <- fl %>% mutate(Positivity = case_when(`Copy #` < LOD ~ "Negative",
                                                    `Copy #` >= LOD ~ "Positive",)) %>%
    mutate(LimitOfDet = LOD)
  
  return(new_table)
}


# Plotting functions ----

# plotting individual replicates
plot_biological_replicates <- function(results_abs, title_text = title_name, xlabel = plot_assay_variable)
{ # Convenient handle for repetitive plotting 'Copy #' vs biological replicate
  
  plt <- results_abs %>% ggplot(aes(x = `Tube ID`, y = `Copy #`, color = Target)) + ylab('Copies/ul RNA extract') +    # Specify the plotting variables 
    geom_point(size = 2) + facet_grid(~`Sample_name`, scales = 'free_x', space = 'free_x') + # plot points and facetting
    ggtitle(title_text) + xlab(xlabel)
  plt.formatted <- plt %>% format_classic(.) %>% format_logscale_y() # formatting plot, axes labels, title and logcale plotting
}
