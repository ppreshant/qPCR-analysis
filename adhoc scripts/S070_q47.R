# S070_q47.R

# Prelim ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
source('./0.5-user-inputs.R') # source the user inputs from a different script / override below
flnm <- 'q47y_S070_Ara 4_14-4-23'
title_name <- 'q47y_S070_Ara'

# select the column name to use for copies -- extrapolated as 2^ (40-Cq) vs absolute quantification
column_for_copies <- quo(Copies.per.ul.template)  # or Copies_proportional for uncalibrated data

# specify targets - as expressions
flipped <- expr(`flipped-v0`)
backbone <- expr(`backbone-v0`)
chromosome <- expr(`chromosome-v0`)

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 
  
  # processing for dose_response (from plate_reader : S070_Ara) / forgot MG1655 in controls..
  mutate(sample_type = if_else(str_detect(assay_variable, 'glu|OFF|ON|MG1655|water'), 'Controls', 'Induction'), # mark controls
       Arabinose = as.numeric(assay_variable), .after = assay_variable) # convert to numeric for plotting properly


# pivot and make ratios ----

grouping_vars_for_ratio <- c('assay_variable', 'Sample_category') # to group after pivoting by target, to take medians etc.

source('scripts_general_fns/22-memory_wrappers_ratio.R')
ratio_data <- calculate_memory_ratios(forplotting_cq.dat)
  

# Plotting ----

source('scripts_general_fns/20-plot_dose_response_and_controls.R') # source the plotter

# plot flipped signal in Cq
plt_flipped <- plot_dose_response_and_controls(.target_to_filter = as_name(flipped)) %>% print

ggsave(plot_as(title_name, '-Cq'), plt_flipped, width = 6, height = 5)

ggplotly(plt_flipped) # only shows the latter plot :(

# plot fraction flipped
plt_fracflip <- plot_dose_response_and_controls(.data = ratio_data, .yvar = flipped_fraction, 
                                                output_all_plots = T)
plt_fracflip[[1]]
ggsave(plot_as(title_name, '-flip_fraction'), width = 6, height = 5)

# save to PDF for paper 
ggsave(str_c('qPCR analysis/', title_name, '.pdf'), width = 5, height = 4)

ggplotly(plt_fracflip[[2]])


# calculations ----

# get summary data : flipped fraction
summary_flipfraction <- 
reframe(ratio_data, #.by = assay_variable,
        median_flipped_fraction) %>% unique

# dynamic range : 1e3 / glu
summary_flipfraction$median_flipped_fraction %>% {.[6] / .[13]}

# dynamic range : ON / OFF
summary_flipfraction$median_flipped_fraction %>% #min(., na.rm = T)
  
  {max(., na.rm = T) / min(., na.rm = T)}


# calculation of flipped copies
flip_summary <- 
  reframe(forplotting_cq.dat, .by = assay_variable,
          median(!!column_for_copies, na.rm = T)) %>% 
  unique %>% print


# adhoc --
# dynamic range ON / glu min
pull(flip_summary, 2) %>% 
  {.[13] / .[4]}

# dynamic range 1e2 (max) / glu min
pull(flip_summary, 2) %>% 
  {max(., na.rm = T) / .[4]}



# extra ----

# chromosome
plt_chromosome <- plot_dose_response_and_controls(.target_to_filter = as_name(chromosome))

ggsave(plot_as(title_name, '-chromosome'), width = 6, height = 5)

plt_chromosome_copies <- plot_dose_response_and_controls(.target_to_filter = as_name(chromosome), .yvar = !!column_for_copies) %>% print
plt_chromosome_copies

# backbone, Cq
plt_bb <- plot_dose_response_and_controls(.target_to_filter = as_name(backbone))
ggsave(plot_as(title_name, '-backbone'), width = 6, height = 5)

# plasmid copy number
plt_copynum <- plot_dose_response_and_controls(.data = ratio_data, .yvar = plasmid_copy_number) %>% print
ggsave(plot_as(title_name, '-copy-num'), width = 6, height = 5)


# Copies : flipped
plt_flipcopy <- plot_dose_response_and_controls(.target_to_filter = as_name(flipped), 
                                                .yvar = !!column_for_copies, layout_widths = c(3,1)) %>% print
ggsave(plot_as(title_name, '-flip_copies'), width = 6, height = 5)

# Copies : backbone
plt_bbcopy <- plot_dose_response_and_controls(.target_to_filter = as_name(backbone), 
                                              .yvar = !!column_for_copies, layout_widths = c(3,1)) %>% print
ggsave(plot_as(title_name, '-bb_copies'), width = 6, height = 5)

# Copies : chromosome
plt_chrcopy <- plot_dose_response_and_controls(.target_to_filter = as_name(chromosome), 
                                              .yvar = !!column_for_copies, layout_widths = c(3,1)) %>% print
ggsave(plot_as(title_name, '-chr_copies'), width = 6, height = 5)

# plot copies of all targets : fig S2
remove_title <- function(pltid) pltid & ggtitle(NULL, subtitle = NULL)

plt_s2 <- 
  (remove_title(plt_flipcopy) + ggtitle('ON state')) / 
  (remove_title(plt_bbcopy) + ggtitle('Total plasmid')) + 
  (remove_title(plt_chrcopy) + ggtitle('Chromosome')) 

# order_assay_variable <- c('OFF', 'ON', 'glu', '0', '1e-1', '1e0', '1e1', '1e2', '1e3', '1e4')
plt_s2

ggsave(str_c('qPCR_analysis/', title_name, '-all copies', '.pdf'), width = 7, height = 7)
