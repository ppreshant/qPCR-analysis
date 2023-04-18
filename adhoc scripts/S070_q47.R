# S070_q47.R

# Prelim ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Load data ----

forplotting_cq.dat <- get_processed_datasets(flnm) %>% 
  
  # processing for dose_response (from plate_reader : S070_Ara) / forgot MG1655 in controls..
  mutate(sample_type = if_else(str_detect(assay_variable, 'glu|OFF|ON|water'), 'Controls', 'Induction'), # mark controls
       Arabinose = as.numeric(assay_variable), .after = assay_variable) # convert to numeric for plotting properly


# pivot and make ratios ----
ratio_data <- select(forplotting_cq.dat, -CT) %>% # remove the non unique columns
  pivot_wider(names_from = Target_name, values_from = Copies_proportional) %>% 
  
  mutate(plasmid_copy_number = backbone/chromosome,
         flipped_fraction = flipped/backbone) %>% 
  
  group_by(assay_variable, Sample_category) %>% 
  mutate(median_flipped_fraction = median(flipped_fraction, na.rm = T), # median to secure from outliers?
         median_copy_number = median(plasmid_copy_number, na.rm = T))

# Plotting ----

source('scripts_general_fns/20-plot_dose_response_and_controls.R') # source the plotter

# plot flipped fraction in Cq
plt_flipped <- plot_dose_response_and_controls() %>% print

ggsave(plot_as(title_name, '-Cq'), plt_flipped, width = 6, height = 5)

ggplotly(plt_flipped) # only shows the latter plot :(

# plot fraction flipped
plt_fracflip <- plot_dose_response_and_controls(.data = ratio_data, .yvar = flipped_fraction) %>% print
ggsave(plot_as(title_name, '-flip_fraction'), width = 6, height = 5)

ggplotly(plt_fracflip)


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
          median(Copies_proportional, na.rm = T)) %>% 
  unique %>% print

# dynamic range ON / glu min
flip_summary$`median(Copies_proportional, na.rm = T)` %>% 
  {.[13] / .[4]}

# dynamic range 1e2 (max) / glu min
flip_summary$`median(Copies_proportional, na.rm = T)` %>% 
  {max(., na.rm = T) / .[4]}



# extra ----

# chromosome
plt_chromosome <- plot_dose_response_and_controls(.target_to_filter = 'chromosome')

ggsave(plot_as(title_name, '-chromosome'), width = 6, height = 5)

plt_chromosome_copies <- plot_dose_response_and_controls(.target_to_filter = 'chromosome', .yvar = Copies_proportional) %>% print
plt_chromosome_copies

# backbone
plt_bb <- plot_dose_response_and_controls(.target_to_filter = 'backbone')
ggsave(plot_as(title_name, '-backbone'), width = 6, height = 5)

# copy number
plt_copynum <- plot_dose_response_and_controls(.data = ratio_data, .yvar = plasmid_copy_number) %>% print
ggsave(plot_as(title_name, '-copies'), width = 6, height = 5)


# Copies_proportional : flipped
plt_flipcopy <- plot_dose_response_and_controls(.target_to_filter = 'flipped', 
                                                .yvar = Copies_proportional) %>% print
ggsave(plot_as(title_name, '-flip_copies'), width = 6, height = 5)

# Copies_proportional : backbone
plt_bbcopy <- plot_dose_response_and_controls(.target_to_filter = 'backbone', 
                                              .yvar = Copies_proportional) %>% print
ggsave(plot_as(title_name, '-bb_copies'), width = 6, height = 5)
