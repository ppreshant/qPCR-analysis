# 19-calibrate_processed_data.R
#' Load processed data, retrieve default std curve parameters and calculate absolute copies

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name


# Load data ----
cq.dat <- get_processed_datasets(flnm)


# Get std curve ----
std_to_retrieve <- default_std.to.retrieve

std_par <- 
  {if(template_source == 'googlesheet')
    googlesheets4::read_sheet(sheeturls$plate_layouts_PK, sheet = 'qPCR Std curves', range = 'A:G', col_types = 'ccnnnnn') else {
      read_csv(file = 'qPCR analysis/Standards/qPCR_Std_curve_parameters.csv')}
  } %>% 
  
  filter(str_detect(ID, std_to_retrieve ))


# Calibration ----
absolute_dat <- calculate_absolute_quant(cq.dat, std_par)


# Plot ----

plt.copies_w.mean <- plot_facetted_assay(.data = absolute_dat, 
                                         .yvar_plot = Copies.per.ul.template, 
                                         .xaxis.label.custom = axislabel.assay_variable, 
                                         points_plt.style = 'jitter', flipped_plot = F) + 
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)

plt.copies_w.mean

# Save data ----

# overwrites existing data
write_csv(absolute_dat,
          str_c('excel files/processed_data/', flnm, '-processed.csv', sep = ''),
          na = '')

