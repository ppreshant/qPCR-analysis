# Attach metadata and plot linreg pcr analyzed data; Also include calibration processing for absolute quantification
# Prashant K
# 5/March/2022


# User input ----
source('./0.5-user-inputs.R') # source the user inputs from a different script

title_name <- paste0(base_title_name, 
                   '_linreg')


# Run linregpcr ----

if(run_linregpcr)
{
  # source the reticulate interface to linregpcr function
  source('scripts_general_fns/17-run_linregpcr.R') 
  
  # run linregpcr
  run_linregpcr(flnm)
  
}


# Post-processing ----


# pre-requisites ----
source('./0-general_functions_main.R') # Source the general_functions file before running this


# Data input ----

# Read the sample names and metadata from google sheet
plate_template <- get_and_parse_plate_layout(flnm)


# Read the rdmlpython output file (tab separated tsv file)
linreg.results <- read_tsv(str_c('excel files/linregpcr/', flnm, '-linreg', '.csv'),
                           col_types = 'nccccccccnnnnniiinnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnccccccccccccc',
                           na = 'nan') %>%
  
  select(-sample) %>% # remove irrelevant columns
  
  left_join(plate_template, by = c('well' = 'Well Position')) %>%  
  # can rename to Well Position if needed later
  
  
  # clean up targets : Giving primacy to the data file, removing the plate_template targets
  mutate(Target_name = target, .keep = 'unused') %>% 
  
  # remove un-named samples
  filter(!is.na(Sample_category)) # remove all samples that don't have a Sample_category - implies plate layout was empty
# This is intended to remove samples whose labels have been removed to prevent analysis



# Labelling translators ----

# Subtitle labeller (for y axis variables)
# variable : yaxis_translation : Check script 16-user_parameters.R


# x axis labeller : Check script 16-user_parameters.R
# attach explanations to assay_variable (plasmid numbers) for interpretability

axislabel.assay_variable <- 'Template name' # printed on the x axis of the graph


# select useful data ----

# define all the metadata columns coming from the plate_layout
metadata_unique_columns <- c('Sample_category', 'assay_variable', 'Target_name')

# change the data with the windows as line equation y - y2 = m(x - x2)
forplotting_cq.dat <- linreg.results %>% 

  # select only relevant columns
  select(any_of(metadata_unique_columns), # metadata unique
         well, biological_replicates, # addition metadata for replicates
         `n in log phase`, `last log cycle`, `n included`,  # window range
         
         # linregpcr actual outputs
         `mean PCR eff`, `N0 (mean eff)`, `Cq (mean eff)`,
         
         # window anchors
         `baseline`, `log lin cycle`, `log lin fluorescence`) %>% 
  
  # rename columns for ease of reference while plotting and for equations
  rename(n = 'n in log phase',
         x0 = 'log lin cycle', # mid point of log linear phase?
         y0 = 'log lin fluorescence', # Already baseline subtracted value
         xend = 'last log cycle',
         n_sml = 'n included',
         eff = 'mean PCR eff',
         N0 = 'N0 (mean eff)',
         Cq = 'Cq (mean eff)') %>%
  
  # make a col for mean of N0 -- for plotting as a line/boxplot
  group_by(across(all_of(metadata_unique_columns))) %>% 
  mutate(mean_N0 = mean(N0, na.rm = TRUE)) %>% 
  
  # attach the explanatory labels to assay variables
  left_join(tbl_assay.vars_translation, by = 'assay_variable') %>%  # attach the assay_var translator
  mutate(assay_var.label = if_else(is.na(assay_var.identifier), # make compound label with translation and original 
                                   assay_variable, # clean up label when no identifier is present 
                                   str_c(assay_var.identifier, assay_variable, sep = '\n')) ) %>% # make compound label
  
  # change label for horizontal plots
  mutate(assay_var.horz_label = str_replace(assay_var.label, '\n', ' ')) %>% 
  
  select(Target_name, Sample_category, assay_var.horz_label, everything())


# Calibrant data saving ----

if('Std' %in% forplotting_cq.dat$Sample_category) # go into the loop if calibrants present
  
{
  
  source('scripts_general_fns/18-process_calibrants_linregpcr.R') # source the script
  process_calibrants_linregpcr(forplotting_cq.dat) # process the Stds within the file
  
}


# Calibration to standards ----

# retrieve linregpcr calibrant data
mean_calibrant.data <- read_csv('qPCR analysis/Standards/linregpcr_calibrants.csv', col_types = 'ccnnn') %>% 
  select(Target_name, mean_N0.per.copy) # select only relevant columns

# --AddFeature: enable filtering of desired standards by ID when necessary, or ask user when redundancy found


# attach to data (remove any Stds in current data)
calibrated_copies.data <- 
  filter(forplotting_cq.dat, !str_detect(Sample_category, 'Std')) %>% # remove standards from current data
  
  left_join(mean_calibrant.data) %>% # attach the calibration data

  mutate(Copies.per.ul.template = N0/`mean_N0.per.copy`) %>% # normalized to calibration - absolute copies
  mutate(mean_Copies.per.ul.template = mean(Copies.per.ul.template, na.rm = TRUE)) # get mean of replicates (within each group)
  


# Plotting ----

# Plot copies along with mean - after calibration
plt.copies_w.mean <- 
  plot_facetted_assay(.data = calibrated_copies.data,
                      .yvar_plot = Copies.per.ul.template,
                      .colourvar_plot = Sample_category,
                      .facetvar_plot = Target_name,
                      .xvar_plot = assay_variable, .xaxis.label.custom = axislabel.assay_variable,
                      points_plt.style = 'jitter') + 
  
  
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)


# Plot N0 along with mean :: 
plt.copies <- 
  plot_facetted_assay(.yvar_plot = N0,
                      .colourvar_plot = Sample_category,
                      .facetvar_plot = Target_name,
                      .xvar_plot = assay_variable, .xaxis.label.custom = axislabel.assay_variable, 
                      points_plt.style = 'jitter') + 
  
  
  geom_boxplot(aes(y = mean_N0), show.legend = FALSE)

# N0 plot with extra labels on x axis
plt.copies_w.mean_straight <- 
  plot_facetted_assay(.yvar_plot = N0,
                      .colourvar_plot = Sample_category,
                      .facetvar_plot = Target_name,
                      .xvar_plot = assay_var.horz_label, .xaxis.label.custom = axislabel.assay_variable,
                      points_plt.style = 'jitter') + 
  
  
  geom_boxplot(aes(y = mean_N0), show.legend = FALSE)


# saving plot
# ggsave(plot_as('S025_linregpcr_N0 log'), width = 6, height = 4)



# plot 40 - Cq to compare

plt.cq <- plot_facetted_assay(.yvar_plot = 40 - Cq,
                              .colourvar_plot = Sample_category,
                              .facetvar_plot = Target_name,
                              .xvar_plot = assay_variable,
                              .xaxis.label.custom = axislabel.assay_variable)

# Cq plot with extra labels on x axis
plt.cq_straight <- plot_facetted_assay(.yvar_plot = 40 - Cq,
                                       .colourvar_plot = Sample_category,
                                       .facetvar_plot = Target_name,
                                       .xvar_plot = assay_var.horz_label,
                                       .xaxis.label.custom = axislabel.assay_variable)


# ggsave(plot_as('q25_linregpcr cq-2'), width = 6, height = 4)


# Save data ----

write.csv(calibrated_copies.data,
          str_c('excel files/linregpcr_w_metadata/', flnm, '-linreg-processed', '.csv', sep = ''),
          na = '', row.names = FALSE)


# Render plots html ----

# calling r markdown file
rmarkdown::render('make_html_plots.rmd', output_file = str_c('./qPCR analysis/', title_name, '.html'))
