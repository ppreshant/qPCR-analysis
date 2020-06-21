# Read in the qPCR file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  June 16 2020

# Loading pre-reqisites ----
# 

# Loading libraries, functions and user inputs
source('./general_functions.R') # Source the general_functions file
source('./inputs_for_analysis.R') # Source the file with user inputs

# Data input ----


# Preperation steps
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
plate_template_raw <- read_sheet('https://docs.google.com/spreadsheets/d/19oRiRcRVS23W3HqRKjhMutJKC2lFOpNK8aNUkC-No-s/edit#gid=478762118', sheet = 'Plate import setup', range = 'G1:S9')

# Read in qPCR data and labels from plate template
fl <- readqpcr(flpath) # read excel file exported by Quantstudio
plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') # convert plate template (sample names) into a single vector, columnwise

sample_order = columnwise_index(fl) # this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order

# Load desired qPCR result sheet and columns
bring_results <- fl$Results %>% select(`Well Position`, `Sample Name`, CT, starts_with('Tm'),`Target Name`) %>% rename(Target = `Target Name`) %>%  .[sample_order,] %>%  # select only the results used for plotting, calculations etc. and arrange them according to sample order
 select(-`Sample Name`) %>% right_join(plate_template, by = 'Well Position') %>%  # Incorporate samples names from the google sheet by matching well position
  filter(!is.na(Target))
   
# Remove unneccesary data
rm(fl, plate_template_raw)  # remove old data for sparsity

# Data polishing ----


# Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)

# isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
polished_results <- bring_results %>% separate(`Sample Name`,c(NA, 'Sample Name'),'-') %>% separate(`Sample Name`,c('Sample Name','Tube ID'),'_') %>% 
  mutate(`Tube ID` = if_else(`Sample Name` == 'NTC', '0', `Tube ID`)) %>% 
  separate(`Tube ID`, c('assay_variable', 'biological_replicates'), remove = F) %>%  # Separate out biological replicates 
  unite('Tube ID', c(assay_variable, biological_replicates), sep = '.', remove = F, na.rm = T) %>% # remaking Tube ID - removes spaces after 'dot'
  arrange(assay_variable, biological_replicates) %>% mutate_if(is.character,as_factor) # Factorise the sample name and rearrange in column order of appearance on the plate (for plotting)

# select samples to plot (or to exclude write a similar command)
results_relevant <- polished_results %>% filter(str_detect(`Sample Name`, paste('^', plot_select_facet, sep = ''))) %>%  # Include only desired facets : str_detect will find for regular expression; ^x => starting with x
  filter(!str_detect(`Sample Name`, plot_exclude_facet)) %>%  # exclude unwanted facets (sample_name) 
  filter(!str_detect(assay_variable, plot_exclude_assay_variable)) # excluding unwanted x axis variables from assay_variable

# Computation ----


# Computing copy number from standard curve linear fit information
results_abs <- results_relevant %>% group_by(Target) %>% do(., absolute_backcalc(., std_par)) %>%  # iteratively calculates copy #'s from standard curve parameters of each Target
  mutate(`Copy #` = `Copy #`/template_volume) # Normalizing copy number per micro litre of template in the reaction

# Finding mean and standard deviation within replicates (both technical and biological)

summary_results <- results_abs %>%  group_by(`Sample Name`, Target, assay_variable) %>% summarise_at(vars(`Copy #`), funs(mean(.,na.rm = T), sd)) # find mean and SD of individual copy #s for each replicate
results_abs$`Copy #` %<>% replace_na(0) # make unamplified values 0 for plotting

plt <- results_abs %>% ggplot(aes(x = `Tube ID`, y = `Copy #`, color = Target)) + ylab('Copies/ul RNA extract') +    # Specify the plotting variables 
  geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') + # plot points and facetting
  ggtitle(title_name) + xlab(plot_assay_variable)
plt.formatted <- plt %>% format_classic(.) %>% format_logscale_y() # formatting plot, axes labels, title and logcale plotting

print(plt.formatted)

# Data output ----
# this is usually commented out to prevent overwriting existing data; turn on only when needed for a single run

write_sheet(results_abs,'https://docs.google.com/spreadsheets/d/1ouk-kCJHERRhOMNP07lXfiC3aGB4wtWXpnYf5-b2CI4/edit#gid=0', sheet = title_name) # save results to a google sheet
# ggsave('qPCR analysis/', WW1_Baylor-bovine_pilot.png', plot = plt.formatted, width = 8, height = 4)