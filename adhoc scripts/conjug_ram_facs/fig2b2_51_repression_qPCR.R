# Read in the final data and make/polish plots for paper

source('./0-general_functions_main.R') # Source the general_functions file before running this

filename <- 'fig2b2_51_repression_qPCR' # filename without the csv suffix
title_name <- 'Ribozyme activity of conjugation plasmid'

qpcr.dat <- read_csv(file = str_c('excel files/paper_data/conjug_ram_facs/', filename, '.csv', sep = ''))

# plot ----
plt_wmean <- {plot_facetted_assay(.data = qpcr.dat, .yvar_plot = Copies.per.ul.template) + 
  geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)} %>% 
  
  format_logscale_y()

# save plot ----
ggsave(filename = str_c('../../Writing/conjug_ram_facs/', title_name, '.png'),
       plt_wmean,
       width = 5, height = 6)



# old processing ----

# proc_data <- get_processed_datasets('q25_S037_RAM repression_14-2-22') # load data from the processed folder
# # Label x axis (assay_variable) : attaches plasmid numbers with informative names for plotting
# assay_var_changer <- c('^51.*' = 'plasmid' , # 'assay_variables' = informative_name
#                        '^67.*' = 'donor')
# 
# 
# p1 <- filter(proc_data, Sample_category != 'lysate',
#              str_detect(assay_variable, '51|^67|ntc')) %>%
# 
#   select(-assay_var.label, -assay_var.horz_label) %>% # remove unnecessary columns
#   select(Target_name, Sample_category, biological_replicates, Copies.per.ul.template, mean_Copies.per.ul.template,
#          everything()) %>%  # set order of columns
#   mutate(assay_var.label = str_replace_all(assay_variable, assay_var_changer), .after = 2)
# 
# 
# write.csv(p1, 'excel files/paper_data/conjug_ram_facs/fig2b2_51_repression_qPCR.csv', na = '')

  