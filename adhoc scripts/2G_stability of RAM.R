# Read in the processed data for specific plotting (papers etc.)
# Author: Prashant Kalvapalle;  Dec 7 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q20b_degradation_5-12-21', 
           'q0201a_degradation ntcs_5-12-21')  # raw filename, without the "-processed" keyword
title_name <-'q20_stability RAM'

# options

# x-axis label
plot_assay_variable <- 'Time (min)' # printed on the x axis of the graph


# Labelling translators ----

# Subtitle labeller (for y axis variables)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template',
                       'Splicing_fraction_ribozyme' = 'Fraction of the ribozyme mRNA spliced',
                       'Splicing_fraction_16s' = 'Fraction of the 16s rRNA spliced',
                       'expression_ratio' = 'Expression of ribozyme mRNA per 16s rRNA')

plasmid_translation <- c('328' = 'Ribo', # regex based translation to change x-axis labels
                         '314' = '(-)',
                         '315' = 'gfp')

# Input the data ----

# reading in file and polishing
.df <- map_dfr(flnms, 
               ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv')) # read excel file exported by Quantstudio
)

# retrieve the Qubit RNA concentrations
RNA_concs_raw <- 
  googlesheets4::read_sheet(
  ss = 'https://docs.google.com/spreadsheets/d/1d4k3pzR2nm8PnBeLntDbnlCXC92hIhrUyvwnsqtrpCA/edit#gid=584815397',
  sheet = 'q20: Degradation expt',
  range = 'A3:E23',
  col_types = 'icinn') 

# Mess with the data ----

# any custom processing goes here
forplot_reduced.data <- 
  
  .df %>% # optional filtering step here
  
  separate(assay_var.label,  # separate plasmid names into a new column
           into = c('plasmid', 'assay_var.label'),
           sep = '/') %>%  # for NTC, fill the right one with NA (default)
  
  # translate to fullish name for plasmid
  mutate(across(plasmid, ~ str_replace_all(.x, plasmid_translation)) ) %>%  
  replace_na(replace = list('plasmid' = '')) %>%  # replace the NA of plasmid with empty string
  
  # convert assay variable (x axis var) to numeric
  mutate(across(assay_var.label, as.numeric)) %>% 
  
  # code ntcs as 0 (otherwise NAs will not be plotted)
  replace_na(list(assay_var.label = 0)) %>% 
  
  # remove extraneous samples
  filter(Sample_category != 'sludgeconjug') %>% 
  
  
  # calculate the mean of the Copies_proportional
  group_by(assay_variable, Target_name) %>% 
  mutate(across(Copies_proportional,
                mean, na.rm = TRUE,
                .names = "mean_{.col}"
  )
  )
  

# process the RNA concentration data
RNA_concs <- 
  
  # fill missing values extrapolating from top to bottom
  fill(RNA_concs_raw, Plasmid) %>% 
  
  # change to meaningful plasmid names
  mutate(across(Plasmid, ~ str_replace_all(.x, plasmid_translation)) ) %>%  
  
  # Make a column for replicates
  group_by(Plasmid, `Timepoint (min)`) %>% 
  mutate(bio.replicate = row_number()) %>% 
  
  # make all concentrations in 1 colum : before and after DNAse becomes a category col
  pivot_longer(cols = contains('DNAse'),
               names_to = 'DNAse status',
               values_to = 'RNA concentration') %>% 
  
  # regroup incl before and after
  group_by(Plasmid, `Timepoint (min)`, `DNAse status`) %>% 

  # make a column for mean
  mutate(across(`RNA concentration`,
                mean, na.rm = TRUE,
                .names = "mean_{.col}"
                )
  )
  

# subset mean of Copies proportional for visual commenting : Will change it to the real absolute data later
mean_copies_prop <- 
  
  ungroup(forplot_reduced.data) %>% # ungroup data 
  
  # select desired cols
  select(plasmid, Target_name, assay_var.label, mean_Copies_proportional) %>% 
  
  # remove duplicates
  unique() %>% 

  # arrange 
  arrange(plasmid, Target_name, assay_var.label) %>% 
  
  # cleanup decimal points
  mutate(across(mean_Copies_proportional,
                ~ formatC(.x, digits = 1)))
  
  
# Plotting ----

# Cq data
plt.cq <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                                          .yvar_plot = 40 - CT, # plot the 40 - Cq values
                                          .colourvar_plot = plasmid, # colour with the plasmid/ntc
                                          points_plt.style = 'jitter') +
  
  scale_x_continuous(breaks = c(0, 30, 60, 120, 180))} %>% # adjust the values on x axis
  
  print()


# theoretical copy number data -- without using std curves
plt.copies <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                               .yvar_plot = Copies_proportional, # plot the fake copy # values
                               .colourvar_plot = plasmid, # colour with the plasmid/ntc
                               points_plt.style = 'jitter') +
    
    geom_line(aes(y = mean_Copies_proportional), show.legend = FALSE) + # add a line connecting the mean
    scale_x_continuous(breaks = c(0, 30, 60, 120, 180))} %>% # adjust the values on x axis
  
  format_logscale_y() %>% # format logscale 
  print()



# RNA concentrations
plt.qubit <- 
  {plot_facetted_assay(.data = RNA_concs,  # plotting function call with defaults
                       .xvar_plot = `Timepoint (min)`,
                       .yvar_plot = `RNA concentration`, # plot the 40 - Cq values
                       .colourvar_plot = `DNAse status`, # colour with DNAse before after
                       .facetvar_plot = Plasmid, # facet by the plasmid
                       points_plt.style = 'jitter') +
    
      scale_x_continuous(breaks = c(0, 30, 60, 120, 180)) +
      
      # plot mean by a line
      geom_line(aes(y = `mean_RNA concentration`),
                show.legend = FALSE) + 
      
      ylab('Qubit RNA concentration (ng/ul)')
    
    } %>% # adjust the values on x axis
  
  print()


# Save plots ----

# Cq plot
ggsave(str_c('qPCR analysis/Archive/', title_name, '.png'),
       plt.cq,
       width = 6,
       height = 4)

# with added lines
ggsave(str_c('qPCR analysis/Archive/', title_name, '-w lines.png'),
       plt.cq + geom_line(),
       width = 6,
       height = 4)

# save Qubit plot
ggsave(str_c('qPCR analysis/Archive/', title_name, '-Qubit.png'),
       plt.qubit,
       width = 6,
       height = 4)

# save Copies proportional plt
ggsave(str_c('qPCR analysis/Archive/', title_name, '-copies abstract.png'),
       plt.copies,
       width = 6,
       height = 4)


# Extra analysis ----

# Estimating the unnormalized copy numbers (if RNA concs were not equal before qPCR)
# Run after 2G

RNA_after <- 
  filter(RNA_concs, str_detect(`DNAse status`, 'Before')) %>% 
  rename('assay_var.label' = 'Timepoint (min)',
         'biological_replicates' = bio.replicate,
         'plasmid' = Plasmid) %>% 
  ungroup() %>% 
  select(1, 2, `RNA concentration`, biological_replicates) %>% 
  mutate(diluted_RNA = min(`RNA concentration`))

# estimate copies in undiluted RNA
copies_and_RNA <- 
  left_join(forplot_reduced.data, RNA_after) %>% 
  mutate(undiluted_copies = Copies_proportional / diluted_RNA * `RNA concentration`,
         mean_undiluted_copies = mean_Copies_proportional / diluted_RNA * `RNA concentration`) #%>% 
# filter(!Target_name == '16s')

# normalize all curves to start from 1
normalized_RNA <- 
  copies_and_RNA %>% 
  group_by(Target_name, plasmid) %>% 
  mutate(across(contains('undiluted_copies'), 
                ~ .x/max(.x)))

# theoretical copy number data -- without using std curves
plt.normalized <- {plot_facetted_assay(.data = normalized_RNA,  # plotting function call with defaults
                                   .yvar_plot = undiluted_copies, # plot the fake copy # values
                                   .colourvar_plot = Target_name, # colour with the primer pair 
                                   .facetvar_plot = plasmid,  # facet by plasmid/ntc
                                   points_plt.style = 'jitter') +
    
    geom_line(aes(y = mean_undiluted_copies), show.legend = FALSE) + # add a line connecting the mean
    scale_x_continuous(breaks = c(0, 30, 60, 120, 180)) +  # adjust the values on x axis
    
    ggtitle(title_name, 
            subtitle = 'Normalized to max of 1, estimated undiluted copies')} %>%
  
  # format_logscale_y() %>% # format logscale
  print()

ggsave(str_c('qPCR analysis/Archive/', title_name, '-normalized.png'),
       plt.normalized,
       width = 6,
       height = 4)

# show dynamic graph
plotly::ggplotly(plt.normalized, dynamicTicks = T)