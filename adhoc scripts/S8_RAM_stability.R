# Read in the processed data for specific plotting (papers etc.)
# Author: Prashant Kalvapalle;  Dec 7 2021

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnms <- c('q33_RNA stability-3_20-6-22')  # raw filename, without the "-processed" keyword
title_name <-'S8_RAM stability'

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

target_translation <- c('16s' = '16S rRNA', # regex to change the target names for publication
                        'gfpbarcode' = 'Total barcode',  # 'unspliced CatRNA' used in 2H data..
                        'U64' = 'Barcoded 16S rRNA')

# Input the data ----

# reading in file and polishing
.df <- get_processed_datasets(flnms) # read processed excel file from the analysis.R script


# Processing ----

# any custom processing goes here
forplot_reduced.data <- 
  
  .df %>% # optional filtering step here
  
  separate(assay_var.label,  # separate plasmid names into a new column
           into = c('plasmid', 'time'),
           sep = '/') %>%  # for NTC, fill the right one with NA (default)
  
  # translate to fullish name for plasmid
  mutate(across(plasmid, ~ str_replace_all(.x, plasmid_translation)) ) %>%  
  replace_na(replace = list('plasmid' = '')) %>%  # replace the NA of plasmid with empty string
  
  # translate to full name for target (could change as paper is flushed out)
  mutate(across(Target_name, ~ str_replace_all(.x, target_translation) )) %>% 
  
  
  # convert assay variable (x axis var) to numeric
  mutate(across(time, as.numeric)) %>% 
  
  # code ntcs as 0 (otherwise NAs will not be plotted)
  replace_na(list(time = 0)) #%>% 
  
  # remove extraneous samples
  # filter(Sample_category != 'sludgeconjug') %>% 
  
  
  # (using copies per ul template now)
  # calculate the mean of the Copies_proportional
  # group_by(assay_variable, Target_name) %>% 
  # mutate(across(Copies_proportional,
  #               mean, na.rm = TRUE,
  #               .names = "mean_{.col}")
  # )
  


# Extra processing ----

# subset mean of Copies proportional for visual commenting
visual_summary_mean <- 
  
  ungroup(forplot_reduced.data) %>% # ungroup data 
  
  # select desired cols
  select(plasmid, Target_name, time, mean_Copies.per.ul.template) %>% 
  
  # remove duplicates
  unique() %>%

  # arrange 
  arrange(plasmid, Target_name, time) %>% 
  
  # cleanup decimal points
  mutate(across(mean_Copies.per.ul.template,
                ~ formatC(.x, digits = 1)))
  
  
# exploratory plotting ----

# Cq data
plt.cq <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                               .yvar_plot = 40 - CT, # plot the 40 - Cq values
                               .xvar_plot = time,
                               .colourvar_plot = plasmid, # colour with the plasmid/ntc
                               points_plt.style = 'jitter') +
    
    scale_x_continuous(breaks = c(0, 30, 60, 120, 180))} %>% # adjust the values on x axis
  
  print()


# theoretical copy number data -- without using std curves
plt.copies <- {plot_facetted_assay(.data = forplot_reduced.data,  # plotting function call with defaults
                               .yvar_plot = Copies.per.ul.template, # plot the fake copy # values
                               .xvar_plot = time,
                               .colourvar_plot = plasmid, # colour with the plasmid/ntc
                               points_plt.style = 'jitter') +
    
    geom_line(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + # add a line connecting the mean
    scale_x_continuous(breaks = c(0, 30, 60, 120, 180))} %>% # adjust the values on x axis
  
  format_logscale_y() %>% # format logscale 
  print()




# Normalization ----

# normalize all curves to start from ~ 1 (dividing by the max mean value in each curve)
normalized_RNA <- 
  forplot_reduced.data %>% 
  group_by(Target_name, plasmid) %>% 
  
  nest() %>% # nest all the other variables : for doing normalization and fitting exponential
  
  # bring out the mean and time 0
  mutate(initial_mean_Copies_per_ul = 
           map_dbl(data, 
                   ~ filter(., time == 0) %>% pull(mean_Copies.per.ul.template) %>% unique()),
         
         # Normalization : within each data frame in the nest
         data = map(data, 
                    ~ mutate(., normalized_Copies_per_ul = 
                               .$Copies.per.ul.template / initial_mean_Copies_per_ul) %>%  # divide such that mean starts at 1
                      
                      group_by(time) %>% # group to calculate mean
                      
                      # calculating mean of the normalized
                      mutate(mean_normalized_Copies_per_ul = mean(normalized_Copies_per_ul, na.rm = TRUE))
         ) 
         
  ) # you can do: unnest(normalized_RNA, data) to see the whole dataset



# Exponential fit ----

# fitting exponential curves 
# BUG :: Causing singular gradient error for 16S rRNA ; using purrr::safely() to fix this

safe_exp_fit <- safely(.f = ~ nls(normalized_Copies_per_ul ~ SSasymp(time, ys, y0, log_alpha), data = .x))
# use as map(data, ~ safe_exp(.x)); tidied = map(.fit, ~ broom::tidy(.x$result)) 
# Source : https://aosmith.rbind.io/2020/08/31/handling-errors/

# TODO : need some way to substitute mean for augmented data when no fit


normalized_with_exponential_fit <- 
  filter(normalized_RNA, plasmid == 'Ribo') %>%  # select only the good curves with decreasing trend
  
  mutate(.fit = # making the exponential fit
           map(data, # SSasymp fitting y ~ ys+(y0-ys)*exp(-exp(log_alpha)*time)
               # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/SSasymp
               
               ~ safe_exp_fit(.x)
               # ~ nls(normalized_Copies_per_ul ~ SSasymp(time, ys, y0, log_alpha),
               #       data = .)
           ),
         
         tidied = map(.fit, ~ broom::tidy(.x$result)), # extracting fitting parameters
         augmented = map(.fit, ~ broom::augment(.x$result)) # extrapolating fitting data, for plotting
  ) %>% 
  
  # Get fitting parameters
  
  # unnest the model parameters
  unnest(tidied) %>% 
  
  # arrange the parameter estimate, std. error and other stuff for each paramameter in each column
  pivot_wider(names_from = term,
              values_from = c(estimate, std.error, statistic, p.value)) %>% 
  
  # produce t1/2 estimates
  mutate(t.half = log(2)* exp(-estimate_log_alpha), 
         std.error_t.half = log(2) * exp(-estimate_log_alpha) * std.error_log_alpha,
         
         t.half.text = str_c( format(t.half, digits = 2), 
                              '+/-', 
                              format(std.error_t.half, digits = 2),
                              sep = ' ')
  ) # using error propagation - https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example
  


# Plot normalized ----

# Plot the normalized data with the exponential curves fit
plt.normalized_fits <-
  
  normalized_with_exponential_fit %>% 
  
  # Get fitted data for plotting lines
  unnest(augmented) %>%   # unpack the fitting data for ease of plotting

  
 {plot_facetted_assay(.data = .,  # plotting function call with defaults
                      .xvar_plot = time,
                      .yvar_plot = normalized_Copies_per_ul, # plot the fake copy # values
                      .colourvar_plot = Target_name, # colour with the primer pair 
                      .facetvar_plot = NULL,  # facet by plasmid/ntc
                      points_plt.style = 'jitter',
                      flipped_plot = FALSE) +
    
     geom_line(aes(y = .fitted), linetype = 2, show.legend = FALSE) + # add a line connecting the mean
     scale_x_continuous(breaks = c(0, 30, 60, 120, 180)) +  # adjust the values on x axis
     
     # show t half estimates
     annotate(geom = 'text', x = 120, y = 1, label = 'Half life') + 
     
     geom_text(data = normalized_with_exponential_fit, # Show the calculated half life 
                mapping = aes(x = 120, 
                              y = 1 - (1:nrow(normalized_with_exponential_fit))/9, # space out the values for readability
                              label = t.half.text),
                direction = 'y', force = 10,
                show.legend = FALSE) +

   
     ggtitle(title_name, 
             subtitle = 'Normalized Copies of RNA template')} %>%
  
  # format_logscale_y() %>% # format logscale
  print()


# Add 16S data that is not fitted

plt.normalized_fits_extended <- 
  
{plt.normalized_fits + 
  geom_point(aes(y = normalized_Copies_per_ul),
             data = filter(normalized_RNA, plasmid == 'Ribo', Target_name == '16S rRNA') %>% 
               unnest(data)) + 
    
    theme(legend.position = 'top') + # put legend on top
    # change titles and y axis label for consistency
    ggtitle(NULL, subtitle = NULL) + ylab('Normalized Copies of RNA') + xlab('Time (min)')} %>% 
  
  format_logscale_y() # make logscale


ggsave(str_c('qPCR analysis/Archive/', title_name, '-normalized_fits-extended.pdf'),
       plt.normalized_fits_extended,
       width = 5.2,
       height = 4)

# show dynamic graph
plotly::ggplotly(plt.normalized_fits, dynamicTicks = T)




# plot the normalized data without fitting

plt.normalized_joined <- filter(normalized_RNA, plasmid == 'Ribo') %>% 
  unnest(data) %>% 
  group_by(Target_name, plasmid) %>%
  
  {plot_facetted_assay(.data = .,  # plotting function call with defaults
                       .xvar_plot = time,
                       .yvar_plot = normalized_Copies_per_ul, # plot the fake copy # values
                       .colourvar_plot = Target_name, # colour with the primer pair 
                       .facetvar_plot = NULL ,  # facet by plasmid/ntc
                       points_plt.style = 'jitter',
                       flipped_plot = F) +
      
      geom_line(aes(y = mean_normalized_Copies_per_ul), linetype = 1, show.legend = FALSE) + # add a line connecting the mean
      scale_x_continuous(breaks = c(0, 30, 60, 120, 180)) +  # adjust the values on x axis
      
      # show t half estimates
      annotate(geom = 'text', x = 120, y = 1, label = 'Half life') + 
      
      geom_text(data = normalized_with_exponential_fit, # Show the calculated half life for fitted data
                mapping = aes(x = 120, 
                              y = 1 - (1:nrow(normalized_with_exponential_fit))/9, # space out the values for readability
                              label = t.half.text),
                direction = 'y', force = 10,
                show.legend = FALSE) +
      
      
      ggtitle(title_name, 
              subtitle = 'Normalized Copies of RNA template')} %>%
  
  # format_logscale_y() %>% # format logscale
  print()

ggsave(str_c('qPCR analysis/Archive/', title_name, '-normalized.png'),
       plt.normalized_joined,
       width = 6,
       height = 4)


format_logscale_y(plt.normalized_joined) + theme(legend.position = 'top')

ggsave(str_c('qPCR analysis/Archive/', title_name, '-normalized-logscale.pdf'),
       # plt.normalized_joined,
       width = 5.5,
       height = 4)


# Statistics ---- 

# Timewise tests bw 16S and the other two targets pairwise -
# filter the tibble and row join the 2 subsets with a new index variable
# nest by time and do t.test() ; can unravel with $p.value

timewise_normalized_RNA <- 
  
  # make two separate groups
  {map2(c('^Barcoded', '^Total'), 1:2,
        ~ filter(normalized_RNA, plasmid == 'Ribo', # retain only the comparison data
                 !str_detect(Target_name, .x)) %>% # and the current Targets
          mutate(comparison_group = .y, .before = 1))} %>% 
  
  bind_rows %>% # join the two comparison groups
  unnest(data) %>% # unnest all data
  nest(data = !c(comparison_group, time)) %>% # nest for each comparison
  
  # one-sided t-test 16S > others
  mutate(testing = map(data, ~ t.test(normalized_Copies_per_ul ~ Target_name, alternative = 'greater',
                                      data = .x)),
         pvalue = map_dbl(testing, ~ .x$p.value)) %>% # extract the p values
  
  mutate(`check pvalue < 0.05` = pvalue < 0.05)
  
  


# ------------------------------------------------
  
# Save data ----

# Save normalized data with exponential fit
exponential_fit_data_formatted <- normalized_with_exponential_fit %>% 
  
  # Get fitted data for plotting lines
  unnest(augmented) %>% 
  
  # remove nested and fits data
  select(-data, -.fit) %>% 
  
  # order important data first
  select(Target_name, time, normalized_Copies_per_ul, .fitted, t.half.text, everything())


write_csv(exponential_fit_data_formatted,
          str_c('excel files/paper_data/', title_name, '-normalized.csv', sep = ''),
          na = '')



# Save all data
forplot_reduced.data %>% 
  
  # arrange important columns first
  select(Target_name, plasmid, time, Copies.per.ul.template, mean_Copies.per.ul.template, everything()) %>% 
  
write_csv(str_c('excel files/paper_data/', title_name, '-raw.csv', sep = ''),
          na = '')


# save important data
concise_data <- 
forplot_reduced.data %>% 
  
  # arrange important columns first
  select(Target_name, plasmid, time, Copies.per.ul.template, mean_Copies.per.ul.template) %>% 
  
  arrange(Target_name, plasmid, time) %>% 
  
  write_csv(str_c('excel files/paper_data/', title_name, '-concise.csv', sep = ''),
            na = '')



# summary mean, SD of concise data
concise_summary <- 
  
  concise_data %>% 
  
  # group and calculate stdev
  group_by(Target_name, plasmid, time) %>% 
  mutate(sd_Copies.per.ul.template = sd(Copies.per.ul.template, na.rm = TRUE)) %>% 
  
  select(-Copies.per.ul.template) %>% unique() %>% # remove the individual replicates
  rename_with(~ str_replace(.x, '_Copies.per.ul.template', '')) %>%  # shorten the names of Copies
  
  pivot_wider(names_from = Target_name, values_from = c(mean, sd)) %>% 
  
  write_csv(str_c('excel files/paper_data/', title_name, '-summary.csv', sep = ''),
            na = '')
 

select(timewise_normalized_RNA, -data) %>% 
  
  write_csv(str_c('excel files/paper_data/', title_name, '-stats.csv', sep = ''),
            na = '') # modified by hand, added key for comparison_groups 1, 2 and moved to writing/Archive

# Save exploratory plots ----

# Cq plot
ggsave(str_c('qPCR analysis/Archive/', flnms, '.png'),
       plt.cq,
       width = 6,
       height = 4)
# 
# # with added lines
# ggsave(str_c('qPCR analysis/Archive/', title_name, '-w lines.png'),
#        plt.cq + geom_line(),
#        width = 6,
#        height = 4)

# save Copies proportional plt
ggsave(str_c('qPCR analysis/Archive/', flnms, '-copies.png'),
       plt.copies,
       width = 6,
       height = 4)

