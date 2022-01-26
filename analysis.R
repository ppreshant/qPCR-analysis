# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

# The last part of the script contains important codes from the console that were used for specific situations : These will be commented

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

flnm <- 'q22_fgfp-triplex-gradient_24-1-22'  
title_name <-'q22_fgfp triplex gradient'

# options
plot_mode <-  'raw_quantification' # Options : ('absolute_quantification' or 'raw_quantification'); 
# absolute_quantification = Will calculate copy #'s based on y_intercept and slope from standard curve - calculated or gathered from old std curves 
# raw_quantification = Cq values are plotted

skip.std.curves_already.exist <- FALSE # If TRUE, will retrieve std curve data from the google sheet
default_std.to.retrieve <-  'Std7' # if the file name doesn't hold any std curve, it will default to this
# This pro-active user setting prevents duplicates being processed into the sheet

# Labelling translators ----

# Subtitle labeller (for y axis variables)
yaxis_translation <- c('40 - CT' = '40 - Cq',
                       'Copies_proportional' = 'Copies proportional (a.u)',
                       'Tm1' = 'Melting temperature - peak 1',
                       'Tm' = 'Melting temperature',
                       'Copies.per.ul.template' = 'Copies per uL template')

# Label x axis (assay_variable) in easy to interpret form 
lst_assay.vars_translation <- list('gfp' = c('89', '315'),
                                'Ribo' = c('328', '295', '297', '298', '299', '300', '186'),
                                'Ribo-P1' = '330',
                                'dead-Ribo' = '54',
                                'empty' = c('314', '103') ) # informative_name -> c('assay_variables' ..)

tbl_assay.vars_translation <- lst_assay.vars_translation %>% # convert the list into tibble
  map2_dfr(., names(.), ~ tibble('assay_variable' = .x, 'assay_var.identifier' = .y))
# Add another column in this tibble for specific conversion of 295, etc. into descriptive names?
# make a lookup named vector; use str_replace() to make this new column

plot_assay_variable <- 'Template name' # printed on the x axis of the graph


# obsolete options ----

experiment_mode <- 'assay' # options ('old_assay' ; 'assay') ; future implementation: 'custom'. Explanation below
# 'assay' =  Plots for Assays (facetted by Target_name, colour by Sample_category = control vs experiment ; 
# naming: primerpairname-overall name_templatename.biologicalreplicatenumber)

# Assay mode features

errorbar_width = 0.1; # width of errorbars - emperically change
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); Determines which variable is chosen for plotting in different colours
std_par <- tibble(                       
  # Input the slope and y_intercept from standard curve of various primer pairs/targets here 
  # Target should match Target field (provided in excel sheet - Sample input reference.csv) 
  target = c('Flipped', 'Unflipped', 'Backbone'),
  slope =  c(-3.36, -3.23, -3.55),
  y_intercept = c(42, 38, 42) # values for primer pairs: Flipped:q4-5. Unflipped:q9-10, Backbone:q12-13
)
plot_select_template <- '' # Options ('' or 'something') ; filters a particular template name to plot 
plot_normalized_backbone <- 'no' # Options: ('yes' or 'no'); plots copy #'s normalized to backbone 
plot_mean_and_sd <- 'yes' # Options: ('yes' or 'no'); plots mean and errorbars instead of each replicate as a point
plot_exclude_category <- '^MHT*' # Regex pattern: 'Controls2', '^MHT*', '^none; exclude Sample_category for plotting; ex: Controls etc.
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; exclude assay_variables for plotting; ex: no template control etc.



# Input the data ----

# reading in file and polishing
fl <- readqpcr(str_c('excel files/',flnm, '.xls')) # read excel file exported by Quantstudio

plate_template <- get_template_for(flnm, sheeturls$plate_layouts_PK) # read samplenames from googlesheets

# this gives a vector to order the samples columnwise in the PCR plate or strip (by default : data is shown row-wise) => This command will enable plotting column wise order
sample_order = columnwise_index(fl) 


Cq_data <- fl$Results %>% 
  select(`Well Position`, CT, starts_with('Tm'), `Target Name`) %>%  # select only the results used for plotting, calculations etc.
  .[sample_order,] %>%  # and arrange them according to sample order
  left_join(plate_template) %>% # join the samples names from the spreadsheet
  
  # useful for multiplexing
  rename('quantstudio_target_name' = 'Target Name') # Target Name from the Quantstudio


# ** Assay mode ----
if(experiment_mode == 'assay')
{
  
  # Parsing names ---- 
  
  polished_cq.dat <- Cq_data %>% 
    
    separate(`Sample_name_bulk`, # Split the components of the sample name bulk by delimiters ('-', '_', '.')
             c('Target_name', 
               'Sample_category',
               'assay_variable',
               'biological_replicates'),
             sep = '-|_|\\.') %>% 
    
    
    mutate(across('assay_variable', as.character)) %>% # useful when plasmid numbers are provided
    mutate(across('biological_replicates', ~str_replace_na(., '')) ) %>% # if no replicates are provided, puts a blank string ('')
    
    filter(!is.na(Target_name)) %>% # remove all samples that don't have a target name - implies plate layout was empty
    # This is intended to remove samples whose labels have been removed to prevent analysis
    
    
    # for TAQMAN (assumed multiplexing)     # replace Target_name from layout with quantstudio
    {if(fl["chemistry_type"] == 'TAQMAN') mutate(., Target_name = quantstudio_target_name) else .} %>% 
    select(-quantstudio_target_name) %>% 
    
    # approx copy #
    mutate('Copies_proportional' = 2 ^ (40-CT)) # copy number is proportional to e^-Cq (e = efficiency ~ close to 2)
   
  # Remove Standards to avoid cluttering plots, and translate assay_variables
  forplotting_cq.dat <- polished_cq.dat %>% 
    filter(Sample_category != 'Std') %>%  # remove standards from the plotted data
    left_join(tbl_assay.vars_translation, by = 'assay_variable') %>%  # attach the assay_var translator
    mutate(assay_var.label = if_else(is.na(assay_var.identifier), # make compound label with translation and original 
                                     assay_variable, # clean up label when no identifier is present 
                                     str_c(assay_var.identifier, assay_variable, sep = '\n')) ) %>% # make compound label
    
    select(Target_name, Sample_category, assay_var.label, CT, everything())
  
  # Cq plot ----
  
  # plot ~ copies (relative quantification)
  plt.copies <- plot_facetted_assay(.yvar_plot = Copies_proportional) 
  # This is a standby plot in the absence of absolute copy number calculation
  
  # plot 40 - Cq
  plt.cq <- plot_facetted_assay(.yvar_plot = 40-CT)
  
  # Tm plots ----
  
  # only plots if it is SYBR_GREEN chemistry ; not TAQMAN
  if(fl["chemistry_type"] != 'TAQMAN') 
  {
    plt.tm1 <- plot_facetted_assay(.yvar_plot = Tm1)
    
    # Gather 4 peaks of melting temperatures into long format
    Tm_data <- forplotting_cq.dat %>% 
      select(Sample_category, assay_variable, biological_replicates, Target_name, assay_var.label,
             starts_with('Tm')) %>%  # select identifiers and Tms
      pivot_longer(cols = starts_with('Tm'), names_to = 'Peak number', values_to = 'Tm') # get all Tms into 1 column
    
    # plot mutliple Tms ; Graph will now show
    plt.alltm <- plot_facetted_assay(.data = Tm_data, .yvar_plot = Tm) + 
      facet_grid(`Peak number` ~ Target_name,  # redo facets
                 scales = 'free_x', space = 'free_x') +
      xlab(plot_assay_variable)
    
    plt.alltm_2 <- Tm_data %>% ggplot(.) + aes(x = assay_variable, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
      theme_classic() + scale_color_brewer(palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
      ggtitle(paste(title_name,': Melting')) + facet_grid(~Sample_category, scales = 'free_x', space = 'free_x') +
      xlab(plot_assay_variable)
  }
  
  # Absolute quantification ----
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    
    # Std curve params ---- 
    
    if(skip.std.curves_already.exist)
    { # if Std curve already exists in sheet, read from google sheet
      
      std_to_retrieve <- str_extract(flnm, 'Std[:alnum:]*') # Find the Std id to retrieve
      std_to_retrieve <- if(is.na(std_to_retrieve)) default_std.to.retrieve else std_to_retrieve # Resort to default if file holds no Std
      
      std_par <- read_sheet(sheeturls$plate_layouts_PK, sheet = 'qPCR Std curves', range = 'A:G', col_types = 'ccnnnnn') %>% 
      filter(str_detect(ID, std_to_retrieve ))
    
      } else {
        
        std_par <- NULL # initialize a dummy std_par
        
        # if it is a standard curve holding file (Stdx), call standard curve processor
        if(str_detect(flnm, 'Std[:digit:]*')) 
        {std_par <- process_standard_curve(flnm, polished_cq.dat) # process the standards within the file
        
        # plot a cq graph with standards included
        plt.cq_w.std <- plot_facetted_assay(.data =  polished_cq.dat, .yvar_plot = 40-CT, .xvar_plot = assay_variable) + # plot 40 - Cq
          theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) # Axis labels vertical
        
        # Plot Tm1 graph with standards included
        plt.tm1_w.std <- plot_facetted_assay(.data = polished_cq.dat, .yvar_plot = Tm1, .xvar_plot = assay_variable) + 
          theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) # Axis labels vertical
        }
        
        
        # if it needs an old standard curve (Stdoldx), bring from google sheet
        if(str_detect(flnm, 'Stdold[:digit:]*')) 
        { # dynamic update of standard curve parameters : default parameters in the inputs_for_analysis.R file
          std_to_retrieve <- str_c('Std', # recreate the Stdxx (xx = number)
                                   str_match(flnm, 'Stdold([:alnum:])*')[2]) # Just take the number
          
          std_par <- read_sheet(sheeturls$plate_layouts_PK, sheet = 'qPCR Std curves', range = 'A:G', col_types = 'ccnnnnn') %>% 
            filter(str_detect(ID, std_to_retrieve ))
        }
    
      }
    
    # If no standard curve is found, then throw an informative error
    if(is.null(std_par)) stop(str_c('No Std curve information found in the file name:',
                                    flnm,
                                    '\n Run with plot_mode <-  raw_quantification if you dont want Std curve information '))
    
    # Use the std curve parameters to back-calculate the absolute copies
    absolute_dat <- forplotting_cq.dat %>%
      select(-Copies_proportional) %>% # remove the dummy copies data
      
      group_by(Target_name) %>%
      nest() %>% # create a new column with data frames for each target
      summarize(w.copy.data = map2(data, Target_name,  # calculate copy number for each dataset, and the mean for replicates
                                    ~ absolute_backcalc(.x, .y, std_par) ) 
      ) %>% 
      unnest(cols = c(w.copy.data)) 
    
      # do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target_name

    
    # plot absolute copies per ul template
    plt.copies <- plot_facetted_assay(.data = absolute_dat, .yvar_plot = Copies.per.ul.template)
    
    plt.copies_w.mean <- plot_facetted_assay(.data = absolute_dat, .yvar_plot = Copies.per.ul.template, points_plt.style = 'jitter') + 
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)
      
  }
  
  # # normalizing copy #s to backbone ----  
  # if(plot_mode == 'absolute_quantification' & plot_normalized_backbone == 'yes')
  # { # computing ratio of copy #s of targets : flipped and unflipped to the backbone
  #   
  #   sel <- absolute_dat %>% select(Sample_category,assay_variable,`Primer pair`,Target,`Copies.per.ul.template`) # select relevant columns (other numeric columns will throw errors)
  #   
  #   sel_b <- sel %>% filter(Target == 'Backbone') # filter out each Target
  #   sel_f <- sel %>% filter(Target == 'Flipped'); sel_u <- sel %>% filter(Target == 'Unflipped');
  #   
  #   sel_f %<>% mutate("Normalized Copies.per.ul.template" = sel_f$`Copies.per.ul.template`/sel_b$`Copies.per.ul.template`); # make ratios to the backbone 
  #   sel_u %<>% mutate("Normalized Copies.per.ul.template" = sel_u$`Copies.per.ul.template`/sel_b$`Copies.per.ul.template`);
  #   
  #   results_ratio <- bind_rows(sel_f, sel_u) # bind results into 1 tibble (for easy plotting)
  #   
  #   # plotting the normalized copy #'s
  #   plt_norm <- results_ratio %>% ggplot(aes(x = `assay_variable`, y = `Normalized Copies.per.ul.template`, color = Target)) +   # plotting
  #     geom_point(size = 2) + facet_grid(~Sample_category, scales = 'free_x', space = 'free_x') # plot points and facetting
  #   
  #   plt_norm.formatted <- plt_norm %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
  #   
  #   print(plt_norm)
  # }
  
}


# ** Old--Assay mode ----

# (facetted by Sample category; naming: 'Sample Name'_variable primer pair)

if (experiment_mode == 'old_assay')
{
  # Separate the sample name into columns and make factors in the right order for plotting (same order as the plate setup)
  
  # isolate the primer pair and assay_variable into 3 columns : Sample name, assay variable and primer pair 
  Cq_data %<>% separate(.,Sample_name,c('Sample Name','Primer pair'),' ') %>% separate(.,Sample_name,c('Sample Name','assay_variable'),'_')
  
  # Factorise the sample name in the order for plotting
  Cq_data %<>% mutate(across(where(is.character),as_factor)) 
  
  # re-arrange the results in same order as the above factors (columnwise order of the plate)
  Cq_data %<>% arrange(`Well Position`) 
  
  # select samples to plot (or to exclude write a similar command)
  Cq_data %<>% filter(str_detect(Sample_name, paste('^', plot_select_template, sep = ''))) # str_detect will find for regular expression; ^x => starting with x
  
  # plot the Tm of multiple peaks in melting curve ; Graph will now show
  
  # Gather the Tm's into another data frame and merge into 1 column
  Tm_data <- Cq_data %>% select(Sample_name, `assay_variable`, `Primer pair`, starts_with('Tm')) %>% gather('Peak number','Tm',-Sample_name, -`Primer pair`, -`assay_variable`)
  
  # plot the Tm ; Graph will now show
  plttm <- Tm_data %>% ggplot(.) + aes(x = `assay_variable`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting')) + facet_grid(~Sample_name, scales = 'free_x', space = 'free_x')
  
  Cq_data %<>% filter(!str_detect(Sample_name, plot_exclude_category)) # exclude unwanted samples categories (sample_name) 
  Cq_data %<>% filter(!str_detect(assay_variable, plot_exclude_assay_variable)) # excluding unwanted samples from assay_variable
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    Cq_data_grouped <- Cq_data %>% group_by(Target) 
    results_abs <- Cq_data_grouped %>% do(., absolute_backcalc(., std_par)) # iteratively calculates copy #'s from standard curve parameters of each Target
    
    if(plot_mean_and_sd == 'yes') {
      y_variable = quo(mean)
      results_abs %<>% group_by(Sample_name, Target, assay_variable) %>% summarise_at(vars(`Copy #`), lst(mean(.,na.rm = T), sd)) 
      # find mean and SD of individual copy #s for each replicate
      } 
    else {y_variable = quo(`Copy #`)}
    
    plt <- results_abs %>% ggplot(aes(x = `assay_variable`, y = !!y_variable, color = !!plot_colour_by)) + ylab('Copy #')    # Specify the plotting variables 

    if(plot_mean_and_sd == 'yes') {plt <- plt + geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd, width = errorbar_width))} # plot errorbars if mean and SD are desired
    
  } 
  
  else plt <- Cq_data %>% ggplot(aes(x = `assay_variable`, y = CT, color = !!plot_colour_by))+ ylab(expression(C[q])) # plot CT values if absolute quantification is not needed
    
  # plot the CT mean and formatting plots
  plt <- plt + geom_point(size = 2) + facet_grid(~Sample_name, scales = 'free_x', space = 'free_x') # plot points and facetting
  plt.formatted <- plt %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
  
  print(plt.formatted)

  # normalizing copy #s to backbone ----  
  if(plot_mode == 'absolute_quantification' & plot_normalized_backbone == 'yes')
  { # computing ratio of copy #s of targets : flipped and unflipped to the backbone
    
    sel <- results_abs %>% select(Sample_name,assay_variable,`Primer pair`,Target,`Copy #`) # select relevant columns (other numeric columns will throw errors)
    
    sel_b <- sel %>% filter(Target == 'Backbone') # filter out each Target
    sel_f <- sel %>% filter(Target == 'Flipped'); sel_u <- sel %>% filter(Target == 'Unflipped');
    
    sel_f %<>% mutate("Normalized copy #" = sel_f$`Copy #`/sel_b$`Copy #`); # make ratios to the backbone 
    sel_u %<>% mutate("Normalized copy #" = sel_u$`Copy #`/sel_b$`Copy #`);
    
    results_ratio <- bind_rows(sel_f, sel_u) # bind results into 1 tibble (for easy plotting)
    
    # plotting the normalized copy #'s
    plt_norm <- results_ratio %>% ggplot(aes(x = `assay_variable`, y = `Normalized copy #`, color = Target)) +   # plotting
    geom_point(size = 2) + facet_grid(~Sample_name, scales = 'free_x', space = 'free_x') # plot points and facetting
    
    plt_norm.formatted <- plt_norm %>% format_classic(., title_name, plot_assay_variable) %>% format_logscale() # formatting plot, axes labels, title and logcale plotting
    
    print(plt_norm)
  }
}

# ggsave('qPCR analysis/S017.png')
# write.xlsx(results_abs, 'excel files/Test.xls', sheetName = 'analysis', append = TRUE, borders = 'surrounding') # saving data table of inferred copy #s to an excel sheet

# Custom plots (Transformed Assay data; Plots copy #; 'Sample Name'_variable 'primer pair') ----

# Save plots manually - copy paste this command to console
# ggsave('qPCR analysis/Chk1.png', dpi = 600)
# ggsave('qPCR analysis/Chk1.png', plot = plttm, dpi = 600)

# assay mode : plotting different colour and facet variables
# plt <- results_abs %>% ggplot(aes(x = `assay_variable`, y = `Copy #`, color = Sample_name)) + ylab('Copy #') +   # plotting
#   scale_y_log10(  # logscale for y axis with tick marks
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x) ))
# 
# plt + geom_point(size = 2, show.legend = T) +
#   theme_classic() + scale_color_brewer(palette="Set1") +
#   theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) +
#   ggtitle(title_name) + xlab(plot_assay_variable) + facet_wrap(~Target, scales = 'free_x')


# Save data output ----

{if(plot_mode == 'absolute_quantification') ( absolute_dat 
)   else forplotting_cq.dat} %>% 
  
  write_csv(.,
            str_c('excel files/processed_data/', flnm, '-processed.csv', sep = ''),
            na = '')

# Plotting into html -----------------------------------------------------------------------


# calling r markdown file
rmarkdown::render('make_html_plots.rmd', output_file = str_c('./qPCR analysis/', title_name, '.html'))
