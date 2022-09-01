# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

source('./0-general_functions_main.R') # Source the general_functions file before running this

# User inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Labelling translators ----

# Subtitle labeller (for y axis variables)
# variable : yaxis_translation : Check script 16-user_parameters.R

# x axis labeller
# attach explanations to assay_variable (plasmid numbers) for interpretability


axislabel.assay_variable <- 'Template name' # printed on the x axis of the graph


# obsolete options ----

experiment_mode <- 'assay' # options ('old_assay' ; 'assay') ; future implementation: 'custom'. Explanation below
# 'assay' =  Plots for Assays (facetted by Target_name, colour by Sample_category = control vs experiment ; 
# naming: primerpairname-overall name_templatename.biologicalreplicatenumber)
# What will custom include?

# Assay mode features

errorbar_width = 0.1; # width of errorbars - emperically change
# Determines which variable is chosen for plotting in different colours
plot_colour_by <- quo(Target) # Options : (quo(Target) or quo(Sample Name); 

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
# exclude Sample_category for plotting; ex: Controls etc.
plot_exclude_category <- '^MHT*' # Regex pattern: 'Controls2', '^MHT*', '^none; 
# exclude assay_variables for plotting; ex: no template control etc.
plot_exclude_assay_variable <- '^none' # Regex pattern: '^N', '^none' or ''; 




# Input the data ----

# reading in qPCR data excel file
fl <- readqpcr(str_c('excel files/',flnm, '.xls')) # read excel file exported by Quantstudio

# Read the sample names and metadata from google sheet
plate_template <- get_and_parse_plate_layout(flnm)

# this gives a vector to order the samples columnwise in the PCR plate or strip 
# (by default : data is shown row-wise) => This command will enable plotting column wise order
sample_order = columnwise_index(fl) 


Cq_data <- fl$Results %>% 
  
  # select only the results used for plotting, calculations etc.
  select(`Well Position`, CT, starts_with('Tm'), `Target Name`) %>%  
  .[sample_order,] %>%  # and arrange them according to sample order
  left_join(plate_template) %>% # join the samples names from the spreadsheet
  
  # Target name readjustment for probe assays (TAQMAN), for multiplexing
  # rename the target name from the quantstudio file
  rename('quantstudio_target_name' = 'Target Name') %>% # Target Name from the Quantstudio
  
  # for TAQMAN (assumed multiplexing)     # replace Target_name from layout with quantstudio
  {if(fl["chemistry_type"] == 'TAQMAN') mutate(., Target_name = quantstudio_target_name) else .} %>% 
  select(-quantstudio_target_name) 



# ** Assay mode ----
if(experiment_mode == 'assay')
{
  
  # Polishing data ---- 
  
  polished_cq.dat <- Cq_data %>% 
  
    # remove un-named samples
    filter(!is.na(Target_name)) %>% # remove all samples that don't have a target name - implies plate layout was empty
    # This is intended to remove samples whose labels have been removed to prevent analysis
    
    # approx copy #
    mutate('Copies_proportional' = 2 ^ (40-CT)) # copy number is proportional to e^-Cq (e = efficiency ~ close to 2)
   
  # Remove Standards to avoid cluttering plots, and translate assay_variables
  forplotting_cq.dat <- polished_cq.dat %>% 
    filter(Sample_category != 'Std') %>%  # remove standards from the plotted data
    left_join(tbl_assay.vars_translation, by = 'assay_variable') %>%  # attach the assay_var translator
    mutate(assay_var.label = if_else(is.na(assay_var.identifier), # make compound label with translation and original 
                                     assay_variable, # clean up label when no identifier is present 
                                     str_c(assay_var.identifier, assay_variable, sep = '\n')) ) %>% # make compound label
    
    # change label for horizontal plots
    mutate(assay_var.horz_label = str_replace(assay_var.label, '\n', ' ')) %>% 
    
    select(Target_name, Sample_category, assay_var.horz_label, CT, everything())
  
  # Cq plot ----
  
  # plot ~ copies (relative quantification)
  plt.copies <- plot_facetted_assay(.yvar_plot = Copies_proportional, .xaxis.label.custom = axislabel.assay_variable) 
  # This is a standby plot in the absence of absolute copy number calculation
  
  # plot 40 - Cq
  plt.cq <- plot_facetted_assay(.yvar_plot = 40-CT, .xaxis.label.custom = axislabel.assay_variable)
  plt.cq_straight <- plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                                         .xaxis.label.custom = axislabel.assay_variable)
  
  # Tm plots ----
  
  # only plots if it is SYBR_GREEN chemistry ; not TAQMAN
  if(fl["chemistry_type"] != 'TAQMAN') 
  {
    plt.tm1 <- plot_facetted_assay(.yvar_plot = Tm1, .xaxis.label.custom = axislabel.assay_variable)
    
    # Gather 4 peaks of melting temperatures into long format
    Tm_data <- forplotting_cq.dat %>% 
      select(Sample_category, assay_variable, biological_replicates, Target_name, assay_var.label,
             starts_with('Tm')) %>%  # select identifiers and Tms
      pivot_longer(cols = starts_with('Tm'), names_to = 'Peak number', values_to = 'Tm') # get all Tms into 1 column
    
    # plot mutliple Tms ; Graph will now show
    plt.alltm <- plot_facetted_assay(.data = Tm_data, .yvar_plot = Tm, .xaxis.label.custom = axislabel.assay_variable) + 
      facet_grid(`Peak number` ~ Target_name,  # redo facets
                 scales = 'free_x', space = 'free_x') +
      xlab(axislabel.assay_variable)
    
    plt.alltm_2 <- Tm_data %>% ggplot(.) + aes(x = assay_variable, y = Tm) + 
      geom_point(aes(color = `Peak number`), size = 2) +
      theme_classic() + scale_color_brewer(palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
      ggtitle(paste(title_name,': Melting')) + facet_grid(~Sample_category, scales = 'free_x', space = 'free_x') +
      xlab(axislabel.assay_variable)
  }
  
  # Absolute quantification ----
  
  if(plot_mode == 'absolute_quantification')
  { # Computing copy number from standard curve linear fit information
    
    # Std curve params ---- 
    
    if(skip.std.curves_already.exist)
    { # if Std curve already exists in sheet, read from google sheet
      
      std_to_retrieve <- str_extract(flnm, 'Std[:alnum:]*') # Find the Std id to retrieve
      std_to_retrieve <- 
        if(force.use_default.std || # if forced to use the default or if Std isn't in the file
           is.na(std_to_retrieve)) default_std.to.retrieve else std_to_retrieve # Resort to default if file holds no Std
      
      std_par <- read_sheet(sheeturls$plate_layouts_PK, sheet = 'qPCR Std curves', range = 'A:G', col_types = 'ccnnnnn') %>% 
      filter(str_detect(ID, std_to_retrieve ))
    
      } else {
        
        std_par <- NULL # initialize a dummy std_par
        
        # if it is a standard curve holding file (Stdx), call standard curve processor
        if(str_detect(flnm, 'Std[:digit:]*')) 
        {std_par <- process_standard_curve(flnm, polished_cq.dat, dilutions_to_truncate) # process the standards within the file
        
        # plot a cq graph with standards included
        plt.cq_w.std <- 
          plot_facetted_assay(.data =  polished_cq.dat, 
                              .yvar_plot = 40-CT, .xvar_plot = assay_variable, 
                              .xaxis.label.custom = axislabel.assay_variable) + # plot 40 - Cq
          theme(plot.title = element_text(hjust = 0.5), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) # Axis labels vertical
        
        # Plot Tm1 graph with standards included
        plt.tm1_w.std <- 
          plot_facetted_assay(.data = polished_cq.dat, 
                              .yvar_plot = Tm1, .xvar_plot = assay_variable, 
                              .xaxis.label.custom = axislabel.assay_variable) + 
          theme(plot.title = element_text(hjust = 0.5), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) # Axis labels vertical
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
      unnest(cols = c(w.copy.data)) %>%  # expand the absolute copy number data list
      
      # append the ID of the std curve and equation used
      mutate('std_curve id' = std_to_retrieve)
    
    
    # plot absolute copies per ul template
    plt.copies <- plot_facetted_assay(.data = absolute_dat, .yvar_plot = Copies.per.ul.template, 
                                      .xaxis.label.custom = axislabel.assay_variable)
    
    plt.copies_w.mean <- plot_facetted_assay(.data = absolute_dat, 
                                             .yvar_plot = Copies.per.ul.template, 
                                             .xaxis.label.custom = axislabel.assay_variable, 
                                             points_plt.style = 'jitter') + 
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)
    
    # same plot with x-labels in a single line
    plt.copies_w.mean_straight <- plot_facetted_assay(.data = absolute_dat, 
                                            .xvar_plot = assay_var.horz_label, 
                                            .xaxis.label.custom = axislabel.assay_variable,
                                            .yvar_plot = Copies.per.ul.template, 
                                            points_plt.style = 'jitter') + 
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)
      
  }
  
  # # normalizing copy #s to backbone ----  
  # if(plot_mode == 'absolute_quantification' & plot_normalized_backbone == 'yes')
  # { # computing ratio of copy #s of targets : flipped and unflipped to the backbone
  #   
  #   sel <- absolute_dat %>% select(Sample_category,assay_variable,`Primer pair`,Target,`Copies.per.ul.template`) 
  #   # select relevant columns (other numeric columns will throw errors)
  #   
  #   sel_b <- sel %>% filter(Target == 'Backbone') # filter out each Target
  #   sel_f <- sel %>% filter(Target == 'Flipped'); sel_u <- sel %>% filter(Target == 'Unflipped');
  #   
  #   # make ratios to the backbone
  #   sel_f %<>% mutate("Normalized Copies.per.ul.template" = sel_f$`Copies.per.ul.template`/sel_b$`Copies.per.ul.template`);  
  #   sel_u %<>% mutate("Normalized Copies.per.ul.template" = sel_u$`Copies.per.ul.template`/sel_b$`Copies.per.ul.template`);
  #   
  #   results_ratio <- bind_rows(sel_f, sel_u) # bind results into 1 tibble (for easy plotting)
  #   
  #   # plotting the normalized copy #'s
  # plt_norm <- results_ratio %>% 
  # ggplot(aes(x = `assay_variable`,
  #            y = `Normalized Copies.per.ul.template`, color = Target)) +   # plotting
  # geom_point(size = 2) + facet_grid(~Sample_category, scales = 'free_x', space = 'free_x') # plot points and facetting
  # 
  #   plt_norm.formatted <- plt_norm %>% format_classic(., title_name, axislabel.assay_variable) %>% format_logscale() 
  #   # formatting plot, axes labels, title and logcale plotting
  #   
  #   print(plt_norm)
  # }
  
}



# Save data output ----

{if(plot_mode == 'absolute_quantification') ( absolute_dat 
)   else forplotting_cq.dat} %>% 
  
  write_csv(.,
            str_c('excel files/processed_data/', flnm, '-processed.csv', sep = ''),
            na = '')

# Plotting into html -----------------------------------------------------------------------


# calling r markdown file
rmarkdown::render('make_html_plots.rmd', output_file = str_c('./qPCR analysis/', title_name, '.html'))
