# Read in the file and do manual analysis and plotting
# Author: Prashant Kalvapalle;  October 10 2018

source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name


# Labeling translators ----

# check 16-user_parameters.R for modifying these options to make plots more readable

# yaxis_translation : (replaces the y axis names to be readable in the plots)
# lst_assay.vars_translation : attach explanations to assay_variable appearing on x axis (ex: plasmid numbers) 


axislabel.assay_variable <- 'Template name' # printed on the x axis of the graph


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
  
  
  # plot 40 - Cq
  plt.cq <- plot_facetted_assay(.yvar_plot = 40-CT, .xaxis.label.custom = axislabel.assay_variable, 
                                flipped_plot = FALSE)
  
  horz.cq <- plot_facetted_assay(.yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                                 .xaxis.label.custom = axislabel.assay_variable)
  
  
  # plot ~ copies (2^40-Cq : not representative copy number)
  plt.copies <- plot_facetted_assay(.yvar_plot = Copies_proportional,
                                    .xaxis.label.custom = axislabel.assay_variable)
  
  # This is a standby plot in the absence of absolute copy number calculation
  
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
      
      std_par <- 
        {if(template_source == 'googlesheet')
          read_sheet(sheeturls$plate_layouts_PK, sheet = 'qPCR Std curves', range = 'A:G', col_types = 'ccnnnnn') else {
            read_csv(file = 'qPCR analysis/Standards/qPCR_Std_curve_parameters.csv')}
        } %>% 
        
      filter(str_detect(ID, std_to_retrieve ))
    
      } else {
        
        std_par <- NULL # initialize a dummy std_par
        
        # if it is a standard curve holding file (Stdx), call standard curve processor
        if(str_detect(flnm, 'Std[:digit:]*')) 
        { # process the standards within the file
          std_par <- process_standard_curve(flnm, polished_cq.dat, dilutions_to_truncate) 
        
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
    
      }
    
    # If no standard curve is found, then throw an informative error
    if(is.null(std_par)) 
      stop(str_c('No Std curve information found in the file name:',
                 flnm,
                 '\n Run with plot_mode <-  raw_quantification if you dont want Std curve information '))

        
    # calculate absolute copies ----
    
    # Use the std curve parameters to back-calculate the absolute copies
    absolute_dat <- forplotting_cq.dat %>%
      select(-Copies_proportional) %>% # remove the dummy copies data
      
      group_by(Target_name) %>%
      nest() %>% # create a new column with data frames for each target
      
      # calculate copy number for each dataset, and the mean for replicates
      summarize(w.copy.data = map2(data, Target_name,  
                                    ~ absolute_backcalc(.x, .y, std_par) ) 
      ) %>% 
      unnest(cols = c(w.copy.data)) %>%  # expand the absolute copy number data list
      
      # append the ID of the std curve and equation used
      mutate('std_curve id' = std_to_retrieve)
    
    
    # plot absolute copies ----
    
    # plot absolute copies per ul template
    plt.copies <- plot_facetted_assay(.data = absolute_dat, .yvar_plot = Copies.per.ul.template, 
                                      .xaxis.label.custom = axislabel.assay_variable, flipped_plot = F)
    
    plt.copies_w.mean <- plot_facetted_assay(.data = absolute_dat, 
                                             .yvar_plot = Copies.per.ul.template, 
                                             .xaxis.label.custom = axislabel.assay_variable, 
                                             points_plt.style = 'jitter', flipped_plot = F) + 
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)
    
    # same plot with x-labels in a single line
    horz.copies_w.mean <- {plot_facetted_assay(.data = absolute_dat, 
                                            .xvar_plot = assay_var.horz_label, 
                                            .xaxis.label.custom = axislabel.assay_variable,
                                            .yvar_plot = Copies.per.ul.template, 
                                            points_plt.style = 'jitter') + 
      geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)} %>% 
      
      format_logscale_y()
      
  }
  
  
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
