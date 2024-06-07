# g.12 facetted assay dot plot
# Prashant, 2/April/21

# Will start with certain specific features. 
# More hard-coded variables can be user-input as use cases arrive

# TODO : add a "not detected" text if no points are present for a particular category? - issues with multi cat placement, colours

plot_facetted_assay <- function(.data = forplotting_cq.dat,  # data.frame or tibble
                                
                                # required variables
                                
                                # variable for y axis
                                .yvar_plot, .xvar_plot = assay_var.label, # variables within the dataframe (not string, use raw name)
                                
                                .colourvar_plot = Sample_category, # variable for colour
                                .filter_colourvar = '.*', # wrapper for plotting a subset of colours
                                
                                # facetting
                                .facetvar_plot = Target_name, # facets align along the x axis, cannot be NULL as of now
                                facet_scale_constraint = 'free_x', # 'free_x' or 'free_y' or 'free' or 'fixed'
                                
                                # optional variables
                                .label.var = NULL, # variable for labels on the plot (for ggplotly interactive)
  
                                .xaxis.label.custom = NULL, # pass in plot_assay_variable if needed
                                .subtitle.plot = NULL, # 'y axis label' supplied as a string to .subtitle.plot
                                
                                points_plt.style = 'jitter', # jitter for scattered effect, else 'point' for lined up points
                                jitter_width = 0.05, # the width of jitter (0 - 0.5), adjust depending on # of assay_variables
                                flipped_plot = TRUE) # if you want the quantity on x axis and labels on y

  
{
  
  # check missing variables ----
  if(is.null(enquo(.yvar_plot))) stop('Missing y axis variable for plot in plot_facetted_assay() call')  
  
  # Skip plotting if column not present in data : BUG : fails for 40 - CT ! 
  # .yvar_string <- quo_name(enquo(.yvar_plot)) # get the .y variable as a string
  # if(!.yvar_string %in% colnames(.data)) {
  #   print(glue::glue('No column with name : "{.yvar_string}" found in the data, skipping the plot'))
  #   return() # stop the plotting function
  #   }
  
  # make labels ----
  # Get plot subtitle (translated from the table if entry exists)
  if(is.null(.subtitle.plot)) .subtitle.plot <- yaxis_translation[deparse(enexpr(.yvar_plot))]
  if(is.na(.subtitle.plot)) .subtitle.plot <- deparse(enexpr(.yvar_plot)) # use full text if not in the translation table
  
  
  # prefilter data for specific colours -- typically sample_categories
  quo_colourvar <- enquo(.colourvar_plot) # get the colour variable as a quosure
  
  if(quo_is_null(quo_colourvar) | .filter_colourvar == '.*')
      filtered_data <- .data # no filtering needed (since there is no variable to colour by / no filter)
  
  else filtered_data <- filter(.data, str_detect({{.colourvar_plot}}, .filter_colourvar)) # prefilter data for specific colours -- typically sample_categories
  
  
  # plotting ----
  {ggplot(filtered_data, 
          
          aes(x = {{.xvar_plot}}, y = {{.yvar_plot}}, colour = {{.colourvar_plot}},
              label = {{.label.var}})) + # label for ggplotly hover
      
      {if(points_plt.style == 'jitter') {geom_jitter(height = 0, width = jitter_width) # plot points in 'jitter' or normal fashion
      } else geom_point() } +
      
      # formatting for flipped plots : Warning x and y still refer to the original orientation (not flipped)
      {if(flipped_plot) {list(
        coord_flip(), 
        
        facet_grid(rows = vars({{.facetvar_plot}}), 
                   scales = facet_scale_constraint, space = facet_scale_constraint), # flipped facets
        
        theme(legend.position = 'top'), # put legend on the top
        xlab(''))
        
        # facetting and ylab for regular plots
      } else list(
        
        facet_grid(cols = vars({{.facetvar_plot}}), 
                   scales = facet_scale_constraint, space = facet_scale_constraint), # facets regular
        
        {if(!is.na(.subtitle.plot)) ylab('')}, # remove y axis label if a subtitle is provided which is more readable
        {if(!is.null(.xaxis.label.custom)) xlab(.xaxis.label.custom)}) # custom label for X axis if specified
      } + 
      
      # titles and axis names
      ggtitle(title_name, subtitle = .subtitle.plot) # custom title and y label as subtitle
      
    
  } %>% print()
}



# Convenience function to plot cq after filtering data

# plots 40-CT data
filter_to_cq_plot <- function(text = 'lysate',  # filtering text
                              filter_var = Sample_category, # the variable to apply filter on
                              horizontal = TRUE) # if you want the plot to show horizontal -- good for long labels
  
{
  plot_facetted_assay(.data = forplotting_cq.dat %>% 
                        filter(str_detect({{filter_var}}, text)),
                      .yvar_plot = 40 - CT) +
    
    {if(horizontal) 
      
      list(coord_flip(), 
           facet_grid(rows = vars(Target_name)),
           theme(legend.position = 'top'))  
      
    }
  
}