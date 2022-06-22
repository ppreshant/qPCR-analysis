# g.12 facetted assay dot plot
# Prashant, 2/April/21

# Will start with certain specific features. 
# More hard-coded variabels can be user-input as use cases arrive

plot_facetted_assay <- function(.data = forplotting_cq.dat,  # data.frame or tibble
                                # required variable .yvar_plot
                                .yvar_plot, .xvar_plot = assay_var.label, # variables within the dataframe (not string, use raw name)
                                .colourvar_plot = Sample_category,
                                .facetvar_plot = Target_name, # facets align along the x axis, cannot be NULL as of now
                                .xaxis.label.custom = plot_assay_variable,
                                .subtitle.plot = NULL, # 'y axis label' supplied as a string to .subtitle.plot
                                points_plt.style = 'point') # jitter for scattered effect, else lined up points
  
{
  
  # check missing variables
  if(is.null(enquo(.yvar_plot))) stop('Missing y axis variable for plot in plot_facetted_assay() call')  
  
  # Get plot subtitle (translated from the table if entry exists)
  if(is.null(.subtitle.plot)) .subtitle.plot <- yaxis_translation[deparse(enexpr(.yvar_plot))]
  if(is.na(.subtitle.plot)) .subtitle.plot <- deparse(enexpr(.yvar_plot)) # use full text if not in the translation table
  
  # plotting
  {ggplot(.data, aes(x = {{.xvar_plot}}, y = {{.yvar_plot}}, colour = {{.colourvar_plot}})) +
      
      {if(points_plt.style == 'jitter') {geom_jitter(height = 0, width = .2) # plot points in 'jitter' or normal fashion
        } else geom_point() } +
      
      facet_grid(cols = vars({{.facetvar_plot}}), scales = 'free_x', space = 'free_x') +  # facets
      
      ggtitle(title_name, subtitle = .subtitle.plot) + # custom title and y label as subtitle
      {if(!is.null(.xaxis.label.custom)) xlab(.xaxis.label.custom)} + # custom label for X axis if specified
      {if(!is.na(.subtitle.plot)) ylab('')} }  %>%  # remove y axis label if a subtitle is provided which is more readable
    print()
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