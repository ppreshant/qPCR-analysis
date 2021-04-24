# g.12 facetted assay dot plot
# Prashant, 2/April/21

# Will start with certain specific features. 
# More hard-coded variabels can be user-input as use cases arrive

plot_facetted_assay <- function(.data = polished_cq.dat,  # data.frame or tibble
                                # required variable .yvar_plot
                                .yvar_plot, .xvar_plot = assay_variable, # variables within the dataframe (not string, use raw name)
                                .colourvar_plot = Sample_category,
                                .facetvar_plot = Target_name, # facets align along the x axis, cannot be NULL as of now
                                .xaxis.label.custom = plot_assay_variable,
                                .subtitle.plot = NULL) # 'y axis label' supplied as a string to .subtitle.plot
{
  
  # check missing variables
  if(is.null(enquo(.yvar_plot))) stop('Missing y axis variable for plot in plot_facetted_assay() call')  
  
  # Get plot subtitle
  if(is.null(.subtitle.plot)) .subtitle.plot <- yaxis_translation[deparse(enexpr(.yvar_plot))]
  
  # plotting
  {ggplot(.data, aes(x = {{.xvar_plot}}, y = {{.yvar_plot}}, colour = {{.colourvar_plot}})) +
      geom_point() +
      facet_grid(cols = vars({{.facetvar_plot}}), scales = 'free_x', space = 'free_x') +  # facets
      ggtitle(title_name, subtitle = .subtitle.plot) + # custom title and y label as subtitle
      xlab(.xaxis.label.custom) + # custom label for X axis
      {if(!is.na(.subtitle.plot)) ylab('')} }  %>%  # remove y axis label if a subtitle is provided which is more readable
    print()
}
