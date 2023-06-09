# 20-plot_dose_response_and_controls.R

# dose_response + controls plotter ----

# dose_response plot / flipped
plot_dose_response_and_controls <- function(.data = forplotting_cq.dat, # use ratio_data for fractions
                                        .target_to_filter = 'flipped', # ignore if plotting fractions
                                        .yvar = 40 - CT, # 40 - CT or flipped_fraction
                                        .xvar_dose = Arabinose, .xlabel = 'Arabinose (uM)',
                                        .xvar_control = assay_variable, # assay_variable/qPCR ; Samples/plate reader
                                        output_all_plots = F)
  
{
  
  # data subset ----
  data_subset <- ungroup(.data) %>% 
    filter({{.xvar_control}} != 'water', # remove water samples
           if_any(any_of('Target_name'), ~ .x == .target_to_filter)) # filter if column is present..!
  
  
  # y-axis limits ----
  yrange <- reframe(data_subset, range({{.yvar}}, na.rm = T))
  ymin <- floor(yrange[[1,1]]) # round down
  ymax <- yrange[[2,1]] * 1.05 # count up by 5% ; 
    # old: ceiling(yrange[[2,1]]) ; plyr::round_any(., 0.5, f = ceiling) # round up to nearest 0.4
  
  # TODO : remove ceiling with user input?
  
  # dose response ----
  ara_plt <- 
    {ggplot(filter(data_subset, sample_type == 'Induction'),
            
            aes({{.xvar_dose}}, {{.yvar}}, label = biological_replicates, 
            label2 = `Well Position`)) + # for interactive troubleshooting
        
        geom_point() + 
        
        geom_line(aes(group = biological_replicates), alpha = 0.2) +
        
        theme(legend.position = 'top') + 
        
        ylim(c(ymin, ymax)) + # set consistant yaxis ranges
        # labels and titles
        xlab(.xlabel) + 
        ggtitle('qPCR memory with arabinose conc.', subtitle = title_name)
      
    } %>% 
    
    format_logscale_x() # format_logscale_y()
  
  
  # controls -----
  control_plt <- 
    ggplot(filter(data_subset, sample_type == 'Controls'),
           
           aes({{.xvar_control}}, {{.yvar}}, label = biological_replicates, 
               label2 = `Well Position`)) + # for interactive troubleshooting
    
    ylim(c(ymin, ymax)) + # set consistant yaxis ranges
    geom_point(position = position_jitter(width = 0.2, height = 0)) + 
    ylab(NULL) + xlab('Controls')
  
  
  # merge plots -----
  library(patchwork)
  
  # attach panels [dose response x 4 + controls x 1 widths]
  combined_plt <- 
  ara_plt + control_plt + 
    plot_layout(widths = c(4, 1))
  
  # output a list of plots if user input asks
  if(output_all_plots) list(combined_plt, ara_plt, control_plt) else combined_plt
  
}

