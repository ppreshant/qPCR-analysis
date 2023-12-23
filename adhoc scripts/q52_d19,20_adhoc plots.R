# q52_d19,20_adhoc plots.R


# D19 data : 1280 series + ntc + MG1655
d19 <- absolute_dat %>% 
  filter(str_detect(assay_variable, '^12..|MG1655|ntc'), # select only the correct ones
         !str_detect(Target_name, '16s') # remove 16S
         )

plt_d19 <- {plot_facetted_assay(.data = d19, 
                                .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = Copies.per.ul.template, 
                                points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + 
    
    ggtitle('q52_D19_new RAM fusions_repression')} %>% 
  
  format_logscale_y()


ggsave(plot_as('q52_D19_zoom'), plt_d19, width = 5, height = 5)


# D19 data : 1280 series + ntc + MG1655
d19 <- absolute_dat %>% 
  filter(str_detect(assay_variable, '^12..|MG1655|ntc'), # select only the correct ones
         !str_detect(Target_name, '16s') # remove 16S
  )

plt_d19 <- {plot_facetted_assay(.data = d19, 
                                .xvar_plot = assay_var.horz_label, 
                                .xaxis.label.custom = axislabel.assay_variable,
                                .yvar_plot = Copies.per.ul.template, 
                                points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE) + 
    
    ggtitle('q52_D19_new RAM fusions_repression')} %>% 
  
  format_logscale_y()


ggsave(plot_as('q52_D19_zoom'), plt_d19, width = 5, height = 5)