source('./0-general_functions_main.R') # Source the general_functions file before running this

# Source user inputs  ----
# choose file name, title for plots (file name starts in the same directory as Rproject)

source('./0.5-user-inputs.R') # source the user inputs from a different script
title_name <- base_title_name

# Load file ----

absolute_dat <- get_processed_datasets(flnm)


# Plots ----


horz.cq <- plot_facetted_assay(.data = absolute_dat,
                               .yvar_plot = 40-CT, .xvar_plot = assay_var.horz_label, 
                               .xaxis.label.custom = axislabel.assay_variable)


# same plot with x-labels in a single line
horz.copies_w.mean <- {plot_facetted_assay(.data = absolute_dat, 
                                           .xvar_plot = assay_var.horz_label, 
                                           .xaxis.label.custom = axislabel.assay_variable,
                                           .yvar_plot = Copies.per.ul.template, 
                                           points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)} %>% 
  
  format_logscale_y()


# output ----

# calling r markdown file
rmarkdown::render('test.rmd', output_file = str_c('./qPCR analysis/', 'test', '.html'))
