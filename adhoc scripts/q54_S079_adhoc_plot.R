# q54_S079_adhoc_plot

# Prereq ----

source('./0-general_functions_main.R') # Source the general_functions file before running this

source('./0.5-user-inputs.R') # source the user inputs from a different script


# load data ----

flnm <- 'q54_S079-RNA2_Travis-1-5_RAM' # overwrites 
title_name <- 'q54_S079_RNA2-RAM' # make title_name for plot and saving

absolute_dat <- get_processed_datasets(flnm)

# labels ----

# regex based translation to change x-axis labels (named vector)
plasmid_translation <- c('502' = '(-)',
                         '553' = 'Ec', 
                         '554' = 'So',
                         '555' = 'Vn',
                         '556' = 'Ab',
                         '557' = 'Pp',
                         '558' = 'So-1'
                         )


# data mix ----

abs2 <- 
  
  filter(absolute_dat, Sample_category != 'Travis') %>% # remove travis data

# translate to fullish name for plasmid
  mutate(across(assay_var.horz_label, 
                ~ str_replace_all(.x, plasmid_translation) %>% # change labels
                  fct_relevel(plasmid_translation) # order them as above!
                ) ) 


## Ratios to 16S ----

# Modified code from 2H_all_organisms.R script of RAM
# Calculating spliced fraction of the Ribozyme/16s
ratio_wider <- 
  
  abs2 %>% 
  
  # select only desired columns (for pivoting)
  select(Target_name,  assay_variable, assay_var.horz_label, Copies.per.ul.template, biological_replicates) %>%  # select only required columns
  rename(organism = assay_var.horz_label) %>%  # rename column to organism
  
  pivot_wider(names_from = Target_name, values_from = Copies.per.ul.template, names_prefix = 'copies_') %>%  # each target to a col
  
  # calculate ratios of each targets per replicate
  mutate(spliced_per_ribozyme = copies_U64 / copies_gfpbarcode,
         spliced_per_16S = copies_U64 / `copies_16s`, 
         expresssion_per_16S = copies_gfpbarcode / `copies_16s`) %>% 
  
  # find the mean across biological replicates for all numeric cols
  group_by(assay_variable) %>% 
  mutate(across(where(is.numeric), 
                \ (x) mean(x, ignore.na = TRUE),
                .names = "mean_{.col}")) 



# plotting ---- 


## all targets ----

horz.copies_w.mean <- {plot_facetted_assay(.data = abs2, 
                                           
                                           .filter_colourvar = 'Neg|Test',
                                           
                                           .xvar_plot = assay_var.horz_label, 
                                           .xaxis.label.custom = axislabel.assay_variable,
                                           .yvar_plot = Copies.per.ul.template, 
                                           points_plt.style = 'jitter') + 
    geom_boxplot(aes(y = mean_Copies.per.ul.template), show.legend = FALSE)} %>% 
  
  format_logscale_y()

# save plot
ggsave(plot_as('q54_S079_RNA2'), width = 7, height = 5)
    

## Expression vs splicing ----


plt.positivesonly_correlation <- 
  {ggplot(
    
    # remove unnecessary data for this plot
    filter(ratio_wider, !str_detect(organism, 'NTC|(-)|552')), 
    
    aes(x = `copies_gfpbarcode`, y = `copies_U64`,
        label = organism)) + 
      
      geom_point() + 
      
      # add labels on the points with organism name
      ggrepel::geom_text_repel(nudge_y = 1) + 
      
      # ggforce::geom_mark_rect(aes(label = organism)) + labelling points does not work
      geom_smooth(method = 'lm', show.legend = FALSE, alpha = 0.2) + # draw a correlation line
      
      
      labs(title = title_name)
    
  } %>%
  
  format_logscale_x() %>%
  format_logscale_y()

# save plot
ggsave(plot_as(title_name, '_expr vs splice'), width = 3, height = 4)


### Ratio to 16S ----

plt_correlation_normalize16s <- 
  {ggplot(
    
    # remove unnecessary data for this plot
    filter(ratio_wider, !str_detect(organism, 'NTC|(-)|552')), 
    
    aes(x = expresssion_per_16S, y = `spliced_per_16S`,
        label = organism)) + 
      
      geom_point() + 
      
      # add labels on the points with organism name
      ggrepel::geom_text_repel(nudge_y = 1) + 
      
      # ggforce::geom_mark_rect(aes(label = organism)) + labelling points does not work
      geom_smooth(method = 'lm', show.legend = FALSE, alpha = 0.2) + # draw a correlation line
      
      
      labs(title = title_name)
    
  } %>%
  
  format_logscale_x() %>%
  format_logscale_y()


ggplotly(plt_correlation_normalize16s, dynamicTicks = T)

# save plot
ggsave(plot_as(title_name, '_expr vs splice_per16S'), width = 3, height = 4)
