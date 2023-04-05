# Making a heatmap for qPCR data by plate layout
# Prashant
# 30-10-22

#' Make a heatmap for qPCR data based on well positions
#' useful for quick spotting of cross contamination
qPCR_data_to_heatmap <- function(.df = forplotting_cq.dat, 
                                 wellposition = `Well Position`,
                                 .target = 'flipped',
                                 .importance_column = `Sample_category`,
                                 title_name = base_title_name)
{
  
  make_rows_cols(.df, {{wellposition}}) %>% # Separate wells into rows and cols and freeze proper order
    
    filter(str_detect(Target_name, .target)) %>% # filter one target -- TODO : vectorize/loop?
    
    # re-order levels : H -> A 
    mutate(across(row, fct_rev)) %>%  # reverse order of rows for plotting
  
    # mutate(across(row, ~ match(.x, LETTERS[1:26]))) # convert the letter to number
    mutate(importance = if_else(str_detect({{.importance_column}}, 
                                           regex('Uninduced|U$|Control|negative|ntc|Water', ignore_case = TRUE)), 
                                1, 0.5)) %>% 
    
    {ggplot(data = ., 
            aes(x = column, y = row, fill = 40-CT, 
                alpha = importance)) + 
        
        geom_tile(width = 0.95, height = 0.95) + # heatmap with spacing
        
        scale_fill_viridis_c() +  # colour scale
        scale_alpha_continuous(breaks = c(0.5, 1), limits = c(0, 1)) + 
        
        geom_text(aes(label = round(40-CT, 0))) + # show the 40-Cq - no decimals 
        
        ggtitle(glue::glue(title_name, " : ", .target)) # add title
        
        } 
}



#' Make a ~heatmap for qPCR data based on well positions, for all targets
#' useful for quick spotting of cross contamination
qPCR_multitarget_heatmap <- function(.df = forplotting_cq.dat, 
                                 wellposition = `Well Position`,
                                 .importance_column = `Sample_category`,
                                 title_name = base_title_name)
{
  
  make_rows_cols(.df, {{wellposition}}) %>% # Separate wells into rows and cols and freeze proper order
      
    # mutate(across(row, ~ match(.x, LETTERS[1:26]))) # convert the letter to number
    mutate(importance = if_else(str_detect({{.importance_column}}, 
                                           regex('Uninduced|Control|negative|ntc|Water', ignore_case = TRUE)), 
                                1, 0.5)) %>% 
    
    {ggplot(data = ., 
            aes(x = 1, y = 40 - CT, fill = Target_name, colour = Target_name,
                alpha = importance)) + 
        
        geom_bar(stat = 'identity', position = position_dodge2(width = 0.9)) + # heatmap with spacing
        
        # scale_fill_viridis_c() +  # colour scale
        scale_fill_brewer(palette = 'Set1') +
        scale_alpha_continuous(breaks = c(0.5, 1), limits = c(0, 1)) + 
        
        geom_text(aes(label = round(40-CT, 0), y = 42 - CT), 
                  size = 2,
                  position = position_dodge2(width = 0.9), vjust = 0, hjust = 0.5) + # show the 40-Cq - no decimals 
        
        theme(axis.text = element_blank(), axis.ticks = element_blank())+ # hide axis labels
        xlab('') + ylab('')+
        
        ggtitle(title_name, subtitle = 'bars showing 40-Cq') +  # add title
        
        facet_grid(rows = vars(row), cols = vars(column), switch = 'y') + 
        
        # Add grey background where data is present (highlights wells with no amplifications)
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), alpha = 0.05, colour = NA)
      
    } 
}

# TODO : move text into the bar (depending on bar height?) and change colour ~ importance -- ensure dodge still works
# TODO : add label for emphasized wells - at the bottom of the bars?

#' Separate qPCR data Wells into rows and columns in alphabetical and numerical order
make_rows_cols <- function(.df = forplotting_cq.dat, 
                           wellposition = `Well Position`)
{
  # separate the well position into rows and columns
  separate(.df,
           {{wellposition}},
           into = c('row', 'column'),
           sep = 1) %>% 
    
    # order columns in numerical order : 1 -> 12
    mutate(across(column, fct_inseq)) %>% # arrange numerical columns sequentially
    
    # order levels : A -> H 
    mutate(across(row, ~ as.factor(.x))) # base R `as.factor` arranges alphabetical order
}