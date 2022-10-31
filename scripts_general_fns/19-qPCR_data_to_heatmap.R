# Making a heatmap for qPCR data by plate layout
# Prashant
# 30-10-22

#' Make a heatmap for qPCR data based on well positions
#' useful for quick spotting of cross contamination
qPCR_data_to_heatmap <- function(.df = forplotting_cq.dat, 
                                 wellposition = `Well Position`,
                                 .target = 'flipped',
                                 .transparency = `Sample_category`,
                                 title_name = base_title_name)
{
  
  # separate the well position into rows and columns
  separate(.df,
           {{wellposition}},
           into = c('row', 'column'),
           sep = 1) %>% 
    
    filter(str_detect(Target_name, .target)) %>% # filter one target -- TODO : vectorize/loop?
    
    # order levels : H -> A 
    mutate(across(row, ~ fct_inorder(.x) %>% fct_rev)) %>%  
  
    # mutate(across(row, ~ match(.x, LETTERS[1:26]))) # convert the letter to number
    mutate(importance = if_else(str_detect({{.transparency}}, 'Uninduced|Control|negative|ntc'), 
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



