# formatting plots

# plot formatting ---- 


# Set theme universally : format as classic, colours = Set1
theme_set(theme_classic()) # theme
scale_colour_discrete <- function(...) { # palette
  scale_colour_brewer(..., palette="Set1")
}

# plot formatting function : format as classic, colours = Set1
format_classic <- function(plt)
{ # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
  plt <- plt +
    theme_classic() + scale_color_brewer(palette="Set1")
}

# plot formatting function : format as classic, colours = Set1
format_classic <- function(plt)
{ # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
  plt <- plt +
    theme_classic() + scale_color_brewer(palette="Set1")
}

# plot formatting function : format as logscale
format_logscale_y <- function(plt)
{ # extra comments
  plt <- {plt +
    scale_y_log10(  # logscale for y axis with tick marks
      labels = fancy_scientific
      #labels = scales::trans_format("log10", scales::math_format(10^.x) )
    )} %>% print()
}

# plot formatting function : format as logscale x
format_logscale_x <- function(plt)
{ # extra comments
  plt <- {plt +
    scale_x_log10(  # logscale for y axis with tick marks
      labels = fancy_scientific
      #labels = scales::trans_format("log10", scales::math_format(10^.x) )
    )} %>% 
    print()
}


# formatting labels in logscale cleanly : a x 10^b
# source: https://stackoverflow.com/a/24241954/9049673

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+","e",l)
  
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  
  # convert 1x10^ or 1.000x10^ -> 10^ 
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  
  # return this as an expression
  parse(text=l)
}

# use as ggplot(df,aes(x,y)) + geom_point() + scale_y_log10(labels = fancy_scientific)



#' Scientific notation with 10^ instead of e
#' @param x numeric vector
#' @return expression vector with 10^ instead of e
#' @examples see below
# source: https://stackoverflow.com/a/18530540/9049673

scientific_10 <- function(x, output_format = 'string') 
{   
  parse(text=gsub("e\\+*", " %*% 10^", 
                  scales::scientific_format()(x))) %>% 
    
    {if(output_format == 'string')
    {map_chr(as.list(.), deparse)  # return as string
    
      } else parse(text = .)  # return as expression
    }
  
    # directly output string : not elegantly formatted
  # gsub("e\\+*", " x 10^", 
  #      scales::scientific_format()(x)) %>%
  
}

# using as geom_text label ; need to map_chr a deparse to convert to string vector for the above two formatting fns

# example
# b <- tibble::tibble(b1 = -2:2, b2 = 10^b1)

# ggplot(b, 
#        aes(x = b1, y = b2)) + 
#   geom_point() + 
#   geom_text(aes(label = scientific_10(b2) #%>% 
#                 # {map_chr(as.list(.), deparse)}
#   ),
#   vjust = -2, parse = TRUE)
