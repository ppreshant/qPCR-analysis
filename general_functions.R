# FUnctions to load qPCR data and manipulate it. The functions can be called from another R file

# read in excel file (.xls) of qPCR exported from Quantstudio 3 
  # Make sure to include raw data as well

# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel)  

# reading files and manipulating columns ----

# read in the excel file (from row 44 onwards)
readqpcr <- function(flnm)
{
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = 43)
  
  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}

order_columnwise <- function(fl)
{
  fl$Results$`Well Position` %>%  str_sub(2) %>% as.integer() %>% order()
}

# standard curve and regressions ----
# Plot Standard curve
plotstdcurve <- function(fl, plttitle, xlabel)
{
  plt <- fl$Results %>% filter(Task == 'STANDARD') %>% 
ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target Name`) + geom_point() +
theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
ggtitle(plttitle) + xlab(xlabel) +
stat_smooth(method ="lm", se = F)
}

# getting regression line and R2 values to put into the standard curve plot
lm_eqn <- function(df, trig = 0){
  x = df %>% pull(Quantity) %>% log10 ; y = df %>% pull(CT)
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x)+ a*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  if(trig == 'coeff') c(round(coef(m)[2], 2), round(summary(m)$r.squared, 2))
  else as.character(as.expression(eq)); 
}

  
  optional1 <- function()
    {# output the difference between consecutive CT values
    tsumrev <- trev %>% group_by(`Sample Name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
    diff(tsumrev$CT) %>% round(2)}

