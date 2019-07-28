# FUnctions to load qPCR data and manipulate it. The functions can be called from another R file

# read in excel file (.xls) of qPCR exported from Quantstudio 3 
  # Make sure to include raw data as well

# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel)  

# reading files and manipulating columns ----

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = 35)
  
  # convert CT values into numeric 
  class(fl$Results$CT) <- 'numeric'
  fl
}

columnwise_index <- function(fl)
{ # this orders the samples columnwise in the PCR plate or strip : order will give indices for A1,B1,C1 ... A2,B2,C2... (wherever samples exist)
  fl$Results$`Well Position` %>%  str_sub(2) %>% as.integer() %>% order() # Take the well position A1 etc., extract the number, read as integer (instead of char), order it (1,1,1.. 2,2,... 12,12...) and regurn indices 
}

# lookup table of primer pairs and respective targets - not useful as same information is in `Target` already
primer_table <- c('q1-3' = 'Flipped', 'q4-5' = 'Flipped', 
                  'q5-11' = 'Unflipped', 'q1-2' = 'Unflipped', 'q9-10' = 'Unflipped', 'q4-2' = 'Unflipped', 'q2-4' = 'Unflipped', 'q5-11' = 'Unflipped', 'q6-7' = 'Unflipped',
                  'q12-13' = 'Backbone')

# function to back-calculate CT using standard curve parameters
absolute_backcalc <- function(df, std_par)
{
  target_current <- df$Target %>% unique()
  std_current <- std_par %>% filter(target == target_current)
  
  df %>% mutate(`Copy #` = 10^( (CT - std_current$intercept)/std_current$slope) )
}

# standard curve and regressions ----
# Plot Standard curve
plotstdcurve <- function(fl, plttitle, xlabel)
{
  plt <- fl$Results %>% filter(Task == 'STANDARD') %>% 
ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target`) + geom_point() +
theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
ggtitle(plttitle) + xlab(xlabel) + ylab(expression(C[q])) +
stat_smooth(method ="lm", se = F) # plots linear regression line
}

# getting regression line and R2 values to put into the standard curve plot
lm_eqn <- function(df, trig = 0){
  x = df %>% pull(Quantity) %>% log10 ; y = df %>% pull(CT)
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x)+ a*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  # if(trig == 'coeff') c(round(coef(m)[2], 2), round(summary(m)$r.squared, 2))
  if(trig == 'coeff') tibble(slope = round(coef(m)[2], 2), r_square = round(summary(m)$r.squared, 2))
  else as.character(as.expression(eq)); 
}

  
  optional1 <- function()
    {# output the difference between consecutive CT values
    tsumrev <- trev %>% group_by(`Sample Name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
    diff(tsumrev$CT) %>% round(2)}

