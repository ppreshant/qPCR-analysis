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
    map(read_excel, path = flnm, skip = 43)
  
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

hill_fit <- function(results_array)
{
  # source: https://github.com/dritoshi/Fitting-Hill-equation/blob/master/bin/hill.r
  # Fittiing Hill equation
  # Itoshi NIKAIDO <dritoshi@gmail.com>
  
  # make demo data
  L  <- results_array$assay_variable
  y  <- results_array$mean
  
  # # conf
  # output <- "results/hill.pdf"
  
  # initial
  y0 <- min(y)
  ymax.init <- 1e10
  n.init  <- 1
  Kd.init <- .005
  
  # fitting Hill equation
  y.nls <- nlsLM(y ~ y0 + (ymax - y0) * L^n / (Kd^n + L^n), start = c(ymax = ymax.init, n = n.init, Kd = Kd.init))
  
  # # extract fitting data
  # y.nls.summary <- summary(y.nls)
  # y.nls.n       <- y.nls.summary$param[1]
  # y.nls.Kd      <- y.nls.summary$param[2]
  # y.nls.predict <- predict(y.nls)
  # results <- cbind(y, y.nls.predict)
}
  
  optional1 <- function()
    {# output the difference between consecutive CT values
    tsumrev <- trev %>% group_by(`Sample Name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
    diff(tsumrev$CT) %>% round(2)}

  # plot formatting ---- 
  
  # plot formatting function : format as classic, colours = Set1
  format_classic <- function(plt, title_name, plot_assay_variable)
  { # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
    plt <- plt +
      theme_classic() + scale_color_brewer(palette="Set1") + 
      theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) +
      ggtitle(title_name) + xlab(plot_assay_variable)
  }
  
  # plot formatting function : format as logscale
  format_logscale <- function(plt)
  { # extra comments
    plt <- plt +
      scale_y_log10(  # logscale for y axis with tick marks
        labels = fancy_scientific
        #labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
  }
  
  # plot formatting function : format x axis as logscale
  format_logscale_x <- function(plt)
  { # extra comments
    plt <- plt +
      scale_x_log10(  # logscale for x axis with tick marks
        labels = fancy_scientific
      )
  }
  
  # formatting labels in logscale cleanly : a x 10^b
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