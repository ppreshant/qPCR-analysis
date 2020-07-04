# FUnctions to load qPCR data and manipulate it. The functions can be called from another R file

# read in excel file (.xls) of qPCR exported from Quantstudio 3 
  # Make sure to include raw data as well

# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel); library(googlesheets4) 

# reading files and manipulating columns ----

# read in the excel file (from row 36 onwards)
readqpcr <- function(flnm)
{
  fl <- flnm %>%  
    excel_sheets() %>% 
    set_names(.,.) %>% 
    map(read_excel, path = flnm, skip = 38)
  
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

read_plate_to_column <- function(data_tibble, val_name)
{ # transforms a plate reader table into a column (named after the top left cell, unless mentioned)
  # eliminates plate row,column numbering ; Select 1 row above the plate (even if it doesn't contain a label)
  
  val_name <- enquo(val_name)
  # colnames(data_tibble) <- data_tibble[1,] # set column names as the first row
  data_tibble[,] %>% gather(key = 'col_num', value = !!val_name, -`<>`) %>% rename(row_num = `<>`) %>% unite('Well Position', c('row_num', 'col_num'), sep = '') %>% drop_na()}

# mutates a subset of data and returns a new array (works for multiple conditions)
mutate_when <- function(data, ...) 
{ # Source: Stackoverflow - https://stackoverflow.com/a/34170176/9049673
  
  dots <- eval(substitute(alist(...)))
  for (i in seq(1, length(dots), by = 2)) 
  {
    condition <- eval(dots[[i]], envir = data)
    mutations <- eval(dots[[i + 1]], envir = data[condition, , drop = FALSE])
    data[condition, names(mutations)] <- mutations
  }
  data
}


# mutates a subset of data and returns a new array (does multiple mutations on same condition)
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) 
{ # Source: Stackoverflow -  https://stackoverflow.com/a/34096575/9049673
  
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# Gets the 96 well layout with template names matching the experiment ID from filename in a google sheet
get_template_for <- function(bait, sheet_url = templates_sheet)
{ # Looking for WWx or Stdx - example WW21 or Std7 within the filename; Assumes plate spans from row B to N (1 row below the matching ID)
 
   # Finding the plate to be read
  plate_names_row <- read_sheet(sheet_url, sheet = 'Plate layouts', range = 'C:C', col_types = 'c')
  m_row <- plate_names_row %>% unlist() %>% as.character() %>% 
    # find the row with standard beginings matching the filename
    str_detect(., bait %>% str_match('^(WW|Std|dd.WW)[:alnum:]*') %>% .[1] ) %>% 
    which() + 1
  range_to_get <- str_c('B', m_row + 1, ':N', m_row + 9)
  
  # read the template corresponding to the file name
  plate_template_raw <- read_sheet(templates_sheet, sheet = 'Plate layouts', range = range_to_get)
  
  # Convert the 96 well into a single column, alongside the Well position
  plate_template <- read_plate_to_column(plate_template_raw, 'Sample Name') # convert plate template (sample names) into a single vector, columnwise
  
}

# standard curve and regressions ----

# Plot Standard curve
plotstdcurve <- function(results_qpcr, plttitle, xlabel)
{
  plt <- results_qpcr %>% 
ggplot(.) + aes(x = log10(Quantity), y = CT, color = `Target`) + geom_point() +
theme_classic() + scale_color_brewer(palette="Set1") + theme(plot.title = element_text(hjust = 0.5)) + 
ggtitle(plttitle) + xlab(xlabel) + ylab(expression(C[q])) +
stat_smooth(data = filter(results_qpcr, Task == 'STANDARD'), method ="lm", se = F) # plots linear regression line
}

# getting regression line and R2 values to put into the standard curve plot 
# Source: https://stackoverflow.com/a/7549819/9049673

# least square fitting - linear regression for standard curve
lm_std_curve <- function(df, trig = 0)
{
  x = df %>% pull(Quantity) %>% log10 ; y = df %>% pull(CT)
  m <- lm(y ~ x, df)
  lm_eqn(m, trig)
}


# linear regression equation
lm_eqn <- function(m, trig = 0){
  
  eq <- substitute(italic(y) == b %.% italic(x)+ a*","~~italic(r)^2~":"~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  # if(trig == 'coeff') c(round(coef(m)[2], 2), round(summary(m)$r.squared, 2))
  if(trig == 'coeff') tibble(slope = round(coef(m)[2], 2), y_intercept = round(coef(m)[1], 2), r_square = round(summary(m)$r.squared, 2))
  else as.character(as.expression(eq)); 
}

  
  optional1 <- function()
    {# output the difference between consecutive CT values
    tsumrev <- trev %>% group_by(`Sample Name`) %>% summarise(CT = mean(CT), Quantity = mean(Quantity), CT_sd = sd(CT))
    diff(tsumrev$CT) %>% round(2)}

# Tm plots ----
# plotting functions for Melting temperature

# plots all the Tm's if samples have multiple peaks in the melting curve
plotalltms <- function(results_relevant)
{ 
  # Gather the Tm's into another data frame and merge into 1 column
  tmfl <- results_relevant %>% select(`Sample Name`, starts_with('Tm')) %>% gather('Peak number','Tm',-`Sample Name`)
  
  # plot the Tm ; Graph will now show
  plttm2 <- tmfl %>% ggplot(.) + aes(x = `Sample Name`, y = Tm) + geom_point(aes(color = `Peak number`), size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}

# plot the first Tm only ; Graph will now show
plottm1 <- function(results_relevant)
{ 
  plttm <- results_relevant %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
    theme_classic() + scale_color_brewer(palette="Set1") + 
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
    ggtitle(paste(title_name,': Melting curves')) + facet_grid(~`Primer pair`, scales = 'free_x', space = 'free_x')
}


# Plotting functions ----


# Plotting mean, sd and individual replicates jitter
plot_mean_sd_jitter <- function(summary_data = summary_results, raw_data = results_abs, long_format = F, measure_var = 'Copy #', sample_var = '.*', exclude_sample = F, colour_var = Target, x_var = assay_variable, y_var = `Copy #`, facet_var = `Sample Name`, title_text = title_name, ylabel = 'Genome copies/ul RNA', xlabel = plot_assay_variable)
{ # Convenient handle for repetitive plotting in the same format; Reads data only in long format or wide (specify in long_format)
  
  # filtering variables by user inputs
  if(long_format) # use long format if not plotting Copy #s - ex. Recovery, % recovery etc.
  {
    summ_relevant <- summary_data %>% filter(Measurement == measure_var, str_detect(`Sample Name`, sample_var, negate = exclude_sample))
    raw_relevant <- raw_data %>% filter(Measurement == measure_var, str_detect(`Sample Name`, sample_var, negate = exclude_sample))
    y_var <- sym('val') # default y variable is val
  } else 
    {
      summ_relevant <- summary_data %>% filter(str_detect(`Sample Name`, sample_var, negate = exclude_sample))
      raw_relevant <- raw_data %>% filter(str_detect(`Sample Name`, sample_var, negate = exclude_sample))
    }
  
  plt1 <- summ_relevant %>% ggplot(aes(x = {{x_var}}, y = mean, colour = {{colour_var}})) +
    geom_point(size = 2) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) +
    
    # Individual data points
    geom_jitter(data = raw_relevant, aes(y = {{y_var}}, alpha = map_chr({{y_var}}, ~. == 0), size = map_chr({{y_var}}, ~. == 0)), width = .2, show.legend = F ) +
    scale_alpha_manual(values = c(.3, 1)) + scale_size_manual(values = c(1, 2)) + # manual scale for emphasizing unamplified samples
    
    # Plotting actual spike ins (only for Recovery plot)
    { if(measure_var == 'Recovered') list(geom_point(data = filter(summary_data, str_detect(Measurement,'Actual'), str_detect(`Sample Name`, sample_var, negate = exclude_sample)), colour = 'black', shape = 21), 
             geom_line(data = filter(summary_data, str_detect(Measurement,'Actual'), str_detect(`Sample Name`, sample_var, negate = exclude_sample)), aes(group = {{colour_var}})))
    } +
    
    # Facetting and labelling
    facet_grid(cols =  vars({{facet_var}}), scales = 'free_x', space = 'free_x') +
    ggtitle(title_text) + ylab(ylabel) + xlab(xlabel)

  plt1.formatted <- plt1 %>% format_classic() # clean formatting
  
}

# plotting individual replicates
plot_biological_replicates <- function(results_abs, title_text = title_name, xlabel = plot_assay_variable)
{ # Convenient handle for repetitive plotting 'Copy #' vs biological replicate
  
  plt <- results_abs %>% ggplot(aes(x = `Tube ID`, y = `Copy #`, color = Target)) + ylab('Copies/ul RNA extract') +    # Specify the plotting variables 
    geom_point(size = 2) + facet_grid(~`Sample Name`, scales = 'free_x', space = 'free_x') + # plot points and facetting
    ggtitle(title_text) + xlab(xlabel)
  plt.formatted <- plt %>% format_classic(.) %>% format_logscale_y() # formatting plot, axes labels, title and logcale plotting
}

# Scatter plot with a linear regression fit and equation
plot_scatter <- function(plot_data = results_abs, long_format = F, measure_var = 'Copy #', sample_var = '.*', exclude_sample = F, colour_var = NULL, x_var = N1_multiplex, y_var = N2_multiplex, title_text = title_name, ylabel = 'Genome copies/ul RNA', xlabel = 'X')
{ # Convenient handle for repetitive plotting in the same format; Reads data only in long format or wide (specify in long_format)
  
  # filtering variables by user inputs
  if(long_format) # use long format if not plotting Copy #s - ex. Recovery, % recovery etc.
  {
    plot_relevant <- plot_data %>% filter(Measurement == measure_var, str_detect(`Sample Name`, sample_var, negate = exclude_sample))
    y_var <- sym('val') # default y variable is val
  } else 
  {
    plot_relevant <- plot_data %>% filter(str_detect(`Sample Name`, sample_var, negate = exclude_sample))
    
  }
  
  # linear regression equation
  lin_reg_eqn <- plot_relevant %>% lm({{y_var}} ~ {{x_var}}, data = .) %>% lm_eqn(.)
  
  plt1 <- plot_relevant %>% ggplot(aes(x = {{x_var}}, y =  {{y_var}}, colour = {{colour_var}})) +
    geom_point(size = 2) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = .1) +
    
    # linear regression
    goem_smooth(method = 'lm') + 
    geom_text(data = . %>% summarise_all(~ max(., na.rm = T)), label = lin_reg_eqn, parse = TRUE, show.legend = F, hjust = 'inward', nudge_x = -5, force = 10)
}

# plot formatting ---- 
 
 
  # plot formatting function : format as classic, colours = Set1
  format_classic <- function(plt)
  { # formats plot as classic, with colour palette Set1, centred title, angled x axis labels
    plt <- plt +
      theme_classic() + scale_color_brewer(palette="Set1")
  }
  
  # plot formatting function : format as logscale
format_logscale_y <- function(plt)
  { # extra comments
    plt <- plt +
      scale_y_log10(  # logscale for y axis with tick marks
        labels = fancy_scientific
        #labels = scales::trans_format("log10", scales::math_format(10^.x) )
      )
  }

# plot formatting function : format as logscale x
format_logscale_x <- function(plt)
{ # extra comments
  plt <- plt +
    scale_x_log10(  # logscale for y axis with tick marks
      labels = fancy_scientific
      #labels = scales::trans_format("log10", scales::math_format(10^.x) )
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