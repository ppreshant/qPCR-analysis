# Script to analyze raw data and multicomponent data of qPCR files by hand
# Author : Prashant
# data = 7 May 2020

# Initialization ----
# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel); library(plotly)

# Enter file name
flnm <- 'Std7_N genes_IDT_25-5-20'  # set the filename : q-S023b GFP 18-2-20
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path
data_range <- c(35,15780) # column range in the sheet where raw data is located
# Functions ----

# Function to choose range for each sheet
make_range <- function(x,y,d_range = data_range)
{
  str_c(x, d_range[1],':',y , d_range[2])
}

# data loading and processing ----
# load file and gather specific data range from specific sheet
rawdat <- read_excel(flpath, sheet = 'Raw Data', range = make_range('A','G')) # reading raw data from 4 channels
multidat <- read_excel(flpath, sheet = 'Multicomponent Data', range = make_range('A','F')) # reading extrapolated data for each fluorophore
procdat <- read_excel(flpath, sheet = 'Amplification Data', range = make_range('A','F')) %>% select(-`Target Name`) # reading extrapolated data for each fluorophore

# Convert data to long form
combdat_all <- list(raw = rawdat, mult = multidat, proc = procdat) %>% map_dfr( ~gather(., key = 'channel', 'signal', -Well, -'Well Position', -Cycle)) # combine all the channels and fluorophores into a pair of name - value columns (long data)
combdat <- list(mult = multidat, proc = procdat) %>% map_dfr( ~gather(., key = 'channel', 'signal', -Well, -'Well Position', -Cycle)) # combine all the channels and fluorophores into a pair of name - value columns (long data)


# plotting ----

# filtering a single well data = A2
# plotting ggplot (raw data and fluorophore data)
# plt1 <- combdat %>% filter(`Well Position` == 'A2') %>%  ggplot(aes(Cycle, signal, colour = channel)) + geom_point() + geom_line()

# interactive plot (raw data and fluorophore data)
iplt1 <- combdat_all %>% filter(`Well Position` == 'A9') %>%  plot_ly(x = ~Cycle, y = ~signal, color = ~channel, type = 'scatter', mode = 'lines+markers')

# interactive plot (raw data and fluorophore data)
iplt2 <- combdat %>% filter(channel == 'ROX') %>%  plot_ly(x = ~Cycle, y = ~signal, color = ~`Well Position`, type = 'scatter', mode = 'lines+markers')

