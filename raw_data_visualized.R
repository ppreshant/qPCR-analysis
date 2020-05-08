# Script to analyze raw data and multicomponent data of qPCR files by hand
# Author : Prashant
# data = 7 May 2020

# calling libraries ; make sure they are installed (install.packages)
library(readxl); library(magrittr); library(tidyverse); library(ggrepel); library(plotly)

# Enter file name
flnm <- 'Pr4 hybrid MHT'  # set the filename : q-S023b GFP 18-2-20
flpath <- str_c('excel files/',flnm,'.xls') # this completes the file path

# load file and gather specific data range from specific sheet
rawdat <- read_excel(flpath, sheet = 'Raw Data', range = 'A36:G15780')
multidat <- read_excel(flpath, sheet = 'Multicomponent Data', range = 'A36:E15780')

# Convert data to long form
combdat <- list(raw = rawdat, mult = multidat) %>% map_dfr( ~gather(., key = 'channel', 'signal', -Well, -'Well Position', -Cycle))

# filtering a single well data = A2

# plotting normal
plt1 <- combdat %>% filter(`Well Position` == 'A2') %>%  ggplot(aes(Cycle, signal, colour = channel)) + geom_point() + geom_line()

# interactive plot
combdat %>% filter(`Well Position` == 'A2') %>%  plot_ly(x = ~Cycle, y = ~signal, color = ~channel, type = 'scatter', mode = 'lines+markers')

