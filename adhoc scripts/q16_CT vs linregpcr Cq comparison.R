# To quality check linregpcr again. Since a few negative samples have ben having suprious Cq calls

# Run the first 50 lines of the 2C_Uxx_variants_copies.R script

.cq <- map_dfr(flnms, 
               # read linregpcr + R processed csv
               ~ read_csv(str_c('excel files/processed_data/', .x , '-processed.csv'),
                          col_types = 'cccncnnnncnnc')
) %>% mutate(across(Sample_category, ~ str_replace(.x, 'negatives', 'negative')))

# get delimiters for hardcoding above :: issue with assay_variable being read as double
tst <- read_csv(str_c('excel files/processed_data/', flnms[2] , '-processed.csv'))
map_chr(tst, class) %>% str_extract('^.') %>% paste(collapse = '')


# merge
.dfall <- left_join(.df, .cq, by = c("Target_name", "Sample_category", "assay_variable")) %>% 
  replace_na(list(CT = Inf))


# plott
pltsctr <- ggplot(.dfall, aes(Cq, CT, colour = Sample_category, label1 = Target_name, label2 = assay_variable)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + annotate(geom = 'text', x = 17, y = 15, label = 'y = x line' ) + 
  theme(legend.position = 'top') #+ 
  # xlim(c(0, 25)) + ylim(c(0,25))
  
pltsctr
ggsave(plot_as('q16_CT vs linregpcr Cq comparison'), width = 4, height = 4)

ggplotly(pltsctr, dynamicTicks = T)

# verify negatives
filter(.dfall, Sample_category == 'negative', !is.na(Cq))
    
filter(.dfall, Sample_category == 'negative') %>% view()


# compare cq extrapolation to N0
{.dfall %>% mutate(inferred_N = 10^(-Cq)) %>% 
  ggplot(aes(inferred_N, N0, colour = Sample_category)) + geom_point()} %>% 
  format_logscale_x() %>%  format_logscale_y()

# Compare Cq to N0 on a logscale
cqn0 <- {ggplot(.dfall, aes(-Cq, N0, colour = Sample_category, label1 = Target_name, label2 = assay_variable)) + geom_point() +
    theme(legend.position = 'top')} %>% 
  format_logscale_y()

cqn0
ggsave(plot_as('q16_linregpcr Cq vs N0 comparison'), width = 4, height = 4)

ggplotly(cqn0, dynamicTicks = T)


# Check all samples that don't have plateau : After the Cq of that of U1/empty vector samples
.dfall %>% filter(Cq >= 34) %>% 
  {ggplot(., aes(-Cq, N0, colour = Sample_category, label1 = Target_name, label2 = assay_variable)) + geom_point() +
      theme(legend.position = 'top')} %>% 
  format_logscale_y()
  
