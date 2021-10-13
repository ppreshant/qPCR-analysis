# subset dead only q18 -------

deaddat <- filter(forplotting_cq.dat,
                  str_detect(assay_variable, '103|295|54|NTC')) 
  

dead.plt <- plot_facetted_assay(.data = deaddat, .yvar_plot = 40-CT)                

ggsave(str_c('qPCR analysis/Archive/', 'q18_dead U64', '.png'), width = 8, height = 4)



# old vs new q17 negatives - 314 ------

flnms <- c('q017a_all organisms_16s_Std17', 
           'q17b_all organisms_duplx U64 gfp_Std8',   # raw filename, without the "-processed" keywork
           'q18_dead U64 and isolated_4-10-21')

d314s <- map_dfr(flnms,
                 ~ read_csv(file = str_c('excel files/processed_data/', .x, '-processed.csv')) %>% 
                   mutate(run_name = str_extract(.x, 'q[:digit:]*')) )  %>% 
   filter(str_detect(assay_variable, '314|NTC'))

title_name <- 'Compare cross runs q17-18: 314s'
plt.314 <- plot_facetted_assay(.data = d314s, .yvar_plot = 40-CT, .colourvar_plot = run_name, )

ggsave(str_c('qPCR analysis/Archive/', '314s - q17-18', '.png'), width = 12, height = 4)
