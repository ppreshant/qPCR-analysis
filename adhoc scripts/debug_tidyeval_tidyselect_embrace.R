.df <- absolute_dat

copies_colname <- 'Copies.per.ul.template' %in% colnames(.df) %>% 
  if(.) 'Copies.per.ul.template' else 'Copies_proportional'

# copies_quocolname <- 'Copies.per.ul.template' %in% colnames(.df) %>% 
#   if(.) quo(Copies.per.ul.template) else quo(Copies_proportional)

copies_quocolname <- 'Copies.per.ul.template' %in% colnames(.df) %>% 
  if(.) expr(Copies.per.ul.template) else expr(Copies_proportional)


copies_quocolname %>% class

givecol <- function(.df, clnm)
{
  select(.df, clnm) %>% head(5)
  
}

givecol(absolute_dat, copies_quocolname)


pivcol <- function(.df = absolute_dat, nmcl, valcl = Copies.per.ul.template)
{
  # clquo <- enquo(valcl)
  
  pivot_wider(.df, 
              id_cols = c(Sample_category, biological_replicates, matches('assay_var')),
              names_from = {{nmcl}}, values_from = any_of(valcl))
  
}

pivcol(nmcl = Target_name, valcl = copies_colname) %>% view
pivcol(copies_quocolname)
