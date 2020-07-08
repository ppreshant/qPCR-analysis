# Biobot ID sheet to draw sample names from
bb_sheets <- c('Week 10 (6/15)', 'Week 11 (6/22)', 'Week 12 (6/29)')


# Bring WWTP names from google sheet: "Biobot Sample IDs"
biobot_lookup <- map_df(bb_sheets , ~ read_sheet('https://docs.google.com/spreadsheets/d/1ghb_GjTS4yMFbzb65NskAlm-2Gb5M4SNYi4FHE4YVyI/edit#gid=233791008', sheet = .x, range = 'H1:J40') %>% 
  rename('Biobot ID' = matches('Biobot|Comments', ignore.case = T), 'WWTP' = contains('SYMBOL', ignore.case = T)) %>% 
  mutate(WWTP = as.character(WWTP)) %>% 
    arrange(`Biobot ID`) %>% 
  mutate(Week = .x)
)

bb_wide <- biobot_lookup %>% pivot_wider(names_from = Week, values_from = `Biobot ID`)
