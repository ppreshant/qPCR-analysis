# Read in the file and do manual analysis and plotting
#---
# Source the general_functions file before running this
# saved important/general analyses codes from console

# choose file name, starts in the same directory as Rproject
flnm <- 'excel files/Chk2 MHT.xls'  

fl <- readqpcr(flnm) # read excel file exported by Quantstudio

# Separate primer pair from sample name
fl$Results <- separate(fl$Results,`Sample Name`,c('Sample Name','Primer pair'),' ')

# Factorise the sample name in the order for plotting: the avg and Stdev is already calculated by quantstudio
fl$Results$`Sample Name` <- fl$Results$`Sample Name` %>% factor(levels = unique(.))

# Factorise the primer pairs in the order for plotting
fl$Results$`Primer pair` <- fl$Results$`Primer pair` %>% factor(levels = unique(.))

# plot the Tm ; Graph will now show
plttm <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 2) +
  theme_classic() + scale_color_brewer(palette="Set1") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
  ggtitle('qPCR on fresh minipreps: Chk1: Melting curves') + facet_grid(~`Primer pair`)

# print(plttm)

# plot the CT mean along with replicates
plt <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 1, show.legend = T) +
  geom_boxplot(aes(x = `Sample Name`, y = `Ct Mean`), show.legend = T) +
  theme_classic() + scale_color_brewer(palette="Set1") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
  ggtitle('qPCR on fresh minipreps: Chk1') + facet_grid(~`Primer pair`)

print(plt)

# Save manually - copy paste this command
# ggsave('qPCR analysis/Chk1.png', dpi = 600)