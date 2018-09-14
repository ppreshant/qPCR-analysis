# Read in the file and do manual analysis and plotting
#---
# Source the general_functions file before running this
# saved important/general analyses codes from console

# choose file name, starts in the same directory as Rproject
flnm <- 'excel files/Chk2c MHT.xls'  

fl <- readqpcr(flnm) # read excel file exported by Quantstudio

# Factorise the sample name in the order to be plotted : the avg and Stdev is already calculated by quantstudio
# separate(fcl,value,c('a','b'), ' ')
fl$Results$`Sample Name` <- fl$Results$`Sample Name` %>% factor(levels = unique(.))

# plot the Tm ; Graph will now show
plttm <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = Tm1) + geom_point(color = 'red', size = 1) +
  theme_classic() + scale_color_brewer(palette="Set1") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
  ggtitle('qPCR on minipreps Chk1: Melting curves') 

# print(plttm)

# plot the CT mean along with replicates
plt <- fl$Results %>% ggplot(.) + aes(x = `Sample Name`, y = CT) + geom_point(color = 'red', size = 1, show.legend = T) +
  geom_boxplot(aes(x = `Sample Name`, y = `Ct Mean`), show.legend = T) +
  theme_classic() + scale_color_brewer(palette="Set1") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 90, hjust = 1, vjust = .3)) + 
  ggtitle('qPCR - minipreps on Gibson reaction: Chk2c') 

print(plt)

# Save manually - copy paste this command
# ggsave('qPCR analysis/Chk1.png', dpi = 600)