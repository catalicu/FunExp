#title: "FunASV_plots"
#author: "Dr CG"
#date:"6/28/2025"

# Description:
# This script creates figures to identify ASV underlying function:
  # * ASV relative abundances vs function
  # * General patterns of ASV distribution and occurrences

# Libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(linewidth = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# Load data - from 'Fig6_FunASV_plotprep.R' 
ASVlist2.div=read.table('output_tables/ASVlist2.txt', header=TRUE)
ASVrel.abund.pa.fun.sorted=read.table('output_tables/ASVrel.abund.pa.txt', header=TRUE)

# a. 
#### are there relationships between abundance, occurrence and function? No
ASVrel.abund.pa_plot_d=ggplot(ASVrel.abund.pa.fun.sorted, aes(log(ASVrel.abund), ASVpa.freq)) + 
  geom_point(aes(fill=Slope), size=3, shape=21) + Theme + 
  scale_fill_manual(values=c('black', 'grey', 'white', 'blue')) + 
  ylab('ASV occurrence') + xlab('ASV Relative Abundance') + 
  annotate('text', x=1, y=1, label='a.')

# b. how do these taxa compare in terms of their average contribution to function?
ASVselect=ASVrel.abund.pa.fun.sorted[c(which(ASVrel.abund.pa.fun.sorted$ASVcode=='ASV35'),
                           which(ASVrel.abund.pa.fun.sorted$ASVcode=='ASV1'),
                           which(ASVrel.abund.pa.fun.sorted$ASVcode=='ASV691')),]
avFunASV=ggplot(ASVselect, aes(ASVcode, ASVfunRat)) + geom_col(color='black', aes(fill=ASVcode)) +
  Theme2 + xlab('') + ylab('average weight loss \n standardized by ASV abundance') +
  scale_fill_manual(values=c('black', 'white', 'white'))+ 
  annotate('text', x=0.7, y=1.1, label='b.')


Fig6=arrangeGrob(ASVrel.abund.pa_plot_d, avFunASV, ncol=2)
quartz(width = 9, height=4.5)
plot(Fig6)
quartz.save('Figures/Fig6.png', type='png', dpi=300)
