#title: "FunASV_plots"
#author: "Dr CG"
#date:"1/13/2023"
# taken from:

# Description:
# This script creates figures to identify ASV underlying function:
  # * ASV relative abundances vs function
  # * General patterns of ASV distribution and occurrences

# Libraries

# Libraries
library(ggplot2)
library(dplyr)
library(reshape2)
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

# A.# does the most functional taxa increase over time? Yes
otumelt2_ASV35.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV35'),]
Fig6a_ASV35=ggplot(otumelt2_ASV35.div[which(otumelt2_ASV35.div$Abundance!=0),], 
       aes(leaf_age_weeks2, log(Abundance+1))) + geom_point() + 
  geom_smooth(method='lm', se=FALSE, color='black') + Theme +
  xlab('Leaf age (weeks)') + ylab('ASV35 standardized abundance') +
  annotate('text', x=1, y=8, label='a.')

# B. 
#### Do the most abundant taxa increase over time? No
otumelt2_ASV1.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV1'),]
Fig6b_ASV1=ggplot(otumelt2_ASV1.div[which(otumelt2_ASV1.div$Abundance!=0),], aes(leaf_age_weeks2, log(Abundance+1))) + geom_point() + 
  geom_smooth(method='lm', se=FALSE, color='black') + Theme +
  xlab('Leaf age (weeks)') + ylab('ASV1 standardized abundance') +
  annotate('text', x=1, y=12, label='b.')


# C. to few datapoints 
#### does the most functional taxa increase over time? Yes
otumelt2_ASV691.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV691'),]
Fig6c_ASV691=ggplot(otumelt2_ASV691.div[which(otumelt2_ASV691.div$Abundance!=0),], 
       aes(leaf_age_weeks2, log(Abundance+1))) + geom_point() + 
  geom_smooth(method='lm', se=FALSE, color='black') + Theme +
  xlab('Leaf age (weeks)') + ylab('ASV691 standardized abundance') +
  annotate('text', x=1, y=4, label='c.')

# D. 
#### are there relationships between abundance, occurrence and function? No
ASVrel.abund.pa_plot_d=ggplot(ASVrel.abund.pa.fun.sorted, aes(log(ASVrel.abund), ASVpa.freq)) + 
  geom_point(aes(fill=Slope), size=3, shape=21) + Theme + 
  scale_fill_manual(values=c('black', 'grey', 'white', 'blue')) + 
  ylab('ASV occurrence') + xlab('ASV Relative Abundance') + 
  annotate('text', x=1, y=1, label='d.')

Fig6=arrangeGrob(Fig6a_ASV35, Fig6b_ASV1, 
                 Fig6c_ASV691, ASVrel.abund.pa_plot_d, ncol=2)
quartz()
plot(Fig6)
