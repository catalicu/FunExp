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

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(linewidth = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))



# A.The most functional species who is related to time
otumelt2_ASV35.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV35'),]
ggplot(otumelt2_ASV35.div, aes(leaf_age_weeks2, (relative_abundance+1))) + geom_point() + 
  geom_smooth(method='lm', se=FALSE, color='black') +
  Theme + annotate('text', x=1, y=150000, label='b.')


# B. 
#### Do the most abundant taxa increase over time? No
otumelt2_ASV1.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV1'),]
ggplot(otumelt2_ASV1.div, aes(richness, (Abundance+1))) + geom_point() + 
  geom_smooth(method='lm', se=FALSE, color='black')

# A. 
#### does the most functional taxa increase over time? Yes
otumelt2_ASV369.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV369'),]
ggplot(otumelt2_ASV369.div, aes(richness, (Abundance+1))) + geom_point() + geom_smooth(method='lm', se=FALSE, color='black')


# B. 
#### Do the most abundant taxa increase over time?
otumelt2_ASV1.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV1'),]
ggplot(otumelt2_ASV1.div, aes(leaf_age_weeks2, (relative_abundance+1))) + 
  geom_point() + geom_smooth(method='lm', se=FALSE, color='black') +
  Theme + annotate('text', x=1, y=150000, label='b.')

# C. 
#### are there relationships between abundance, occurrence and function? No
ASVrel.abund.pa.fun.sorted=read.table('output_tables/ASVrel.abund.pa.txt', header=TRUE)

ASVrel.abund.pa_plot_c=ggplot(ASVrel.abund.pa.fun.sorted, aes(log(ASVrel.abund), ASVpa.freq)) + 
  geom_point(aes(fill=Slope), size=3, shape=21) + Theme + 
  scale_fill_manual(values=c('black', 'grey', 'white', 'blue')) + 
  ylab('ASV occurrence') + xlab('ASV Relative Abundance') + 
  annotate('text', x=1, y=1, label='c.')