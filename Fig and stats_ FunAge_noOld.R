###
# Figure and stats for W_Change ~ age leaf
# REMOVE OLD LEAVES

# Libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(lme4)
library(nlme)
library(gridExtra)
library(car)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# load data
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)
div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
colnames(div.table)[2]='Treatment'
fundiv.table=left_join(fun.table, div.table)
fundiv.table2= fundiv.table[-which(fundiv.table$leaf_age_weeks==24),]

# Statistics
mod.noold=lm(W_Change ~ leaf_age_weeks, data=fundiv.table2)
Anova(mod.noold)

# plot
fun.plot=ggplot(fundiv.table2, aes(leaf_age_weeks, W_Change)) + 
  geom_jitter(size=3, color='black', shape=21, aes(fill=treatment)) + 
  Theme2 +  scale_fill_manual(values=c('black', 'white')) + 
  ylab('Weight change (g)') + xlab('Leaf age (weeks)')  + 
  geom_smooth(method='gam', se=FALSE, color='black', linetype='dashed')

quartz(width=4, height=3.5)
plot(fun.plot)
quartz.save('Figures/Fig_time_fun_noold.png', type='png', dpi=300)
