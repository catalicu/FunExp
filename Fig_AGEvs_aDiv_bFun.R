###
# Figure leaf age vs (a) diversity  and (b) function

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
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)
fundiv.table=left_join(fun.table, div.table)

## format data
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))
# Remove controls from metadata
meta.table2=meta.table[-c(which(meta.table$treatment=='F50'),
                          which(meta.table$treatment=='water_control')),]
# Remove controls from ASVtable 
ASVtable2=ASVtable[-c(which(meta.table$treatment=='F50'),
                      which(meta.table$treatment=='water_control')),]

# Calculate diversity metrics
richness=specnumber(ASVtable)
shannon=diversity(ASVtable, index='shannon')
simpson=diversity(ASVtable, index='simpson')
evenness=shannon/log(richness)
meta.div=data.frame(meta.table, richness, shannon, simpson, evenness)

# remove controls
meta.div2=meta.div[-c(which(meta.div$treatment=='F50'), 
                      which(meta.div$treatment=='water_control'), 
                      which(meta.div$treatment=='broth_control'), 
                      which(meta.div$treatment=='TAI_control'), 
                      which(meta.div$treatment=='TAO_control')),]
# check that the column was eliminated by looking at the dimensions of the datasets
data.frame(dim(meta.div), dim(meta.div2)) 
# Merge function and diversity table
fundiv.table2=left_join(fundiv.table, meta.table2)

# Figures
rich.plot=ggplot(meta.div2, aes(leaf_age_weeks, log(richness))) + 
  geom_jitter(size=3, color='black', shape=21, aes( fill=treatment))+  
  geom_smooth(method='gam', se=FALSE, color='black', linetype='dashed') + 
  ylab('log(ASV richness') + xlab('Leaf age (weeks)')  + 
  xlab('') + Theme2 + scale_fill_manual(values=c('black', 'white')) + 
  theme(legend.position = 'none') +
  annotate('text', x=0, y=5.5, label='a.')

fun.plot=ggplot(fundiv.table2, aes(leaf_age_weeks, W_Change)) + 
  geom_jitter(size=3, color='black', shape=21, aes(fill=treatment)) + 
  Theme2 +  scale_fill_manual(values=c('black', 'white')) + 
  ylab('Weight change (g)') + xlab('Leaf age (weeks)')  + 
  geom_smooth(method='gam', se=FALSE, color='black', linetype='dashed') +
  annotate('text', x=0, y=0.035, label='b.')

quartz(width=8, height=3.7)
Age_divfun=arrangeGrob(rich.plot, fun.plot, ncol=2)
plot(Age_divfun)
quartz.save('Figures/Fig_time_ADiv_BFun.png', type='png', dpi=300)
