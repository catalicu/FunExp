#title: "Fig1_BEF_funexp_main_bDiv_cFun.R"
#author: "Dr CG"
#date edits:"6/27/2025"

# Description:
# This script creates the three panels for figure 1 including a) the BEF
# relationship, b) the change of diversity with leaf age and c) the change
# in function with leaf age. 

# Libraries
library(ggplot2)
library(dplyr)
library(outliers)
library(effects)
library(vegan)
library(patchwork)

#########################
## BEF plots (a) ##
#########################

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)
div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
  colnames(div.table)[2]='Treatment' # create columns with same names for left_join

# Merge diversity and functional data and make sure weeks are numeric
fundiv.table=left_join(fun.table, div.table)
fundiv.table$leaf_age_weeks2=as.numeric(fundiv.table$leaf_age_weeks)

# remove rows with no diversity data for BEF analyses
fundiv.table_BEF1=fundiv.table[-which(is.na(fundiv.table$SampleID)),]
fundiv.table_BEF2=fundiv.table_BEF1[-which(fundiv.table_BEF1$W_Change==min(fundiv.table_BEF1$W_Change)),]

# is the highest diversity point an outlier?
boxplot(fundiv.table_BEF1$richness, main = "Boxplot of MPG", col = "lightblue")
  grubbs.test(fundiv.table_BEF1$richness)
fundiv.table_BEF=fundiv.table_BEF2[-which(fundiv.table_BEF2$richness==(max(fundiv.table_BEF2$richness))),]
fundiv.table_BEF$log_richness=log(fundiv.table_BEF$richness)

# model selected out of options 
mod.rich7=glm(W_Change~log_richness*(leaf_age_weeks), data=fundiv.table_BEF, family=Gamma(link='log'))
effect_data=data.frame(effect('log_richness', mod.rich7)) # predict on one variable of a model

BEFmain_fig= ggplot(fundiv.table_BEF, aes(log_richness, W_Change)) + 
  geom_jitter(color='black',shape=21, aes(size=leaf_age_weeks, fill=treatment))+ 
  geom_line(data=effect_data, aes(x=log_richness, y=fit), 
            linetype='dashed', linewidth=1) +
  Theme2 + xlab('log(number of ASV)') + 
  ylab('Weight change (g)') + 
  annotate(geom='text', y=0.0325, x=2.7, label='a.')+
  scale_fill_manual(values=c('black','white')) 

plot(BEFmain_fig)

### evenness - for sup mat

# model selected was not significantly different than the null model.
mod.eve.null=lm(W_Change~1, data=fundiv.table_BEF)
mod.eve4=lm(W_Change~(evenness)*leaf_age_weeks, data=fundiv.table_BEF)

#########################
## Age plots (b and c) ##
#########################

## format data
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))
# Remove controls from metadata
meta.table2=meta.table[-c(which(meta.table$treatment=='F50'),
                          which(meta.table$treatment=='water_control')),]
# Remove controls from ASVtable 
ASVtable2=ASVtable[-c(which(meta.table$treatment=='F50'),
                      which(meta.table$treatment=='water_control')),]

# Re-Calculate diversity metrics
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
  geom_smooth(method='lm', se=FALSE, color='black', linetype='dashed') + 
  ylab('log(ASV richness)') + xlab('Leaf age (weeks)')  + 
  Theme2 + scale_fill_manual(values=c('black', 'white')) + 
  theme(legend.position = 'none') +
  annotate('text', x=0, y=5.5, label='b.')

fun.plot=ggplot(fundiv.table2, aes(leaf_age_weeks, W_Change)) + 
  geom_jitter(size=3, color='black', shape=21, aes(fill=treatment)) + 
  Theme2 +  scale_fill_manual(values=c('black', 'white')) + 
  ylab('Weight change (g)') + xlab('Leaf age (weeks)')  + 
  geom_smooth(method='gam', se=FALSE, color='black', linetype='dashed') +
  annotate('text', x=0, y=0.035, label='c.')

## Merge all plots together and save
quartz(width=6, height=5)
layout2 = (BEFmain_fig | plot_spacer()) / (rich.plot | fun.plot) 
layout2
quartz.save('Figures/Fig_BEF_time_ADiv_BFun.png', type='png', dpi=300)
