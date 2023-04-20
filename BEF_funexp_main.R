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

ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)

div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
colnames(div.table)[2]='Treatment'

fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)

# Merge diversity and functional data and make sure weeks are numeric
fundiv.table=left_join(fun.table, div.table)
fundiv.table$leaf_age_weeks2=as.numeric(fundiv.table$leaf_age_weeks)

# remove rows with no diversity data for BEF analyses
fundiv.table_BEF1=fundiv.table[-which(is.na(fundiv.table$SampleID)),]
fundiv.table_BEF=fundiv.table_BEF1[-which(fundiv.table_BEF1$W_Change==min(fundiv.table_BEF1$W_Change)),]

## Figure BEF main
# create lines for each treatment separately:
mod.TAI=lm(W_Change~(richness), data=fundiv.table_BEF[which(fundiv.table_BEF$treatment=='TAI'),])
mod.TAO=lm(W_Change~(richness), data=fundiv.table_BEF[which(fundiv.table_BEF$treatment=='TAO'),])

ggplot(fundiv.table_BEF, aes(log(richness), W_Change)) + geom_jitter(color='black',shape=21, aes(size=leaf_age_weeks, fill=treatment))+ Theme2 + xlab('log(number of ASV)') + ylab('Weight change (g)') + scale_fill_manual(values=c('black','white')) 
#  geom_line(aes(log(richness), predict(mod.TAI)), linetype='dashed', data=fundiv.table2[which(fundiv.table2$treatment=='TAI'),]) + 
#  geom_line(aes(log(richness), predict(mod.TAO)), linetype='dotted', data=fundiv.table2[which(fundiv.table2$treatment=='TAO'),])
