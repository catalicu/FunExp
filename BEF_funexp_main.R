# Libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(lme4)
library(nlme)
library(gridExtra)
library(car)

# original code in FunExp_BEFunctAnalysis_v2.Rmd

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

BEFmain_fig= ggplot(fundiv.table_BEF, aes(log(richness), W_Change)) + geom_jitter(color='black',shape=21, aes(size=leaf_age_weeks, fill=treatment))+ Theme2 + xlab('log(number of ASV)') + ylab('Weight change (g)') + scale_fill_manual(values=c('black','white')) 
quartz(height=4, width=4.5)
plot(BEFmain_fig)
quartz.save('Figures/Fig_BEFmain.png', type='png', dpi=300)


### stats
# model selected out of options 
# (transformed-log, fixed factors: age, treatment, distributions: gaussian, gamma)
mod.rich.null=lm(W_Change~1, data=fundiv.table_BEF)
mod.rich.null.gamma=glm(W_Change~1, data=fundiv.table_BEF, family=Gamma(link='log'))
mod.rich7=glm(W_Change~(richness)*(leaf_age_weeks), data=fundiv.table_BEF, family=Gamma(link='log'))

# lowest AIC
AIC(mod.rich.null, mod.rich.null.gamma, mod.rich7)
# test against the null
Anova(mod.rich.null.gamma, mod.rich7)
# test the best model against the null
Anova(mod.rich7)


### evenness - for sup mat

# model selected was not significantly different than the null model.
mod.eve.null=lm(W_Change~1, data=fundiv.table_BEF)
mod.eve4=lm(W_Change~(evenness)*leaf_age_weeks, data=fundiv.table_BEF)
