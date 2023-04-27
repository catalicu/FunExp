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

# STATS
# first model
mod.rich.null= lm(W_Change ~ 1, data=fundiv.table_BEF)
mod.rich.null.gamma= glm(W_Change ~ 1, data=fundiv.table_BEF, family=Gamma(link='log'))
mod.rich= lm(W_Change ~ log(richness+1) *leaf_age_weeks*treatment, 
             data=fundiv.table_BEF)
mod.rich2= lm(W_Change ~ (richness) *leaf_age_weeks*treatment, 
              data=fundiv.table_BEF)
mod.rich3= lm(W_Change ~ (richness) *treatment, 
              data=fundiv.table_BEF)
mod.rich3b= lm(W_Change ~ log(richness) *treatment, 
               data=fundiv.table_BEF)
mod.rich4= lm(W_Change ~ (richness) *leaf_age_weeks, 
              data=fundiv.table_BEF)
mod.rich4b= lm(W_Change ~ log(richness) *leaf_age_weeks, 
               data=fundiv.table_BEF)
mod.rich5= lm(W_Change ~ (richness), 
              data=fundiv.table_BEF)
mod.rich6= lm(W_Change ~ log(richness), 
              data=fundiv.table_BEF)
mod.rich7= glm(W_Change ~ (richness) *(leaf_age_weeks), 
               data=fundiv.table_BEF, family=Gamma(link='log'))
mod.rich8= glm(W_Change ~ log(richness) *(leaf_age_weeks)*treatment, 
               data=fundiv.table_BEF, family=Gamma(link='log'))
mod.rich9= glm(W_Change ~ log(richness) *(leaf_age_weeks)+treatment, 
               data=fundiv.table_BEF, family=Gamma(link='log'))
#Test the first models:
# lowest AIC
AIC(mod.rich.null, mod.rich, mod.rich2, mod.rich3,  mod.rich3b, mod.rich4,mod.rich4b, 
    mod.rich5, mod.rich6, mod.rich7, mod.rich8)
# test against null
anova(mod.rich.null.gamma, mod.rich7, test='Chisq')
# test the best model against the null
Anova(mod.rich7)
