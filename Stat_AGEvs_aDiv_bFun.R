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

# Load data
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)

## edit tables
# create leaf age column
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))

## Calculate diversity metrics
richness=specnumber(ASVtable)
shannon=diversity(ASVtable, index='shannon')
simpson=diversity(ASVtable, index='simpson')
evenness=shannon/log(richness)

# create data frame
meta.div=data.frame(meta.table, richness, shannon, simpson, evenness)

# remove controls
meta.div2=meta.div[-c(which(meta.div$treatment=='F50'), 
                      which(meta.div$treatment=='water_control'), 
                      which(meta.div$treatment=='broth_control'), 
                      which(meta.div$treatment=='TAI_control'), 
                      which(meta.div$treatment=='TAO_control')),]
data.frame(dim(meta.div), dim(meta.div2)) # check that the column was eliminated by looking at the dimensions of the datasets

# Statistical models for diversity
mod.rich.null=lm(log(richness)~1, data=meta.div2)
mod.rich=lm(log(richness)~leaf_age_weeks, data=meta.div2)
mod.rich1a=lm(log(richness)~treatment, data=meta.div2)
mod.rich2=lm(log(richness)~leaf_age_weeks*treatment, data=meta.div2)
mod.rich3=lm(log(richness)~leaf_age_weeks+treatment, data=meta.div2)
mod.rich4=lm(log(richness)~treatment, data=meta.div2, family=Gamma(link='log'))
AIC(mod.rich.null, mod.rich, mod.rich1a,  mod.rich2, mod.rich3, mod.rich4)
anova(mod.rich.null, mod.rich)
Anova(mod.rich)

#Test the first models:
# lowest AIC
AIC(mod.rich.null, mod.rich, mod.rich2, mod.rich3,  mod.rich3b, mod.rich4,mod.rich4b, 
    mod.rich5, mod.rich6, mod.rich7, mod.rich8)
# test against null
anova(mod.rich.null.gamma, mod.rich7, test='Chisq')
# test the best model against the null
Anova(mod.rich7)

#### Function
# load data
# full diversity
div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
colnames(div.table)[2]='Treatment'
# function diversity
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)

# edit tables
fundiv.table=left_join(fun.table, div.table)
fundiv.table$leaf_age_weeks2=as.numeric(fundiv.table$leaf_age_weeks)
fundiv.table=fundiv.table[-which(fundiv.table$W_Change==min(fundiv.table$W_Change)),]
fundiv.table2=fundiv.table[-which(is.na(fundiv.table$SampleID)),]

## Statistical models for function
mod2.age.null=lm(W_Change~1, data=fundiv.table2)
mod2.age1=lm(W_Change~leaf_age_weeks, data=fundiv.table2)
mod2.age2=lm(W_Change~treatment, data=fundiv.table2)
mod2.age3=lm(W_Change~leaf_age_weeks*treatment, data=fundiv.table2)
mod2.age4=lm(W_Change~leaf_age_weeks+treatment, data=fundiv.table2)
mod2.age5=glm(W_Change~leaf_age_weeks, data=fundiv.table2, family=Gamma(link='log'))
mod2.age6=glm(W_Change~leaf_age_weeks*treatment, data=fundiv.table2, family=Gamma(link='log'))
mod2.age7=glm(W_Change~leaf_age_weeks+treatment, data=fundiv.table2, family=Gamma(link='log'))

#Test the models:
# lowest AIC
AIC(mod2.age.null, mod2.age1, mod2.age2, mod2.age3, mod2.age4, mod2.age5, mod2.age6, mod2.age7)
# test against each other
anova(mod2.age.null, mod2.age5)

# test the best model against the null
Anova(mod2.age5, test.statistic = 'Wald')
