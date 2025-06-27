#title: "Fig1_BEF_funexp_main_bDiv_cFun.R"
#author: "Dr CG"
#date:"6/28/2025"

# Description:
# This script performs statistical testing to support patterns portrayed in 
# figure 1, including a) the BEF relationship, b) change in diversity with leaf
# age and c) change in function with leaf age. 

#########################
## Stats ##
#########################

# Load libraries
library(dplyr)
library(vegan)

# Load tables
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
#TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)

div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
colnames(div.table)[2]='Treatment'

fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)

# Merge diversity and functional data and make sure weeks are numeric
fundiv.table=left_join(fun.table, div.table)
fundiv.table$leaf_age_weeks2=as.numeric(fundiv.table$leaf_age_weeks)

# remove rows with no diversity data for BEF analyses
fundiv.table_BEF1=fundiv.table[-which(is.na(fundiv.table$SampleID)),]
fundiv.table_BEF2=fundiv.table_BEF1[-which(fundiv.table_BEF1$W_Change==min(fundiv.table_BEF1$W_Change)),]

# is the highest diversity point an outlier?
boxplot(fundiv.table_BEF1$richness, main = "Boxplot of MPG", col = "lightblue")
  library(outliers)
  grubbs.test(fundiv.table_BEF1$richness)
fundiv.table_BEF=fundiv.table_BEF2[-which(fundiv.table_BEF2$richness==(max(fundiv.table_BEF2$richness))),]

#### GLM div ~ richness * leaf age * protozoan treatment
# The stats can be assessed as a linear model on leaf age and protozoan 
# treatment as a factor.
mod.rich.null=lm(W_Change~1, data=fundiv.table_BEF2)
mod.rich.null.gamma=glm(W_Change~1, data=fundiv.table_BEF2, family=Gamma(link='log'))
mod.rich=lm(W_Change~log(richness+1)*leaf_age_weeks*treatment, data=fundiv.table_BEF)
mod.rich2=lm(W_Change~(richness)*leaf_age_weeks*treatment, data=fundiv.table_BEF)
mod.rich3=lm(W_Change~(richness)*treatment, data=fundiv.table_BEF)
mod.rich3b=lm(W_Change~log(richness)*treatment, data=fundiv.table_BEF)
mod.rich4=lm(W_Change~(richness)*leaf_age_weeks, data=fundiv.table_BEF)
mod.rich4b=lm(W_Change~log(richness)*leaf_age_weeks, data=fundiv.table_BEF)
mod.rich5=lm(W_Change~(richness), data=fundiv.table_BEF)
mod.rich6=lm(W_Change~log(richness), data=fundiv.table_BEF)
mod.rich7=glm(W_Change~(richness)*(leaf_age_weeks), data=fundiv.table_BEF2, family=Gamma(link='log'))
mod.rich8=glm(W_Change~log(richness)*(leaf_age_weeks)*treatment, data=fundiv.table_BEF, family=Gamma(link='log'))
mod.rich9=glm(W_Change~log(richness)*(leaf_age_weeks)+treatment, data=fundiv.table_BEF, family=Gamma(link='log'))

#Test the models:
# lowest AIC
AIC(mod.rich.null, mod.rich, mod.rich2, mod.rich3,  mod.rich3b, mod.rich4,mod.rich4b, mod.rich5, mod.rich6, mod.rich7, mod.rich8)
# test the best modelmod.rich7# test the best model against the null
Anova(mod.rich7)

#### GLM div ~ evenness * leaf age * protozoan treatment
# same model for evenness - for sup mat
mod.eve.null=lm(W_Change~1, data=fundiv.table_BEF)
mod.eve=lm(W_Change~log(evenness+1)*leaf_age_weeks*treatment, data=fundiv.table_BEF)
mod.eve2=lm(W_Change~(evenness)*leaf_age_weeks*treatment, data=fundiv.table_BEF)
mod.eve3=lm(W_Change~(evenness)*treatment, data=fundiv.table_BEF)
mod.eve3b=lm(W_Change~log(evenness+1)*treatment, data=fundiv.table_BEF)
mod.eve4=lm(W_Change~(evenness)*leaf_age_weeks, data=fundiv.table_BEF)
mod.eve4b=lm(W_Change~log(evenness+1)*leaf_age_weeks, data=fundiv.table_BEF)
mod.eve5=lm(W_Change~(evenness), data=fundiv.table_BEF)
mod.eve6=lm(W_Change~log(evenness), data=fundiv.table_BEF)

#Test the models:
# lowest AIC
AIC(mod.eve.null, mod.eve, mod.eve2, mod.eve3,  mod.eve3b, mod.eve4,mod.eve4b, mod.eve5, mod.eve6)
# the null model had lowest AIC

#########################
## Change in a) diversity and b) function with leaf age ##
#########################

## DIVERSITY
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
mod.rich4=glm(log(richness)~treatment, data=meta.div2, family=Gamma(link='log'))

#Test the models:
AIC(mod.rich.null, mod.rich, mod.rich1a,  mod.rich2, mod.rich3, mod.rich4)
Anova(mod.rich)

## FUNCTION
# load data
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
# test the best model against the null
Anova(mod2.age5, test.statistic = 'Wald')

