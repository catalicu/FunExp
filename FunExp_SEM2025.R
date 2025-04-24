######################################
## SEM - Multivariate analyses
## 
#######################################

# Preparation

# Libraries
library(lavaan)
library(lavaanPlot)
library(dplyr)
library(ggplot2)
library(semPaths)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))
# load tables: 
div.table=read.table('input_data/MetaDiv_table_FunExp12021-03-30.txt', header=TRUE)
colnames(div.table)[2]='Treatment'
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)
fundiv.table=left_join(fun.table, div.table)
# edit:
## leaf ages as numbers
fundiv.table$leaf_age_weeks2=as.numeric(fundiv.table$leaf_age_weeks)

# remove rows with no diversity data for BEF analyses
fundiv.table_BEF1=fundiv.table[-which(is.na(fundiv.table$SampleID)),]
fundiv.table_BEF2=fundiv.table_BEF1[-which(fundiv.table_BEF1$W_Change==min(fundiv.table_BEF1$W_Change)),]

# scale the variables
fundiv.table_BEF3 =fundiv.table_BEF2
fundiv.table_BEF3[,c(9,15)] =apply(fundiv.table_BEF3[,c(9,15)], 2, scale)

# remove outlier - based on histogram and percentile method
hist(fundiv.table_BEF3$richness)
lower=quantile(fundiv.table_BEF3$richness, 0.025)
upper=quantile(fundiv.table_BEF3$richness, 0.975)
outlier = which(fundiv.table_BEF3$richness < lower | fundiv.table_BEF3$richness > upper)
outlier = c(outlier)
fundiv.table_BEF4=fundiv.table_BEF3[-outlier,]

## Structural Equation modeling
### full model
m1.full = '
  # regressions
  W_Change ~ leaf_age_weeks2 + richness
  richness ~ leaf_age_weeks2
  W_Change ~~ 0*richness
  W_Change ~~ 0*leaf_age_weeks2
  '
fit1.full = sem(m1.full, data=fundiv.table_BEF4)
summary(fit1.full, fit.measures=TRUE)

### model selection approach

# mod2 is the simplest model with all the relationships we care about.
mod2 = '
  # regressions
  W_Change ~ leaf_age_weeks2 + richness
  richness ~ leaf_age_weeks2
'
fit2 = sem(mod2, data=fundiv.table_BEF4)
summary(fit2, fit.measures=TRUE)


# mod3 is the simplest model without richness influence on function
mod3 = '
  # regressions
  W_Change ~ leaf_age_weeks2
  richness ~ leaf_age_weeks2
'
fit3 = sem(mod3, data=fundiv.table_BEF4)
summary(fit3, fit.measures=TRUE)

# mod4 is the simplest model without leaf influence on function
mod4 = '
  # regressions
  W_Change ~ richness
  richness ~ leaf_age_weeks2
'
fit4 = sem(mod4, data=fundiv.table_BEF4)
summary(fit4, fit.measures=TRUE)
#AIC220
#BIC227

## compare models
chi_square_diff <- lavTestLRT(fit4, fit2)
chi_square_diff
chi_square_diff2 <- lavTestLRT(fit1.full, fit2)
chi_square_diff2

#Plot
labels = list(leaf_age_weeks2= 'Leaf age (weeks)', 
              richness='ASV richness', 
              W_Change='Weight loss (mg)')
lavaanPlot(model=fit1.full, labels=labels, node_options= list(shape='box', fontname='Helvetica'), edge_options=list(color='grey'), coefs=TRUE)

