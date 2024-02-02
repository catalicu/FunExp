#title: "FunExp_MultivAnalysis_v1"
#author: "Dr CG"
#date: "2/2/2024"

# Multivariate Diversity analysis for Fun Exp
## Libraries
library(ggplot2)
library(vegan)
library(plyr)
library(lme4)
library(nlme)
library(gridExtra)
library(reshape2)
library(car)
library(MASS)

## Plot themes
# With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
# Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# Load data
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)

## edits:
#Account for leaf_age as weeks 
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))

# Remove from metadata - replicate what was done in divAnalysis
meta.table2=meta.table[-c(
  which(meta.table$leaf_age=='A00'), 
  which(meta.table$leaf_age=='A08'), 
  which(meta.table$leaf_age=='A09'),
  which(meta.table$leaf_age=='A10'), 
  which(meta.table$leaf_age=='BIB'),
  which(meta.table$treatment=='F50'), 
  which(meta.table$treatment=='water_control'),
  which(meta.table$treatment=='broth_control'), 
  which(meta.table$treatment=='TAI_control'), 
  which(meta.table$treatment=='TAO_control')
),]
data.frame(dim(meta.table), dim(meta.table2)) # check that the column was eliminated by looking at the dimensions of the datasets

# Remove from ASVtable 
ASVtable2=ASVtable[-c(
  which(meta.table$leaf_age=='A00'), 
  which(meta.table$leaf_age=='A08'), 
  which(meta.table$leaf_age=='A09'),
  which(meta.table$leaf_age=='A10'), 
  which(meta.table$leaf_age=='BIB'),
  which(meta.table$treatment=='F50'), 
  which(meta.table$treatment=='water_control'),
  which(meta.table$treatment=='broth_control'), 
  which(meta.table$treatment=='TAI_control'), 
  which(meta.table$treatment=='TAO_control')
),]
data.frame(dim(ASVtable), dim(ASVtable2)) 
# check that the column was eliminated by looking at the dimensions of the datasets

# perMANOVA
adonis1=adonis2(ASVtable2~leaf_age_weeks*treatment, data=meta.table2) 
(adonis1)
adonis2=adonis2(ASVtable2.t.wtzero~leaf_age_weeks, data=meta.table2) 
(adonis2)

# calculate NMDS scores per sample
NMDS.mod1=metaMDS(ASVtable2, k=3, distance ='bray', trymax = 50)
meta.table2.nmds=data.frame(meta.table2, scores(NMDS.mod1)$sites)

## ordielipse function
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

## plot NMDS - simple
ggplot(meta.table2.nmds, aes(NMDS1, NMDS2)) + geom_point(size=3, aes(fill=(leaf_age_weeks), shape=treatment)) + Theme + scale_shape_manual(values=c(21,22)) +scale_fill_gradient(low='blue', high='red') 

## plot NMDS - with ordiellipses
# Prepare for the ordihull
NMDS.table = data.frame(NMDS1 = meta.table2.nmds$NMDS1, 
                        NMDS2 = meta.table2.nmds$NMDS2,
                        group=factor(meta.table2.nmds$leaf_age))
NMDS.mean=aggregate(NMDS.table[,1:2],list(group=NMDS.table$group),mean)
df_ell <- data.frame()

for(g in levels(NMDS.table$group)){
  #g='A10'
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS.table[NMDS.table$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}

ggplot(data = NMDS.table, aes(NMDS1, NMDS2)) + geom_point(size=3,shape=21, aes(fill=(group))) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.5, linetype=2)+
  annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=2) + Theme + 
  ggtitle('automatic colors')

ggplot(meta.table2.nmds, aes(NMDS1, NMDS2)) + 
  geom_point(size=3, aes(fill=(leaf_age_weeks), shape=treatment)) + Theme + 
  scale_shape_manual(values=c(21,22)) +scale_fill_gradient(low='blue', high='red') + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.5, linetype=2) + 
  annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=2) +
  ggtitle('color formatting from first figure')
