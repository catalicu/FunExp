# Libraries
library(ggplot2)
library(vegan)
library(plyr)
library(lme4)
library(nlme)
library(gridExtra)
library(reshape2)
library(car)
library(MASS)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# load datasets
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
TAXtable=read.table('input_data/TAXtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)

# edit datasets for multivariate analyses
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))

# Remove baseline data:
baseline.rows=c(which(meta.table$leaf_age=='BIB'),
  which(meta.table$treatment=='F50'), 
  which(meta.table$treatment=='water_control'),
  which(meta.table$treatment=='broth_control'), 
  which(meta.table$treatment=='TAI_control'), 
  which(meta.table$treatment=='TAO_control'))
# from metadata
meta.table2=meta.table[-baseline.rows,]
# check that the column was eliminated by comparing tables
data.frame(dim(meta.table), dim(meta.table2)) 

# from ASVtable 
ASVtable2=ASVtable[-baseline.rows,]
# check that the column was eliminated by comparing tables
data.frame(dim(ASVtable), dim(ASVtable2)) 

# identify and remove taxa with 0 abundances (present in baselines). 
ASVtable2.t=data.frame(t(ASVtable2))
class(ASVtable2.t)

ASVtable2.t.wtzero=ASVtable2.t[-c(which(rowSums(ASVtable2.t.n)==0)), ]
range(rowSums(ASVtable2.t.wtzero))

# perMANOVA
adonis1=adonis2(ASVtable2~leaf_age_weeks*treatment, data=meta.table2) 
adonis2=adonis2(t(ASVtable2.t.wtzero)~leaf_age_weeks*treatment, data=meta.table2) 

# calculate NMDS scores per sample
NMDS.mod1=metaMDS(ASVtable2.t.wtzero, k=3, distance ='bray', trymax = 50)
meta.table2.nmds=data.frame(meta.table2, scores(NMDS.mod1))

# prepare for ordination plots: code for ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Plot
ggplot(meta.table3.nmds, aes(NMDS1, NMDS2)) + geom_point(size=3, aes(fill=(leaf_age_weeks), shape=treatment)) + Theme + scale_shape_manual(values=c(21,22)) +scale_fill_gradient(low='blue', high='red') 

# Prepare for the ordihull
NMDS.table = data.frame(NMDS1 = meta.table3.nmds$NMDS1, NMDS2 = meta.table3.nmds$NMDS2,group=factor(meta.table3.nmds$leaf_age))
NMDS.mean=aggregate(NMDS.table[,1:2],list(group=NMDS.table$group),mean)
df_ell <- data.frame()

for(g in levels(NMDS.table$group)){
  #g='A10'
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS.table[NMDS.table$group==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}

# full plot
ggplot(data = NMDS.table, aes(NMDS1, NMDS2)) + geom_point(size=3,shape=21, aes(fill=(group))) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.5, linetype=2)+
  annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=2) + Theme + ggtitle('automatic colors')

ggplot(meta.table3.nmds, aes(NMDS1, NMDS2)) + geom_point(size=3, aes(fill=(leaf_age_weeks), shape=treatment)) + 
  Theme + scale_shape_manual(values=c(21,22)) +scale_fill_gradient(low='blue', high='red') + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.5, linetype=2) + 
  annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=2) + 
  ggtitle('color formatting from first figure')

# distance analysis
dists=vegdist(ASVtable2, method='bray')
groups.age = meta.table2$leaf_age
groups.treat = meta.table2$treatment

dists.age=betadisper(dists,  groups.age)
anova(dists.age)

dists.treat=betadisper(dists,  groups.treat)
anova(dists.treat)