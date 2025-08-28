#title: "Fig4_NMDS_adonis.R"
#author: "Dr CG"
#date edits:"6/27/2025"

# Description:
# This script creates the three panels for figure 1 including a) the BEF
# relationship, b) the change of diversity with leaf age and c) the change
# in function with leaf age. 

# Libraries
library(ggplot2)
library(vegan)
library(dplyr)
library(gridExtra)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# load datasets
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
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

ASVtable2.t.wtzero=ASVtable2.t[-c(which(rowSums(ASVtable2.t)==0)), ]
range(rowSums(ASVtable2.t.wtzero))

# perMANOVA
adonis1=adonis2(ASVtable2~leaf_age_weeks*treatment, data=meta.table2) 
adonis2=adonis2(t(ASVtable2.t.wtzero)~leaf_age_weeks*treatment, data=meta.table2) 

# calculate NMDS scores per sample
set.seed(123)
NMDS.mod1=metaMDS(ASVtable2.t.wtzero, k=3, distance ='bray', trymax = 50)
meta.table2.nmds=data.frame(meta.table2, scores(NMDS.mod1)$species)

# prepare for ordination plots: code for ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

# Plot
ggplot(meta.table2.nmds, aes(NMDS1, NMDS2)) + 
  geom_point(size=3, aes(fill=(leaf_age_weeks), shape=treatment)) + 
  Theme + scale_shape_manual(values=c(21,22)) +
  scale_fill_gradient(low='blue', high='red') 

# Prepare for the ordihull
meta.table2.nmds2=meta.table2.nmds[-which(meta.table2.nmds$leaf_age=='A08'),]
NMDS.table2 = data.frame(NMDS1 = meta.table2.nmds2$NMDS1, 
                        NMDS2 = meta.table2.nmds2$NMDS2,
                        group=factor(meta.table2.nmds2$leaf_age))

NMDS.mean=aggregate(NMDS.table2[,1:2],list(group=NMDS.table2$group),mean)
#check for groups with few samples:
NMDS.table3=NMDS.table2[-which(NMDS.table2$group==names(
  which(table(NMDS.table2$group)==min(table(NMDS.table2$group))))),]
NMDS.table3$group=factor(NMDS.table3$group)

df_ell <- data.frame()

for(g in levels(NMDS.table3$group)){
  #g='A10'
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS.table3[NMDS.table3$group==g,],
              veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),
              length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,group=g))
}

# full plot
ggplot(meta.table2.nmds, aes(NMDS1, NMDS2)) + 
  geom_point(size=3, aes(fill=factor(leaf_age_weeks), shape=treatment)) + Theme2 + 
  scale_shape_manual(values=c(21,22)) +#scale_fill_gradient(low='blue', high='red') + 
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=0.5, linetype=2) + 
  annotate("text",x=NMDS.mean$NMDS1,y=NMDS.mean$NMDS2,label=NMDS.mean$group, size=4) +
  scale_fill_manual(values=c('black', 'dark blue', 'light blue', 'dark green', 'green', 'white', 'grey','yellow', 'orange', 'red')) +
  scale_color_manual(values=c('black', 'dark blue', 'light blue', 'dark green', 'green','white', 'grey', 'yellow', 'orange', 'red')) 


# Plot nmds against function:Fig S4
## import the function data and merge with the metadata table
fun.table=read.table('input_data/FunExp_metadata.txt', header=TRUE)
colnames(fun.table)[3]=  colnames(meta.table2.nmds)[2]
funmeta.table=left_join(fun.table, meta.table2.nmds)
funmeta.table2=funmeta.table[-which(is.na(funmeta.table[,10])),]

divfun_sup_nmdsplo1=ggplot(funmeta.table2, aes(x=NMDS1, y=W_Change)) + 
  geom_point(aes(color=factor(leaf_age_weeks))) + Theme2 + ylab('Weigth loss (g)') +
  scale_color_manual(values=c('black', 'dark blue', 'light blue', 'dark green', 'green', 'white', 'grey','yellow', 'orange', 'red')) 


divfun_sup_nmdsplot2=ggplot(funmeta.table2, aes(x=NMDS2, y=W_Change)) + 
  geom_point(aes(color=factor(leaf_age_weeks))) + Theme2 + ylab('Weigth loss (g)') +
  scale_color_manual(values=c('black', 'dark blue', 'light blue', 'dark green', 'green', 'white', 'grey','yellow', 'orange', 'red')) 

quartz(width=7.5, height=3.5)
grid.arrange(divfun_sup_nmdsplo1, divfun_sup_nmdsplot2, ncol=2)
#quartz.save('/Figures/FigS4_NMDSaxis.pdf')

summary(lm(W_Change ~ NMDS1 * leaf_age_weeks, data=funmeta.table2))
summary(lm(W_Change ~ NMDS2 * leaf_age_weeks, data=funmeta.table2))

# no significant trends - Fig S4

# distance analysis
dists=vegdist(ASVtable2, method='bray')
groups.age = meta.table2$leaf_age
groups.treat = meta.table2$treatment

dists.age=betadisper(dists,  groups.age)
anova(dists.age)

dists.treat=betadisper(dists,  groups.treat)
anova(dists.treat)