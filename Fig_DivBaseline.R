# Libraries
library(ggplot2)
library(vegan)
library(plyr)
library(nlme)
library(gridExtra)
library(dplyr)
library(car)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
    theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# set wd
 setwd("/Users/cata/Documents/CGLab/GitHub/FunExp_Rstudio/FunExp")

#Color friendly: The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Load data
ASVtable=read.table('input_data/ASVtable_FunExp12021-03-24.txt', header=TRUE)
meta.table=read.table('input_data/METAtable_FunExp12021-03-24_sampleNamesFIX.txt', header=TRUE)
# reformatting:
meta.table$leaf_age_weeks=as.numeric(substr(meta.table$leaf_age, 2,3))
meta.table[which(is.na(meta.table$leaf_age_weeks)),6]=0 # replace NA with zero
# to calculate leaves per treatment
leavs_per_age=read.csv('input_data/Leavs_per_age.csv', header=TRUE)
# remove leaves with no fluid
leavs_per_age2= leavs_per_age[-which(is.na(leavs_per_age$Volume)),]
leaves_count=ddply(leavs_per_age, ('AGE_WEEKS'), summarize, count=length(BAND))
colnames(leaves_count)[1]='leaf_age_weeks'

# calculate richness metrics
richness=specnumber(ASVtable)
shannon=diversity(ASVtable, index='shannon')
simpson=diversity(ASVtable, index='simpson')
evenness=shannon/log(richness)
meta.div=data.frame(meta.table, richness, shannon, simpson, evenness)
head(meta.div)

#extract data to focus on pre-experiemnt diversity
Baseline = meta.div[which(meta.div$treatment=='F50'),]
# merge with baseline dataframe
Baseline2=left_join(Baseline, leaves_count, by = "leaf_age_weeks")
# create weighted version of richness
Baseline2$w_richness=Baseline2$richness/Baseline2$count

# plot richness panels
rich_regional=ggplot(Baseline, aes(leaf_age_weeks, richness))  + 
  geom_point(size=2)+ Theme + ylab('Richness') + 
  xlab('Leaf age (weeks)') + 
  geom_smooth(method=lm, formula=y ~ (x), se=FALSE, color='black', size=0.5) + 
  annotate('text', label='a.', x=1, y=115)
rich_local=ggplot(Baseline2, aes(leaf_age_weeks, w_richness))  + 
  geom_point(size=2)+ Theme + ylab('Weighted Richness per leaf') + 
  xlab('Leaf age (weeks)') + 
  geom_smooth(method='lm', formula= y ~ (x), se=FALSE, color='black', size=0.5)+ 
  annotate('text', label='b.', x=1, y=4)


Fig_divControl=arrangeGrob(rich_regional, rich_local, ncol=2)
quartz(width=6, height=3)
plot(Fig_divControl)
quartz.save('Figures/Fig_div_time_controls.png', type='png', dpi=300)

# stats
# rich regional
model_baseline=glm(richness ~ leaf_age_weeks, data=Baseline)
summary(model_baseline)
Anova(model_baseline)

# rich_local
model_baseline2a=glm(w_richness ~ leaf_age_weeks, data=Baseline2)
model_baseline2b=glm(w_richness ~ poly(leaf_age_weeks,2), data=Baseline2)
AIC(model_baseline2a, model_baseline2b)
summary(model_baseline2a)
Anova(model_baseline2a) # wald test

### repeat for evenness
# plot evenness panels
eve_regional=ggplot(Baseline, aes(leaf_age_weeks, evenness))  + 
  geom_point(size=2)+ Theme + ylab('Richness') + 
  xlab('Leaf age (weeks)') + 
  geom_smooth(method=lm, formula=y ~ (x), se=FALSE, color='black', size=0.5) + 
  annotate('text', label='a.', x=0.8, y=0.7)


plot(eve_regional)
quartz.save('Figures/Fig_eve_time_controls_supmat.png', type='png', dpi=300)

# stats
# eve regional
model_eve_baseline=glm(evenness ~ leaf_age_weeks, data=Baseline)
summary(model_eve_baseline)
Anova(model_eve_baseline)
