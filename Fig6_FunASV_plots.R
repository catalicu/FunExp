#title: "FunASV_plots"
#author: "Dr CG"
#date:"1/13/2023"
# taken from:

# Description:
# This script creates figures to identify ASV underlying function:
  # * ASV relative abundances vs function
  # * General patterns of ASV distribution and occurrences

# Libraries

# Libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(linewidth = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# Load data
ASVtable.fundiv=read.table('input_data/ASVtable_forGLM_FunExp12023-01-04.txt', header=TRUE)
ASVtable.fundiv1=ASVtable.fundiv[,-1]
taxatable=read.table('input_data/taxtable_FunExp12023-01-04.txt', header=TRUE)
ASVslope_sorted=read.table('output_tables/ASVslope_sorted52024-02-21.txt', header=TRUE)

# list of ASV codes to go through:
ASVfun_names=names(ASVtable.fundiv1)[1:1375]

# Before melting, reduce the dataset to taxa that we were able to relate to function (in scritp v8)
# select the columns that will go into the long format:
# include treatment, leaf age, and taxa that were realted to function
ASVfor_melt=ASVtable.fundiv[,c(which(colnames(ASVtable.fundiv)%in%(ASVslope_sorted$ASVcode)),1377:1394)]
# Then use melt to go from the wide table format to the long table format. 
otumelt_fun.div=melt(ASVfor_melt, id=c("richness"))
length(names(otumelt_fun.div))
names(otumelt_fun.div)=c('richness',"ASVcode", "Abundance") # the last two columns will represent the ASV code an the abundance value 
head(otumelt_fun.div)

#Clean the tax table to add tax information to the long version of the ASV table with metadata (otumelt)
otumelt2.div=left_join(otumelt_fun.div, taxatable,  by='ASVcode')

# make sure the relative abundance column is numeric
otumelt2.div$Abundance=as.numeric(as.character(otumelt2.div$Abundance))


#### Generate list of ASVs with taxonomic and abundance data
#Generate a list of unique family names in order of decreasing abundance	 (object: familylist3)
#Also create a table with abundance per family (object: meanlist.fam3).

ASVlist.div=unique(otumelt2.div$ASVcode)	#866 taxa
# calculate means for abundance per taxon
ASVlist_means.div=ddply(otumelt2.div, .(ASVcode), summarize,  Abund=mean(Abundance, na.rm=TRUE))
ASVlist_means.div=ASVlist_means.div[order(ASVlist_means.div$Abund, decreasing=TRUE),]
ASVlist2.div=left_join(ASVlist_means.div, taxatable)

#### Do the most abundant taxa increase over time?
otumelt2_ASV1.div=otumelt2.div[which(otumelt2.div$ASVcode=='ASV1'),]
ggplot(otumelt2_ASV1.div, aes(richness, (Abundance+1))) + geom_point() + geom_smooth(method='lm', se=FALSE, color='black')

#### Are there general trends of taxa over time?
ggplot(otumelt2.div, aes(richness, log(Abundance+1))) + geom_jitter(aes(color=ASVcode), alpha=0.4)  + geom_smooth(method='lm', se=FALSE, color='black') + theme(legend.position='none') + geom_smooth(method='lm', aes(color=ASVcode), se=FALSE)
