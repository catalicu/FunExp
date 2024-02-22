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

# edit tables
## table with ASVs exclusively
ASVtable=ASVtable.fundiv[,2:1376]
## metadata exclusive
metatable=ASVtable.fundiv[,1377:1394]


# from v10 (unedited:
### Calculate: mean abundance, occurrence and mean functional output per ASV

#presence absence: occurrence
ASVpa=ASVtable
ASVpa[(ASVpa!=0)]=1
ASVpa=data.frame(ASVpa)
ASVpa=sapply(ASVpa, as.numeric)
ASVpa.freq=colSums(ASVpa)/(dim(ASVpa)[1])

#relative abundance
ASVrel.abund=colSums(ASVtable)

### Combine these values into a table

#### Relative abundance vs occupancy.
#To better visualize this relationship:
 # * create a data frame with the relative abundance and the occupancy data
 # * classify the taxa into 3 groups based on their occurrence and abundance
 # * save this table for later use
# data frame with abundance and freq
ASVrel.abund.pa=data.frame(ASVrel.abund, ASVpa.freq)
# highlight frequency taxa
ASVrel.abund.pa$Frequency=ASVrel.abund.pa$ASVpa.freq
ASVrel.abund.pa[which(ASVrel.abund.pa$ASVpa.freq>0.035),3]='Frequent'
ASVrel.abund.pa[which(ASVrel.abund.pa$ASVpa.freq<0.035),3]='Unfrequent'
ASVrel.abund.pa[which(ASVrel.abund.pa$ASVrel.abund>1000000),3]='Most Frequent'
ASVrel.abund.pa$Frequency=factor(ASVrel.abund.pa$Frequency, levels = c('Unfrequent', 'Frequent', 'Most Frequent'))

#### Add the Functional Ratio to the abundance and occupancy
ASVrel.abund.pa.fun=data.frame(ASVrel.abund.pa) #, ASVfunRat) removed funRat because it was a weird metric
ASVrel.abund.pa.fun$ASVcode=rownames(ASVrel.abund.pa.fun)

# Join relative abundance, occurrence and slope classification
ASVrel.abund.pa.fun.sorted=left_join(ASVrel.abund.pa.fun, ASVslope_sorted)
# label the taxa that were not analyzed for slope as 'undetermined' slopes.
ASVrel.abund.pa.fun.sorted[which(is.na(ASVrel.abund.pa.fun.sorted$Slope)),6]='Undetermined'

# Save this table:
write.table(ASVrel.abund.pa.fun.sorted, file='output_tables/ASVrel.abund.pa.txt', sep='\t')

ASVrel.abund.pa_plot_c=ggplot(ASVrel.abund.pa.fun.sorted, aes(log(ASVrel.abund), ASVpa.freq)) + 
  geom_point(aes(fill=Slope), size=3, shape=21) + Theme +ggtitle('Relationship between ASV occupancy (frequency) and relative abundance') + scale_fill_manual(values=c('black', 'grey', 'white', 'blue')) + ylab('ASV occurrence') + xlab('ASV Relative Abundance')

# Before melting, reduce the dataset to taxa that we were able to relate to function (in scritp v8)
# select the columns that will go into the long format:
# include treatment, leaf age, and taxa that were realted to function
ASVfor_melt=ASVtable.fundiv[,c(which(colnames(ASVtable.fundiv)%in%(ASVslope_sorted$ASVcode)),1377:1394)]
# Then use melt to go from the wide table format to the long table format. 
otumelt_fun.div=melt(ASVfor_melt)
length(names(otumelt_fun.div))
names(otumelt_fun.div)=c('leaf_age_weeks2',"ASVcode", "Abundance") # the last two columns will represent the ASV code an the abundance value 
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
