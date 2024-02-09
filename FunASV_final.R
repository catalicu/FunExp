#title: "FunExp_v11"
#author: "Dr CG"
#date:"1/13/2023"

# Libraries
library(ggplot2)
#library(vegan)
library(dplyr)
#library(plyr)
#library(lme4)
#library(nlme)
#library(gridExtra)
library(reshape2)
#library(car)
#library(MASS)

# Plot themes
## With legend
Theme=theme_classic(base_size=11, base_family="Helvetica") +
  theme(axis.line = element_line(size = 1, colour = "black", linetype = "solid")) +theme(plot.title = element_text(size = 12))
## Without legends
Theme2=Theme+ theme(legend.position="none") + theme(panel.border=element_rect(fill=NA))

# Load data
ASVtable.fundiv=read.table('input_data/ASVtable_forGLM_FunExp12023-01-04.txt', header=TRUE)
ASVtable.fundiv1=ASVtable.fundiv[,-1]
taxatable=read.table('input_data/taxtable_FunExp12023-01-04.txt', header=TRUE)

# list of ASV codes to go through:
ASVfun_names=names(ASVtable.fundiv1)[1:1375]

# from FunASV_v8: create an iterative GLM function
## ### To Iterate the GLM:
# * glm (gaussian distribution) vs null model
# * Store the model's outputs and the comparison with the null. 

#Input:
# * ASVfun = wide format dataset with ASV, fundiv data: ASVtable.fundiv1
# * taxatable = list of taxa with number, DNA sequence, taxonomic affiliation and organized from more abundant to least abundant: taxatable
# * minocc = minimum number of occurrences for taxa to be run through the model
IterativeGLM_Fun <- function(ASVfun, taxatable, minocc){
  ASVfun_names=names(ASVfun)[1:1375]
  # create objects to store model results
  model.glmFun.results=c() # store model coefficients and output
  model.glmFun.evaluation=c() # store comparison between models
  
  for (i in 1:length(ASVfun_names)) {
    ASVdataFun=(ASVfun[,c(i,1376:1393)])
    names(ASVdataFun)[1]='ASVrelabun'
    
    # skip columns with 0 to 4 appearances
    ASVdataFun_presence=ASVdataFun[1]
    ASVdataFun_presence[which(ASVdataFun_presence>1),]=1
    if (sum(ASVdataFun_presence)<minocc) {print(paste('low occurrence (<)',minocc))} else {
      
      tax_rownum=which(taxatable$ASVcode==ASVfun_names[i])
      
      # data clean up steps   
      if (length(which(ASVdataFun$ASVrelabun==0))>0) {
        ASVdataFun=ASVdataFun[-which(ASVdataFun$ASVrelabun==0),] } else {print('nope')}
      if (length(which(is.na(ASVdataFun$W_Change)))>0) {
        ASVdataFun=ASVdataFun[-which(is.na(ASVdataFun$W_Change)),] } else {print('nope2')}
      
      # Check progress
      print(paste('Progress:', (ASVfun_names[i]), 
                  (taxatable$Family[tax_rownum]),
                  (taxatable$Genus[tax_rownum]), sep=' '))
      
      # Model runs
      # GLM with gaussian distribution
      Mod.fun1=lm(W_Change ~ ASVrelabun, data=ASVdataFun)
      Mod1.summary=summary(Mod.fun1)
      Mod1.coeff=Mod1.summary$coefficients
      # Null model
      Mod.null=glm(W_Change ~ 1, data=ASVdataFun)
      
      # Model comparison
      Mod1.aov=anova(Mod.fun1, Mod.null, test='Chisq')
      aic.model=AIC(Mod.fun1,Mod.null)		
      
      # Store coefficients
      coeff.full.model = data.frame((Mod1.summary$coefficients), 
                                    adj.p=p.adjust(Mod1.summary$coefficients[,2]), 
                                    Family=rep(taxatable$Family[tax_rownum],2),
                                    Genus=rep(taxatable$Genus[tax_rownum],2),
                                    ASVcode=rep(ASVfun_names[i],2), 
                                    ASVseq=rep(taxatable$ASVseq[tax_rownum], 2)) 
      # Store comparison
      coeff.aov = data.frame(model=c(paste(ASVfun_names[tax_rownum], 'null.model', sep='_'),
                                     paste(ASVfun_names[tax_rownum], 'full.model', sep='_')), 
                             resid.df=Mod1.aov[1], resid.dev=Mod1.aov[2], deviance=Mod1.aov[[4]], 
                             AIC=aic.model$AIC, chisq.p=Mod1.aov[[5]], adj.p=p.adjust(Mod1.aov[[5]]), 
                             Family=rep(taxatable$Family[tax_rownum],2) ,
                             Genus=rep(taxatable$Genus[tax_rownum],2), 
                             ASVcode=rep(ASVfun_names[i],2),
                             ASVseq=rep(taxatable$ASVseq[tax_rownum],2))		
      
      # format results table
      model.glmFun.results=(rbind(model.glmFun.results, coeff.full.model)) 		
      model.glmFun.evaluation=(rbind(model.glmFun.evaluation, coeff.aov))		
    }
  }
  # Finalize output:
  dim(model.glmFun.results) # check for expected dimensions
  # Add names to the columns of output tables:
  names(model.glmFun.results)=c('Estimate', 'Serror', 't.value', 'p.value', 'p.adjust','Family', 'Genus', 'ASVcode','ASVseq')
  names(model.glmFun.evaluation)=c('Model', 'df', 'Resid.dev', 'Deviance',  'AIC', 'p.value','adjust.p', 'Family', 'Genus', 'ASVcode','ASVseq')
  
  # Store both objects within a list
  model.results.output.glmFun=list(coefficients=model.glmFun.results, evaluation=model.glmFun.evaluation)
  
  return(model.results.output.glmFun)
}

# run the function by keeping:
  # taxa in 0.14% of samples - 2 occurrences or more - does not run well, low abundances
#model.results.trial1=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,2)
  # taxa in 0.35% of samples - 5 occurrences or more
model.results.trial5=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,5)
  # taxa in 0.71% of samples - 10 occurrences or more
model.results.trial10=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,10)
  # taxa in 1% of samples - 15 occurrences or more
model.results.trial15=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,15)



# sort the output tables
# Input:
# * ASVfun_names = ASV names (codes), the function will go one by one
# * models.coefficients = the coefficient output of iterativeGLM function
SortingASVtables=function(ASVfun_names, models.coefficients){
  ## create lists to store potential ASVs correlated with function
    list_asv_slopePos.Fun=c()
    list_asv_noslope.Fun=c()
    list_asv_slopeNeg.Fun=c()

for (i in 1:length(ASVfun_names)){
  # check that the ASV is present in the dataset
  if (sum(which(models.coefficients$ASVcode==ASVfun_names[i]))==0) { print(paste(ASVfun_names[i], 'ASV not found'))
  } else {
    
    # extract the slope estimate:
    ASVabund.estimate= models.coefficients[which(models.coefficients$ASVcode==ASVfun_names[i]),1][2]
    # extract the p value:
    ASVabund.p= models.coefficients[which(models.coefficients$ASVcode==ASVfun_names[i]),5][2]
    ASVcode= models.coefficients[which(models.coefficients$ASVcode==ASVfun_names[i]),8][2]
    ASVGenus= models.coefficients[which(models.coefficients$ASVcode==ASVfun_names[i]),6][2]
    
    # inform of progress:  
    print(paste('Progress:', i, ASVcode, ASVGenus))
    
    # check that output differs from NA
    if (is.na(ASVabund.p)) { 
      print(paste(ASVfun_names[i], 'p value was NA'))
      list_asv_noslope.Fun=rbind(list_asv_noslope.Fun, c(ASVfun_names[i], ASVabund.estimate, ASVabund.p))
    } else {
      
      # Sort the taxa
      # if the taxa has a positive slope estimate and a significant adjusted p value -> accept and store ASVcode
      if (ASVabund.estimate > 0 && ASVabund.p < 0.05 ) {
        list_asv_slopePos.Fun=rbind(list_asv_slopePos.Fun, c(ASVfun_names[i], ASVabund.estimate, ASVabund.p))
      } else   { 
        # if the taxa has a negative slope estimate and a significant adjusted p value -> accept and store ASVcode
        if (ASVabund.estimate < 0 && ASVabund.p < 0.05 ) {
          list_asv_slopeNeg.Fun=rbind(list_asv_slopeNeg.Fun, c(ASVfun_names[i], ASVabund.estimate, ASVabund.p))
        } else {
          # store non significant and NA results ASV in noslope list
          list_asv_noslope.Fun=c(list_asv_noslope.Fun, paste(ASVfun_names[i], ASVabund.estimate, ASVabund.p)) 
        }
      } 
    }
  }
}

#### Edit the sorted tables and save
#Change from matrix to data frame and add column names.
#Focus on negative and positive slopes from here on. 
list_asv_slopePos.Fun.df=data.frame(list_asv_slopePos.Fun)
  colnames(list_asv_slopePos.Fun.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listPos.rank=list_asv_slopePos.Fun.df[order(list_asv_slopePos.Fun.df$ASVestimate),]

list_asv_slopeNeg.Fun.df=data.frame(list_asv_slopeNeg.Fun)
  colnames(list_asv_slopeNeg.Fun.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listNeg.rank=list_asv_slopeNeg.Fun.df[order(list_asv_slopeNeg.Fun.df$ASVestimate),]

list_asv_noslope.Fun.df=data.frame(list_asv_noslope.Fun)
  colnames(list_asv_noslope.Fun.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listnoslope.rank=list_asv_noslope.Fun.df[order(list_asv_noslope.Fun.df$ASVestimate),]
  
# Merge the lists with their taxonomic affiliations
  listPos.rank_taxa=left_join(listPos.rank, taxatable)
  listnoslope.rank_taxa=left_join(listnoslope.rank, taxatable)
  listNeg.rank_taxa=left_join(listNeg.rank, taxatable)
  
# Merge all lists
  ASVslope_sorted=rbind(
    data.frame(Slope=rep('Positive',dim(listPos.rank_taxa)[1]), listPos.rank_taxa),
    data.frame(Slope=rep('Neutral',dim(listnoslope.rank_taxa)[1]), listnoslope.rank_taxa),
    data.frame(Slope=rep('Negative',dim(listNeg.rank_taxa)[1]), listNeg.rank_taxa))
  return(ASVslope_sorted)
}

ASVslope_sorted5=SortingASVtables(ASVfun_names, model.results.trial5$coefficients)
ASVslope_sorted10=SortingASVtables(ASVfun_names, model.results.trial10$coefficients)
ASVslope_sorted15=SortingASVtables(ASVfun_names, model.results.trial15$coefficients)

######
# prepare dataset
  # select the columns that will go into the long format:
  # include treatment, leaf age, and taxa that were realted to function
  ASVfor_melt=ASVtable.fundiv[,c(which(colnames(ASVtable.fundiv)%in%(ASVslope_sorted$ASVcode)),1387,1389)]

# Then use melt to go from the wide table format to the long table format. 
  otumelt_fun=melt(ASVfor_melt, id=c("treatment", "leaf_age_weeks"))
  length(names(otumelt_fun))
  names(otumelt_fun)[(length(names(otumelt_fun))-1):length(names(otumelt_fun))]=c("ASVcode", "Abundance") # the last two columns will represent the ASV code an the abundance value 
  head(otumelt_fun)
  
  otumelt2=left_join(otumelt_fun, taxatable,  by='ASVcode')
  
  # make sure the relative abundance column is numeric
  otumelt2$Abundance=as.numeric(as.character(otumelt2$Abundance))
  
  # remove rows with no information in the treatment and leaf age (where did these come from?)
  otumelt3=otumelt2[-which(is.na(otumelt2$treatment)),]
  head(otumelt3)

#### Generate list of ASVs with taxonomic and abundance data
  ASVlist=unique(otumelt2$ASVcode)	
  length(ASVlist)
  # calculate means for abundance per taxon
  ASVlist_means=ddply(otumelt2, .(ASVcode), summarize,  Abund=mean(Abundance, na.rm=TRUE))
  ASVlist_means=ASVlist_means[order(ASVlist_means$Abund, decreasing=TRUE),]
  ASVlist2=left_join(ASVlist_means, taxatable)
  
#### research questions:
  #### Do the most functional taxa increase over time?
  #ASV691 was identified as a low abundance taxon with the highest performance and classified as Rodobacter sp. 
  otumelt3_ASV691=otumelt3[which(otumelt3$ASVcode=='ASV691'),]
  otumelt3_ASV691_w0=otumelt3_ASV691[-which(otumelt3_ASV691$Abundance==0),]
  ggplot(otumelt3_ASV691_w0, aes(leaf_age_weeks, (Abundance+1))) + geom_point() + 
  ylim(0,20) +xlim(0, 24) +  Theme + 
    xlab('Leaf age (weeks)') + ylab('log(Abundance +1)')
  
  #### Do the most abundant taxa increase over time?
  #ASV1 Enterobacteraceae - It was identified as having a significant negative relationship with function.
  otumelt3_ASV1=otumelt3[which(otumelt3$ASVcode=='ASV1'),]
  otumelt3_ASV1_w0=otumelt3_ASV1#[-which(otumelt3_ASV1$Abundance==0),]
  ggplot(otumelt3_ASV1_w0, aes(leaf_age_weeks, (Abundance+1))) + geom_point() +xlim(0, 24) + 
    geom_smooth(method='lm', se=FALSE, color='black') + Theme + 
    xlab('Leaf age (weeks)') + ylab('log(Abundance +1)')
  
  #### Are there general trends of taxa over time?
  #Plot them all together to see general trends. this suggests:
  #  * most abundant ones decrease a few remain constant. 
  #  * some of the less aboundant increase with time
  ggplot(otumelt3, aes(leaf_age_weeks, log(Abundance+1))) + 
    geom_jitter(aes(color=ASVcode), alpha=0.4) +xlim(0, 24) + 
    geom_smooth(method='lm', se=FALSE, color='black') + 
    theme(legend.position='none') + geom_smooth(method='lm', aes(color=ASVcode), se=FALSE) + 
    Theme + theme(legend.position='none')
  
  ### Iterative GLM model
 # Now build an iterative model function:
  #  * glm (poisson distribution) vs null model
  # * Store the model's outputs and the comparison with the null. 
# Input:
#* long format dataset with ASV and metadata: otumel2
#* list of taxa with number, DNA sequence, taxonomic affiliation and organized from more abundant to least abundant: ASVlist2
  IterativeGLM1 <- function(long_data, ASVlist2){
    model.results.otuA=c() # store model coefficients and output
    model.evaluation.otuA=c() # store comparison between models
    
    for (i in 1:length(ASVlist2$ASVcode)){
      # extract data
      OTUdata=otumelt2[which(otumelt2$ASVcode%in%ASVlist2$ASVcode[i]),] 
      OTUdata2=OTUdata 
      # GLM with Poisson distribution
      Mod1=glm.nb(Abundance ~ leaf_age_weeks, data=OTUdata2)
      Mod1.summary=summary(Mod1)
      Mod1.coeff=Mod1.summary$coefficients
      # Null model
      Mod1.null=glm.nb(Abundance ~ 1, data=OTUdata2)
      Mod1.aov=anova(Mod1, Mod1.null, test='Chisq')
      # Comparison	
      aic.model=AIC(Mod1,Mod1.null)
      # store coefficients
      coeff.full.model = data.frame((Mod1.summary$coefficients), 
                                    adj.p=p.adjust(Mod1.summary$coefficients[,2]), 
                                    ASVcode=rep(unique(OTUdata2$ASVcode), 2), 
                                    Family=rep(OTUdata2$Family[i],2) ,
                                    Genus=rep(OTUdata2$Genus[i],2), 
                                    ASVcode=rep(ASVlist2$ASVcode[i],2)) 
      # store comparison
      
      coeff.aov = data.frame(model=c(paste(unique(OTUdata2$ASVcode), 'null.model', sep='_'),
                                     paste(unique(OTUdata2$ASVcode), 'full.model', sep='_')), 
                             resid.df=Mod1.aov[3], 
                             resid.theta=Mod1.aov[2], 
                             #deviance=Mod1.aov[[4]], 
                             AIC=aic.model$AIC, 
                             p.value=c('NA', (Mod1.aov$`Pr(Chi)`)[2]), 
                             adj.p=c('NA', p.adjust(Mod1.aov$`Pr(Chi)`)[2]), 
                             Family=rep(OTUdata2$Family[i],2) ,
                             Genus=rep(OTUdata2$Genus[i],2), 
                             ASVcode=rep(OTUdata2$ASVcode[i],2))		
      
      # format result table
      model.results.otuA=(rbind(model.results.otuA, coeff.full.model)) 			
      model.evaluation.otuA=(rbind(model.evaluation.otuA, coeff.aov))					
      
      # print taxa number, family and genus as a way to track progress of the function: remove if >100 taxa
      print(paste(OTUdata$Family[i], OTUdata$Genus[i], i))
    }
    dim(model.results.otuA) # check the dimensions
    
    # add names to the columns of the output tables
    names(model.results.otuA)=c('Estimate', 'Sdev', 'z.value', 'p.value', 'p.adjust', 'ASV2','Family', 'Genus', 'ASVcode')
    names(model.evaluation.otuA)=c('Model', 'Resid.dev', 'theta', 'AIC', 'p.value','adjust.p', 'Family', 'Genus', 'ASVcode')
    
    # store both data frames within a list
    model.results.output = list(coefficients=model.results.otuA, evaluation=model.evaluation.otuA)
    
    return(model.results.output)
  } 
  # run it:
  model.results.trial1=IterativeGLM1(otumelt2, ASVlist2[1:30,])
  models.coefficients1 = model.results.trial1$coefficients
  models.evals1 = model.results.trial1$evaluation  
  
  
#### research questions:
  ### ID taxa that increase with time from Iterative model
  #Can I identify which taxa are increasing with age relably to later test whether they are correlated with function?
  # create lists to store potential ASVs correlated with function
  list_asv_slopePos.time=c()
  list_asv_noslope.time=c()
  list_asv_slopeNeg.time=c()
  
  # list of ASV codes to go through:
  ASVtime_names=ASVlist2$ASVcode
  
  for (i in 1:length(ASVtime_names)){
    # check that the ASV is present in the dataset
    if (sum(which(models.coefficients1$ASVcode==ASVtime_names[i]))==0) { print(paste(ASVtime_names[i], 'ASV not found'))
    } else {
      # extract the slope estimate:
      ASVabund.estimate= models.coefficients1[which(models.coefficients1$ASVcode==ASVtime_names[i]),1][2]
      # extract the p value:
      ASVabund.p= models.coefficients1[which(models.coefficients1$ASVcode==ASVtime_names[i]),5][2]
      ASVcode= models.coefficients1[which(models.coefficients1$ASVcode==ASVtime_names[i]),8][2]
      ASVGenus= models.coefficients1[which(models.coefficients1$ASVcode==ASVtime_names[i]),6][2]
      
      # inform of progress:  
      print(paste('Progress:', i, ASVcode, ASVGenus))
      
      # check that output differs from NA
      if (is.na(ASVabund.p)) { 
        print(paste(ASVtime_names[i], 'p value was NA'))
        list_asv_noslope.time=rbind(list_asv_noslope.time, c(ASVtime_names[i], ASVabund.estimate, ASVabund.p))
      } else {
        
        # Sort the taxa
        # if the taxa has a positive slope estimate and a significant adjusted p value -> accept and store ASVcode
        if (ASVabund.estimate > 0 && ASVabund.p < 0.05 ) {
          list_asv_slopePos.time=rbind(list_asv_slopePos.time, c(ASVtime_names[i], ASVabund.estimate, ASVabund.p))
        } else   { 
          # if the taxa has a negative slope estimate and a significant adjusted p value -> accept and store ASVcode
          if (ASVabund.estimate < 0 && ASVabund.p < 0.05 ) {
            list_asv_slopeNeg.time=rbind(list_asv_slopeNeg.time, c(ASVtime_names[i], ASVabund.estimate, ASVabund.p))
          } else {
            # store non significant and NA results ASV in noslope list
            list_asv_noslope.time=rbind(list_asv_noslope.time, c(ASVtime_names[i], ASVabund.estimate, ASVabund.p)) 
          }
        } 
      }
    }
  }
  
  #### Edit the sorted tables and save
  list_asv_slopePos.time.df=data.frame(list_asv_slopePos.time)
  colnames(list_asv_slopePos.time.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listPos.rank.time=list_asv_slopePos.time.df[order(list_asv_slopePos.time.df$ASVestimate),]
  
  list_asv_slopeNeg.time.df=data.frame(list_asv_slopeNeg.time)
  colnames(list_asv_slopeNeg.time.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listNegslope.rank.time=list_asv_slopeNeg.time.df[order(list_asv_slopeNeg.time.df$ASVestimate),]
  
  list_asv_noslope.time.df=data.frame(list_asv_noslope.time)
  colnames(list_asv_noslope.time.df)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  listnoslope.rank.time=list_asv_noslope.time.df[order(list_asv_noslope.time.df$ASVestimate),]
  
  # merge with taxa
  listPos.rank_taxa.time=left_join(listPos.rank.time, taxatable)
  listnoslope.rank_taxa.time=left_join(listnoslope.rank.time, taxatable)
  listNeg.rank_taxa.time=left_join(listNegslope.rank.time, taxatable)
  
  
  ASVslope.sorted.time=rbind(
    data.frame(Slope.time=rep('Positive',dim(listPos.rank_taxa.time)[1]), listPos.rank_taxa.time),
    data.frame(Slope.time=rep('Neutral',dim(listnoslope.rank_taxa.time)[1]), listnoslope.rank_taxa.time),
    data.frame(Slope.time=rep('Negative',dim(listNeg.rank_taxa.time)[1]), listNeg.rank_taxa.time))
  
  
  ### Compare ASV listsnames(ASVslope.sorted.time)[3:4]=c('ASVestimate.time', 'ASVp.value.time')
  Compare_asv_list=left_join(ASVslope_sorted, ASVslope.sorted.time)
  Compare_asv_list$ASVestimate.time=as.numeric(Compare_asv_list$ASVestimate.time)
  Compare_asv_list$ASVp.value.time=as.numeric(Compare_asv_list$ASVp.value.time)
  
  