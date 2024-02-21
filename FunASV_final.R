#title: "FunASV_final.R"
#author: "Dr CG"
#date:"1/13/2023"
# taken from 'FunExp_v11.R'

# Description:
# This script runs iterative GLMs on individual ASVs, performs correction
# for multiple testing and prints the output in tables that identify 
# significant slopes and their direction. 
# It is set up to run two models:
  # * one with Function (weight loss) ~ ASV abundance 
  # * one with Function (weight loss) ~ ASV abundance * leaf age
# *******continue work at Sorting function line 301


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

# list of ASV codes to go through:
ASVfun_names=names(ASVtable.fundiv1)[1:1375]

# from FunASV_v8: create an iterative GLM function for:
# Model: Weight loss ~ ASV relative abundance
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

# run the IterativeGLM_Fun function by keeping:
  # taxa in 0.14% of samples - 2 occurrences or more - does not run well, low abundances
#model.results.trial1=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,2)
  # taxa in 0.35% of samples - 5 occurrences or more
model.results.trial5=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,5)
  # taxa in 0.71% of samples - 10 occurrences or more
model.results.trial10=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,10)
  # taxa in 1% of samples - 15 occurrences or more
model.results.trial15=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,15)

# from FunASV_v13: create an iterative GLM function for:
# Model: Weight loss ~ ASV relative abundance * time
## ### To Iterate the GLM:
# * glm (gaussian distribution) vs null model
# * Store the model's outputs and the comparison with the null. 
#Input:
# * ASVfun = wide format dataset with ASV, fundiv data: ASVtable.fundiv1
# * taxatable = list of taxa with number, DNA sequence, taxonomic affiliation and organized from more abundant to least abundant: taxatable
# * minocc = minimum number of occurrences for taxa to be run through the model

IterativeGLM_Fun_time <- function(ASVfun, taxatable, minocc){
  ASVfun_names=names(ASVfun)[1:1375]
  # create objects to store model results
  model.glmFun.results=c() # store model coefficients and output
  model.glmFun.evaluation=c() # store comparison between models
  
  for (i in 1:length(ASVfun_names)) {
    ASVdataFun=(ASVfun[,c(i,1376:1393)])
    names(ASVdataFun)[1]='ASVrelabun'
    
    tax_rownum=which(taxatable$ASVcode==ASVfun_names[i])
    
    # Check progress
    print(paste('Progress:', (ASVfun_names[i]), 
                (taxatable$Family[tax_rownum]),
                (taxatable$Genus[tax_rownum]), sep=' '))
    
    # data clean up steps   
    ## remove rows with zero abundance
    if (length(which(ASVdataFun$ASVrelabun==0))>0) {
      ASVdataFun=ASVdataFun[-which(ASVdataFun$ASVrelabun==0),] } else {print('nope')}
    ## remove rows with NA in the function 
    if (length(which(is.na(ASVdataFun$W_Change)))>0) {
      ASVdataFun=ASVdataFun[-which(is.na(ASVdataFun$W_Change)),] } else {print('nope2')}
    
    # skip columns with 0 to 4 appearances
    ASVdataFun_presence=ASVdataFun[1]
    ASVdataFun_presence[which(ASVdataFun_presence>1),]=1
    if (sum(ASVdataFun_presence)<minocc) {print(paste('low occurrence (<',minocc))} else {
      
      # Model runs
      # GLM with gaussian distribution
      Mod.fun1=lm(W_Change ~ ASVrelabun + leaf_age_weeks, data=ASVdataFun)
      Mod1.summary=summary(Mod.fun1)
      Mod1.coeff=Mod1.summary$coefficients
      # Null model
      Mod.null=glm(W_Change ~ 1, data=ASVdataFun)
      
      # Model comparison
      Mod1.aov=anova(Mod.fun1, Mod.null, test='Chisq')
      aic.model=AIC(Mod.fun1,Mod.null)		
      
      # Store coefficients
      coeff.full.model = data.frame((Mod1.summary$coefficients), 
                                    adj.p=p.adjust(Mod1.summary$coefficients[,4]),
                                    Family=rep(taxatable$Family[tax_rownum],3),
                                    Genus=rep(taxatable$Genus[tax_rownum],3),
                                    ASVcode=rep(ASVfun_names[i],3), 
                                    ASVseq=rep(taxatable$ASVseq[tax_rownum],3)) 
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

# run the IterativeGLM_Fun_time function by keeping:
# taxa in 0.14% of samples - 2 occurrences or more - does not run well, low abundances
#model.results.trial1=IterativeGLM_Fun(ASVtable.fundiv1, taxatable,2)
# taxa in 0.35% of samples - 5 occurrences or more
model.results_time.trial5=IterativeGLM_Fun_time(ASVtable.fundiv1, taxatable,5)
# taxa in 0.71% of samples - 10 occurrences or more
model.results_time.trial10=IterativeGLM_Fun_time(ASVtable.fundiv1, taxatable,10)
# taxa in 1% of samples - 15 occurrences or more
model.results_time.trial15=IterativeGLM_Fun_time(ASVtable.fundiv1, taxatable,15)

# sort the output tables
# this function identifies positive and negative significant slopes
# Input:
# * ASVfun_names = ASV names (codes), the function will go one by one
# * models.coefficients = the coefficient output of iterativeGLM function
models.coefficients=model.results_time.trial15$coefficients

SortingASVtables=function(ASVfun_names, models.coefficients){
  ## create lists to store potential ASVs correlated with function
    list_asv_slopePos.Fun=c()
    list_asv_noslope.Fun=c()
    list_asv_slopeNeg.Fun=c()

#i=1
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
          list_asv_noslope.Fun=rbind(list_asv_noslope.Fun, c(ASVfun_names[i], ASVabund.estimate, ASVabund.p)) 
        }
      } 
    }
  }
}

#### Edit the sorted tables:
#Change from matrix to data frame and add column names.
list_asv_slopePos.Fun.df=data.frame(list_asv_slopePos.Fun)
if (nrow(list_asv_slopePos.Fun.df)>0) {
  listPos.rank=list_asv_slopePos.Fun.df[order(list_asv_slopePos.Fun.df[,2]),]
                                    } else {
                                      listPos.rank=rbind(c('A000', rep(1, 2)),c('A000', rep(1, 2)))
                                    }
colnames(listPos.rank)=c('ASVcode', 'ASVestimate', 'ASVp.value')

list_asv_slopeNeg.Fun.df=data.frame(list_asv_slopeNeg.Fun)
if (nrow(list_asv_slopeNeg.Fun.df)>0) {  
  listNeg.rank=list_asv_slopeNeg.Fun.df[order(list_asv_slopeNeg.Fun.df[,2]),]
                                    } else {
                                      listNeg.rank=rbind(c('A000', rep(1, 2)),c('A000', rep(1, 2)))
                                    }
colnames(listNeg.rank)=c('ASVcode', 'ASVestimate', 'ASVp.value')

list_asv_noslope.Fun.df=data.frame(list_asv_noslope.Fun)
if (nrow(list_asv_noslope.Fun.df)>0) {
  listnoslope.rank=list_asv_noslope.Fun.df[order(list_asv_noslope.Fun.df[,2]),]
                                    } else {
                                      listnoslope.rank=rbind(c('A000', rep(1, 2)),c('A000', rep(1, 2)))
                                    }
colnames(listnoslope.rank)=c('ASVcode', 'ASVestimate', 'ASVp.value')
  
# Merge the lists with their taxonomic affiliations
if (sum(as.numeric(listPos.rank[,2]))==2) { listPos.rank=data.frame(listPos.rank)
      } else { 
        listPos.rank=listPos.rank  
      }
  listPos.rank_taxa=left_join(listPos.rank, taxatable)

if (sum(as.numeric(listnoslope.rank[,2]))==2) { listnoslope.rank=data.frame(listnoslope.rank)
      } else {
        listnoslope.rank=listnoslope.rank
      }
  listnoslope.rank_taxa=left_join(listnoslope.rank, taxatable) 
  
if (sum(as.numeric(listNeg.rank[,2]))==2) { listNeg.rank=data.frame(listNeg.rank)
      } else {
  listNeg.rank=listNeg.rank
      }
  listNeg.rank_taxa=left_join(listNeg.rank, taxatable)
  
# Merge all lists
  ASVslope_sorted=rbind(
    data.frame(Slope=rep('Positive',dim(listPos.rank_taxa)[1]), listPos.rank_taxa),
    data.frame(Slope=rep('Neutral',dim(listnoslope.rank_taxa)[1]), listnoslope.rank_taxa),
    data.frame(Slope=rep('Negative',dim(listNeg.rank_taxa)[1]), listNeg.rank_taxa))
  return(ASVslope_sorted)
}

# from models: weight loss ~ asv abundance
ASVslope_sorted5=SortingASVtables(ASVfun_names, model.results.trial5$coefficients)
ASVslope_sorted10=SortingASVtables(ASVfun_names, model.results.trial10$coefficients)
ASVslope_sorted15=SortingASVtables(ASVfun_names, model.results.trial15$coefficients)

# from models: weight loss ~ asv abundance * time
ASVslope_sorted_time5=SortingASVtables(ASVfun_names, model.results_time.trial5$coefficients)
ASVslope_sorted_time10=SortingASVtables(ASVfun_names, model.results_time.trial10$coefficients)
ASVslope_sorted_time15=SortingASVtables(ASVfun_names, model.results_time.trial15$coefficients)

# summary table:
data.frame(
  model=c('basic','basic','basic', 'time', 'time','time'),
  occurrence=c(5, 10, 15, 5, 10, 15),  
  NoASV=c(dim(model.results.trial5$evaluation)[1]/2, 
    dim(model.results.trial10$evaluation)[1]/2, 
    dim(model.results.trial15$evaluation)[1]/2, 
    dim(model.results_time.trial5$evaluation)[1]/2, 
    dim(model.results_time.trial10$evaluation)[1]/2, 
    dim(model.results_time.trial15$evaluation)[1]/2))
