# FunExp
Code that supports figures and stats from FunExp project.

## Contents and narrative
The central figure in this study depicts the relationship betwen ASV richness and Degradation, as they demonstrate a negative diversity - function relationship. [BEF_funexp_main.R](https://github.com/catalicu/FunExp/blob/main/BEF_funexp_main.R) contains the figure and stats to represent this relationship. The stats include model selection as conducted in Stat_BEFmain_models.R.
To better understand how ASV richness and degradation change with leaf age, we built figures representing this diversity in 
Fig_AGEvs_aDiv_bFun.R and Stat_AGEvs_aDiv_bFun.R. We repeated this analysis while removing the oldest leaves to assess their role in 
Fig_time_fun_noold.R in stats_FunAge_noOld. 

**Controls**  
To establish the diversity and function baselines before running the degradation experiment, we built a Fig_DivBaseline.R. 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in qPCR_FunExp.R.  

**Overal projects**   
FunASV_v11.Rmd
FunExp.Rproj


