# FunExp
Code that supports figures and stats from FunExp project submitted for publication with the title 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves'. 

## Contents and narrative
This story depicts the changes in bacterial diversity and function along the life of a pitcher plant leaf. A central figure in this study depicts the relationship betwen ASV richness and Degradation, showing a negative rela diversity - function relationship. Albeit not causal because we did not manipulate diversity directly, this is one of the first studies showing this negative trend.  [BEF_funexp_main.R](https://github.com/catalicu/FunExp/blob/main/BEF_funexp_main.R) contains the figure and stats to represent this relationship. The stats include model selection as conducted in [Stat_BEFmain_models.R](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main.R).
To better understand how ASV richness and degradation change with leaf age, we built figures representing this diversity in 
[Fig and stats_ FunAge_noOld.R](https://github.com/catalicu/FunExp/blob/main/Fig_and_stats_FunAge_noOld.R). 

**Controls**  
To establish the diversity and function baselines before running the degradation experiment, we built a [Fig_DivBaseline.R](https://github.com/catalicu/FunExp/blob/main/Fig_DivBaseline.R). 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in [qPCR_FunExp.R](https://github.com/catalicu/FunExp/blob/main/qPCR_FunExp.Rmd).  

**Overal projects**   
FunASV_v11.Rmd
FunExp.Rproj


