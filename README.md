# FunExp
Code that supports figures and stats from FunExp project submitted for publication with the title 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves'. 

## Contents and narrative
This story depicts the changes in bacterial diversity and function along the life of a pitcher plant leaf. 
A central figure in this study depicts the relationship betwen ASV richness and Degradation, showing a negative relationship between diversity and in pitcher plant bacterial communities. Albeit not causal because we did not manipulate diversity directly, this is one of the first studies showing this negative trend due to diversity decreasing and function increasing with leaf age. 

[Fig1_BEF_funexp_main.R](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main.R) contains the figure and stats to represent this general diversity-funciton relationship. The stats include model selection as conducted in [Stat_BEFmain_models.R](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main.R).

Other files and figures reflect the remaining figures and statistics included in the manuscript. 

**Controls**  
To establish the diversity and function baselines before running the degradation experiment, we built a [Fig_DivBaseline.R](https://github.com/catalicu/FunExp/blob/main/Fig_DivBaseline.R). 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in [qPCR_FunExp.R](https://github.com/catalicu/FunExp/blob/main/qPCR_FunExp.Rmd).  


