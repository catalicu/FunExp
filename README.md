This README.txt file was updated on 04/16/2025 
A. Paper associated with this archive
Title: 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves'
Citation: To be updated upon acceptance
Brief abstract: 
How biodiversity and ecosystem functions change with succession has proven to be difficult to predict. We hypothesize that community diversity and function may respond in opposite ways to successional drivers such as nutrient availability, species interactions or abiotic stress. The microbial communities within Sarracenia purpurea leaves perform degradation functions, providing essential nutrients to the plant, but we know little about how succession within the leaf influences bacterial diversity and degradation. We collected pitcher plant fluid from leaves aged 2 to 24 weeks to use in microcosm experiments. We added a common bacterivore to half of the replicated microcosms to establish whether predation is a central successional driver. We used amplicon sequencing and a degradation assay to quantify diversity and ecosystem function.
Names, institutions of all authors
C. Contact information
Name: To be updated upon acceptance
Address:To be updated upon acceptance
email:To be updated upon acceptance
D. Dates of data collection
To be updated upon acceptance
E. Geographic Location(s) of data collection
To be updated upon acceptance
F. Funding Sources
To be updated upon acceptance
OMIT this General Information (above) for double-blind review but include it on final acceptance

DATA & CODE FILE OVERVIEW
This data repository consist of data files, code scripts, and this README document, with the following data and code filenames and variables
Data files are stored in the input_data folder
    1. TAXtable Taxonomy table
    2. ASVtable Species list
    3. METAtable Metadata1: Leaf age code, Treatment code, Replicate code.
    4. FundExp_metadata: Metadata2, contains weight loss data: Tube_name, Treatment age code, Treatment code, Age code, Protozoan treatment, replicate code, initial worm weight, final worm weight, weigth loss
    5. MetaDiv_table: Metadata3, contains diversity calculations.
Code scripts and workflow 
Figures and statistics
    1. Fig1_BEF_funexp
    2. Fig2
    3. Fig3
    4. Fig4
    5. Fig5
    6. Fig6
    7. FigS2
    8. FigS3
    9. FigS3

SOFTWARE VERSIONS
[provide the version numbers of software (R, Julia, Python, Mathematica, etc) and loaded packages that you used to analyze your data files or run your simulations. If you used software that does not provide scripts (e.g. some popular statistical applications), please provide detailed information on the steps you used to perform the analyses and obtain the results reported in your paper]
REFERENCES
[references to papers referred to in this repository, if any]



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


