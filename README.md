This README.txt file was updated on 04/16/2025. 

## A. Paper associated with this archive 
Title: 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves' 
Citation: To be updated upon acceptance 

Brief abstract: 
How biodiversity and ecosystem functions change with succession has proven to be difficult to predict. We hypothesize that community diversity and function may respond in opposite ways to successional drivers such as nutrient availability, species interactions or abiotic stress. The microbial communities within Sarracenia purpurea leaves perform degradation functions, providing essential nutrients to the plant, but we know little about how succession within the leaf influences bacterial diversity and degradation. We collected pitcher plant fluid from leaves aged 2 to 24 weeks to use in microcosm experiments. We added a common bacterivore to half of the replicated microcosms to establish whether predation is a central successional driver. We used amplicon sequencing and a degradation assay to quantify diversity and ecosystem function.

## B. Author information: 
Names, institutions of all authors: To be updated upon acceptance 

## C. Contact information
Name: To be updated upon acceptance
Address:To be updated upon acceptance
email:To be updated upon acceptance 

# DATA & CODE FILE OVERVIEW
This data repository consist of data files, code scripts, and this README document, with the following data and code filenames and variables. 

## Contents and narrative
This manuscript depicts the changes in bacterial diversity and function along the life of a pitcher plant leaf. 

A central figure in this study depicts the relationship betwen ASV richness and Degradation, showing a negative relationship between diversity and in pitcher plant bacterial communities. Figure 1 in the manuscript includes the relationship between biodiversity and function (Fig1a) as well as the change in diversity (Fig1b) and function (Fig1c) over time, as the pitcher plant ages. Here are the scripts that were used to create these figures and the associated statistics to test for significance:

* [Fig1 script](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun.R) contains the figure and stats to represent this general diversity-funciton relationship. 
* The stats include model selection as conducted in [Fig1 stats script](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun_Stats.R). The final figure was edited in power point for the position of the panels [Figure 1](https://github.com/catalicu/FunExp/blob/main/Figures/Fig1_wlegend.png). 

Other files and figures reflect the remaining figures and statistics included in the manuscript: 

* We complemented the analysis using a Structural Equation Modeling (SEM) approach. The script includes the modeling and plotting steps: [SEM script](https://github.com/catalicu/FunExp/blob/main/FunExp_SEM2025.R). Figure 2 can be seen [HERE](https://github.com/catalicu/FunExp/blob/main/Figures/Fig2_SEM.png). 
* We estimated bacterial abundance using 16S rRNA gene copy numbers. Abundance relationships with leaf age (Fig3a), ASV richness (Fig3b) and weight change (Fig3c) are assessed in script for [Fig3](https://github.com/catalicu/FunExp/blob/main/Fig3_qPCR_FunExp.R). The script includes steps to create the plot and calculate statistics. Figure 3 can be seen [HERE](https://github.com/catalicu/FunExp/blob/main/Figures/Fig3_qPCR.tiff). 
* We illustrated the change in community composition due to leaf age using an NMDS and a perMANOVA. The script to perform these steps can be found in the script for [Figure 4](https://github.com/catalicu/FunExp/blob/main/Fig4_NMDS_adonis.R). Figure 4 can be seen [HERE](https://github.com/catalicu/FunExp/blob/main/Figures/Fig4_NMDS_small.pdf). 
* We explored the role of specific ASVs in contributing to diversity and function. The script [AVS plotprep](https://github.com/catalicu/FunExp/blob/main/Fig6_FunASV_plotprep.R) sets up the data table in long format while the script [ASV plot](https://github.com/catalicu/FunExp/Fig6_FunASV_plotprep.R) performs the plotting and statistical calculations. The resulting plots can be seen in [Figure 6](https://github.com/catalicu/FunExp/blob/main/Figures/Fig6.png). 

## Controls and supplementary materials 
To establish the diversity and function baselines before running the degradation experiment, we built a [Fig_DivBaseline.R](https://github.com/catalicu/FunExp/blob/main/Fig_DivBaseline.R). 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in [qPCR_FunExp.R](https://github.com/catalicu/FunExp/blob/main/qPCR_FunExp.Rmd).  

## Data file descriptions
Data files are stored in the input_data folder
    1. TAXtable Taxonomy table
    2. ASVtable Species list
    3. METAtable Metadata1: Leaf age code, Treatment code, Replicate code.
    4. FundExp_metadata: Metadata2, contains weight loss data: Tube_name, Treatment age code, Treatment code, Age code, Protozoan treatment, replicate code, initial worm weight, final worm weight, weigth loss
    5. MetaDiv_table: Metadata3, contains diversity calculations.

SOFTWARE VERSIONS
The software and package versions can be found in the [session information](https://github.com/catalicu/FunExp/blob/main/sessionInfo_funexp04242025.txt) file. 

REFERENCES
[references to papers referred to in this repository, if any]






