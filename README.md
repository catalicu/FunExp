This README.md file was updated on 06/27/2025. 

## A. Paper associated with this archive 
Title: 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves' 
Citation: To be updated upon acceptance 

Brief abstract: 
How biodiversity and ecosystem functions change with succession has proven to be difficult to predict. We hypothesize that community diversity and function may respond in opposite ways to successional drivers such as nutrient availability, species interactions or abiotic stress. The microbial communities within Sarracenia purpurea leaves perform degradation functions, providing essential nutrients to the plant, but we know little about how succession within the leaf influences bacterial diversity and degradation. We collected pitcher plant fluid from leaves aged 2 to 24 weeks to use in microcosm experiments. We added a common bacterivore to half of the replicated microcosms to establish whether predation is a central successional driver. We used amplicon sequencing and a degradation assay to quantify diversity and ecosystem function.

## B. Author information: 
Names, institutions of all authors: To be updated upon acceptance. 

## C. Contact information
Name: To be updated upon acceptance  
Address:To be updated upon acceptance  
Email:To be updated upon acceptance 

# DATA & CODE FILE OVERVIEW
This data repository consist of data files, code scripts, and this README document, with the following data and code filenames and variables. 

## Contents and narrative
This manuscript depicts the changes in bacterial diversity and function along the life of a pitcher plant leaf. 

A central figure in this study depicts the relationship betwen ASV richness and Degradation, showing a negative relationship between diversity and in pitcher plant bacterial communities. Figure 1 in the manuscript includes the relationship between biodiversity and function (Fig1a) as well as the change in diversity (Fig1b) and function (Fig1c) over time, as the pitcher plant ages. Here are the scripts that were used to create these figures and the associated statistics to test for significance:

* [Fig1 script](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun.R) contains the figure and stats to represent this general diversity-funciton relationship. 
* The stats include model selection as conducted in [Fig1 stats script](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun_Stats.R). The final figure was edited in power point for the position of the panels.  

Other files and figures reflect the remaining figures and statistics included in the manuscript: 

* We complemented the analysis using a Structural Equation Modeling (SEM) approach. The script includes the modeling and plotting steps: [SEM script](https://github.com/catalicu/FunExp/blob/main/FunExp_SEM2025.R). 
* We estimated bacterial abundance using 16S rRNA gene copy numbers. Abundance relationships with leaf age (Fig3a), ASV richness (Fig3b) and weight change (Fig3c) are assessed in script for [Fig3](https://github.com/catalicu/FunExp/blob/main/Fig3_qPCR_FunExp.R). The script includes steps to create the plot and calculate statistics. 
* We illustrated the change in community composition due to leaf age using an NMDS and a perMANOVA. The script to perform these steps can be found in the script for [Figure 4](https://github.com/catalicu/FunExp/blob/main/Fig4_NMDS_adonis.R).  
* We explored the role of specific ASVs in contributing to diversity and function. The script [AVS plotprep](https://github.com/catalicu/FunExp/blob/main/Fig6_FunASV_plotprep.R) sets up the data table in long format while the script [ASV plot](https://github.com/catalicu/FunExp/Fig6_FunASV_plotprep.R) performs the plotting and statistical calculations. 

## Controls and supplementary materials 
To establish the diversity and function baselines before running the degradation experiment, we built a [Fig_DivBaseline.R](https://github.com/catalicu/FunExp/blob/main/Fig_DivBaseline.R). 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in [qPCR_FunExp.R](https://github.com/catalicu/FunExp/blob/main/qPCR_FunExp.Rmd).  

## Data file descriptions
Data files are stored in the input_data folder: 
1. TAXtable_FunExp12021-03-24.txt: Taxonomy table with ASV sequence, Phylum, Class, Order, Family, Genus. It lists the full dataset with 1375 ASVs.  
2. ASVtable_FunExp12021-03-24.txt: Sample vs ASV table showing read abundances. It contains 97 samples and 1375 ASVs.  
3. METAtable_RunExp1230-03-24_sampleNamesFIX.txt: Table listing each sample's Leaf age code, Treatment code, and Replicate code. It contains 97 samples.  
4. FundExp_metadata.txt: Contains weight loss data including Tube_name, Treatment age code, Treatment code, Age code, Protozoan treatment, replicate code, initial worm weight, final worm weight, weigth loss. It contains 97 samples.  
5. MetaDiv_table_FunExp12021-03-30.txt: Table listing each sample's diversity calculations, including Richness, Shannon and Simpson's index and Pielou's evenness. It contains 97 samples.  
6. FunEx1_calculations.csv and FunEx2_calculations.csv: Datatables with calculated 16S rRNA gene copy counts from qPCR. 
7. ASVtable_forGLM_FunExp12023-01-04.txt: This table includes ASV abundance data from data file 1 and Metadata table from data file 3 for each sample. It contains 97 samples.  
8. Leavs_per_age.csv: This table lists information from individual leaves sampled from the field, before being pooled per age group. The dataset contains information on 457 leaves and includes each leave's age, band number (identifier tag), and fluid volume.  

Data files product of data analysis are stored in the output_table/ folder:
9. ASVslope_sorted52024-02-21.txt:
10. ASVrel.abund.pa.txt: This table is generated with script 6. It includes each ASV's average relative abundance, frequency, the slope of its relationship with function and its phylogenetic classification.  
11. ASVlist2.txt: This table is generated with script 6. It includes each ASV's code, average relative abundance, phylogenetic affiliations and its identifier sequence.



## Script descriptions
These scripts process data to generate the figures and statistical tests for the main text:
1. Fig1_BEF_funexp_main_bDiv_cFun.R: This script creates the three panels for Fig 1. The input for this script includes data files 2, 3, 4 and 5.  
2. Fig1_BEF_funexp_main_bDiv_cFun_Stats.R: Performs statistical testing to support patterns displayed in Fig1. Input used includes data files 2, 3, 4, and 5.  
3. Fig2_FunExp_SEM2025.R: Performs SEM analysis and generates Figure 2. Input files are data files 3, and 4.  
4. Fig3_qPCR_FunExp.Rmd: This script is in Rmarkdown. It creates figure 3 and performs statistical analyses associated with our qPCR data. It uses data file 6, as well as 2 and 7 for comparisons.  
5. Fig4_NMDS_adonis.R: Performs multivariate analyses including perMANOVA, NMDS ordination and homogeneity of distances. This code also generates Figure S4 and runs stats associated with it. Data files 2 and 3 are used here.  
6. Fig5_FunASV_plotprep.R: This script prepares data to create figures to identify ASV underlying function. It relies on data file 1, 7 and 9. This scripts generates data files 9 and 10, required for script 7.
7. Fig5_FunASV_plot.R: This script generates Figure 5. Using data files 9, 10, and 11. Data files 10 and 11 are generated with script 6, make sure to run it before script 7.

These scripts process data to generate supplementary materials:
8. FigS1and2_DivBaseline.R: This script generates figure S1 and 2. It performs the associated statistical testing. Requires data files 2, 3 and 8.
9. FigS3_and stats_FunAge_noOld.R: This script generates figure 3 and performs the associated statistical tests. It requires data files 4 and 5. 

SOFTWARE VERSIONS
The software and package versions can be found in the [session information](https://github.com/catalicu/FunExp/blob/main/sessionInfo_funexp04242025.txt) file. 

REFERENCES
[references to papers referred to in this repository, if any]






