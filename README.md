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

* [Fig1_BEF_funexp_main_bDiv_cFun.R](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun.R) contains the figure and stats to represent this general diversity-funciton relationship. 
* The stats include model selection as conducted in [Fig1_BEF_funexp_main_bDiv_cFun Stats.R](https://github.com/catalicu/FunExp/blob/main/Fig1_BEF_funexp_main_bDiv_cFun Stats.R).
The final figure was edited in power point for the position of the panels [Figure 1](https://github.com/catalicu/FunExp/blob/main/Figures/Fig1_wlegend.png)

Other files and figures reflect the remaining figures and statistics included in the manuscript. 

**Controls**  
To establish the diversity and function baselines before running the degradation experiment, we built a [Fig_DivBaseline.R](https://github.com/catalicu/FunExp/blob/main/Fig_DivBaseline.R). 
We assessed bacterial abundance using qPCR on the 16S rRNA gene and established whether abundance influences diversity and/or function in [qPCR_FunExp.R](https://github.com/catalicu/FunExp/blob/main/qPCR_FunExp.Rmd).  

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






