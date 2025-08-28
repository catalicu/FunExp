This README.md file was updated on 06/27/2025. 

## A. Paper associated with this archive 
Title: 'Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves' 
Citation: Cuellar-Gempeler, C., C. terHorst and T.E. Miller. Opposing effects of succession on bacterial diversity and function within pitcher plant (Sarracenia purpurea) leaves. The American Naturalist. In press.

Brief abstract: 
How biodiversity and ecosystem functions change with succession has proven to be difficult to predict. We hypothesize that community diversity and function may respond in opposite ways to successional drivers such as nutrient availability, species interactions or abiotic stress. The microbial communities within Sarracenia purpurea leaves perform degradation functions, providing essential nutrients to the plant, but we know little about how succession within the leaf influences bacterial diversity and degradation. We collected pitcher plant fluid from leaves aged 2 to 24 weeks to use in microcosm experiments. We added a common bacterivore to half of the replicated microcosms to establish whether predation is a central successional driver. We used amplicon sequencing and a degradation assay to quantify diversity and ecosystem function.

## B. Author information: 
Authors: Cuellar-Gempeler, Catalin (1), C. terHorst(2) and T. Miller(3). 
Affiliations: 
1 California State Polytechnic University, Humboldt, ORCID: 0000-0002-2634-9599
2 California State University, Northridge, ORCID: 0000-0001-7261-2178
3 Florida State University, ORCID: 0000-0002-2654-4079


## C. Contact information
Name: Catalina Cuellar-Gempeler
Address: 1st Harpst Street, Arcata CA. 95521. 
Email: ccg@humboldt.edu 

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
1. TAXtable_FunExp12021-03-24.txt: Taxonomy table with ASV sequence, Phylum, Class, Order, Family, Genus. ASV sequence is the full DNA sequence after denoising and merging sequences in dada2 with functions dada() and mergePairs() respectively. Taxonomic groupings were obtained by matching ASV sequences to GreenGenes using the assignTaxonomy() funciton in dada2. This table lists the full dataset with 1375 ASVs.  
2. ASVtable_FunExp12021-03-24.txt: Sample vs ASV table showing read abundances. It contains 97 samples and 1375 ASVs. Sample names reflect sample identity as described in METAtable_RunExp1230-03-24_sampleNamesFIX.txt. ASV sequence in the columns results from denoising and merging sequences in dada2 with functions dada() and mergePairs() respectively. Read abundance data within each cell is calculated with the makeSequenceTable() function in dada2.
3. METAtable_RunExp1230-03-24_sampleNamesFIX.txt: Table listing each sample's SampleID, sample_names_list1, Leaf age code, Treatment code, and Replicate code. It contains 97 samples. sample_names_list1 includes three pieces of information information: 1) on the age of the leaf in weeks, where A02 represents 2 weeks old leaves, and A11 represents 11 weeks old leaves, 2) on the treatment with TA0 meaning no predator added, TAI meaning predator added, and F50 representing the samples taken after filtering (pre-experimentation),  3) on the replicate identity, labeled R1, R2 or R3, B or C for the controls. This information is then parsed into the leaf_age (ex: A02, A04, A11, A12, etc), treatment (TA0, TAI or F50) and replicate columns. 
4. FundExp_metadata.txt: Contains weight loss data including Tube_name, Treatment age code, Treatment code, Age code, Protozoan treatment, replicate code, initial worm weight, final worm weight, weigth loss. Weight  data, including initial (T0), final (T1) and weight change (W_Change), are given in grams (g). Tube names, age codes and other column information reflect description in data file 3. This table contains 97 samples.  
5. MetaDiv_table_FunExp12021-03-30.txt: Table listing each sample's diversity calculations, including Richness, Shannon and Simpson's index and Pielou's evenness. Diversity calculations were made with functions Richness=specnumber(), Shannon and Simpson=diversity() (with index='shannon' or 'simpson' respectively) and Evenness=diversity()/log(specnumber()) in the package vegan. See script 1 for the code used in generating these calculations. Tube names, age codes and other column information reflect description in data file 3. This table contains 97 samples.  

6. ASVtable_forGLM_FunExp12023-01-04.txt: This table includes ASV abundance data from data file 1 and Metadata table from data file 3 for each sample. Column names and descriptions can be found within each of those entries. This table contains 97 samples.  
7. Leavs_per_age.csv: This table lists information from individual leaves sampled from the field, before being pooled per age group. The dataset contains information on 457 leaves and includes each leave's age (AGE_WEEKS), band number from the identifier tag (BAND)), and fluid volume in mL.  

Data files in the subfolder qPCR
8. FunEx1_calculations.csv and FunEx2_calculations.csv: Datatables with calculated 16S rRNA gene copy counts from qPCR. Columns include the Sample name (Sample), leaf age code (agetreatment), age, treatment, replicates, all of which follow the descriptions in data file 3. The last three columsn rely information regarding the qPCR output including Quantity (number of copies as obtained from fluorescence data calculated against the standard curve), Cells (number of cells calculated assuming 7 16S rRNA copies per cell - data not used in the manuscript), and melt.shape (describing the melt curve). The melt.shape column is used to convey and nterpret the shape of the melt curve as we interpret qPCR data outputs. 

Data files product of data analysis are stored in the output_table/ folder:  
9. ASVslope_sorted52024-02-21.txt (A) and ASVslope_sorted_time52024-02-21.txt (B): These tables present the output of GLMs that include  function as response variable and ASV relative abundance as explanatory variable for (A), and with leaf age (time) as an additional explanatory varialbe in (B). The data recorded from these model outputs are the direction of the Slope (positive, neutral or negative), the ASVcode (name generated to refer to each ASV see script 6), the slope estimate (ASVestimate), p value (ASVp.value), taxonomic affiliation for each taxa (Kingdom, Phylum, Class, Order, Family, Genus), and ASV sequence (ASVseq). See script 6 for details on the models and how this table is generated. 
10. ASVrel.abund.pa.txt: This table is generated with script 6. It includes each ASV's average relative abundance (ASVrel.abund), frequency (ASVpa.freq), frequency classification (Frequency), the slope of its relationship with function and its phylogenetic classification. Frequency is calculated as the average percentage samples each ASV appears in: 1 represents ASV present in 100% samples, while 0.03 represents ASV present in 3% samples. Frequency classification includes most frequent (the most frequent ASV), frequent (between 1 and 0.04), unfrequent (below 0.3). The slope estimate and p value are obtained from the GLMs, as described in data file 9. 
11. ASVlist2.txt: This table is generated with script 6. It includes each ASV's code (ASVcode), average relative abundance (Abund), phylogenetic affiliations (Kingdom, Phylum, Class, Order, Family, Genus) and its identifier sequence (ASVseq). 

## Script descriptions
These scripts process data to generate the figures and statistical tests for the main text:
1. Fig1_BEF_funexp_main_bDiv_cFun.R: This script creates the three panels for Fig 1. The input for this script includes data files 2, 3, 4 and 5.  
2. Fig1_BEF_funexp_main_bDiv_cFun_Stats.R: Performs statistical testing to support patterns displayed in Fig1. Input used includes data files 2, 3, 4, and 5.  
3. Fig2_FunExp_SEM2025.R: Performs SEM analysis and generates Figure 2. Input files are data files 3, and 4.  
4. Fig3_qPCR_FunExp.Rmd: This script is in Rmarkdown. It creates figure 3 and performs statistical analyses associated with our qPCR data. It uses data file 8, as well as 2 and 6 for comparisons.  
5. Fig4_NMDS_adonis.R: Performs multivariate analyses including perMANOVA, NMDS ordination and homogeneity of distances. This code also generates Figure S4 and runs stats associated with it. Data files 2 and 3 are used here.
6. FunASV_abund-funct_B.Rmd (A) and FunASV_time-abund-funct.Rmd (B): These scripts evaluate each ASV for its contributions to function by building GLMs that include function as response variable and ASV relative abundance as and explanatory variable. Script (B) includes time as an explanatory factor as well. These scripts generate data files in 9. 
7. Fig5_FunASV_plotprep.R: This script prepares data to create figures to identify ASV underlying function. It relies on data file 1, 6 and 9. This scripts generates data file 10 and 11.
8. Fig5_FunASV_plot.R: This script generates Figure 5. Using data files 10, and 11. Data files 10 and 11 are generated with script 6, make sure to run it before script 7.


These scripts process data to generate supplementary materials:  
9. FigS1and2_DivBaseline.R: This script generates figure S1 and 2. It performs the associated statistical testing. Requires data files 2, 3 and 8.  
10. FigS3_and stats_FunAge_noOld.R: This script generates figure 3 and performs the associated statistical tests. It requires data files 4 and 5. 

# SOFTWARE VERSIONS
The software and package versions can be found in the [session information](https://github.com/catalicu/FunExp/blob/main/sessionInfo_funexp04242025.txt) file. 

R version 4.5.0 (2025-04-11)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] effects_4.2-2        outliers_0.15        lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
 [6] purrr_1.0.4          readr_2.1.5          tidyr_1.3.1          tibble_3.2.1         tidyverse_2.0.0     
[11] FSA_0.9.6            car_3.1-3            carData_3.0-5        gridExtra_2.3        nlme_3.1-168        
[16] lme4_1.1-37          plyr_1.8.9           phyloseq_1.52.0      metagenomeSeq_1.50.0 RColorBrewer_1.1-3  
[21] glmnet_4.1-8         Matrix_1.7-3         limma_3.64.0         Biobase_2.68.0       BiocGenerics_0.54.0 
[26] generics_0.1.3       dada2_1.36.0         Rcpp_1.0.14          dplyr_1.1.4          knitr_1.50          
[31] ggplot2_3.5.2        vegan_2.6-10         lattice_0.22-6       permute_0.9-7       

loaded via a namespace (and not attached):
  [1] rstudioapi_0.17.1           jsonlite_2.0.0              shape_1.4.6.1              
  [4] magrittr_2.0.3              dunn.test_1.3.6             estimability_1.5.1         
  [7] farver_2.1.2                nloptr_2.2.1                rmarkdown_2.29             
 [10] vctrs_0.6.5                 multtest_2.64.0             minqa_1.2.8                
 [13] Rsamtools_2.24.0            htmltools_0.5.8.1           S4Arrays_1.8.0             
 [16] survey_4.4-2                Rhdf5lib_1.30.0             SparseArray_1.8.0          
 [19] Formula_1.2-5               rhdf5_2.52.0                KernSmooth_2.23-26         
 [22] GenomicAlignments_1.44.0    igraph_2.1.4                lifecycle_1.0.4            
 [25] iterators_1.0.14            pkgconfig_2.0.3             R6_2.6.1                   
 [28] fastmap_1.2.0               GenomeInfoDbData_1.2.14     rbibutils_2.3              
 [31] MatrixGenerics_1.20.0       digest_0.6.37               colorspace_2.1-1           
 [34] ShortRead_1.66.0            S4Vectors_0.46.0            GenomicRanges_1.60.0       
 [37] hwriter_1.3.2.1             labeling_0.4.3              timechange_0.3.0           
 [40] httr_1.4.7                  abind_1.4-8                 mgcv_1.9-1                 
 [43] compiler_4.5.0              withr_3.0.2                 BiocParallel_1.42.0        
 [46] DBI_1.2.3                   gplots_3.2.0                MASS_7.3-65                
 [49] DelayedArray_0.34.1         biomformat_1.36.0           gtools_3.9.5               
 [52] caTools_1.18.3              Wrench_1.26.0               tools_4.5.0                
 [55] ape_5.8-1                   nnet_7.3-20                 glue_1.8.0                 
 [58] rhdf5filters_1.20.0         grid_4.5.0                  cluster_2.1.8.1            
 [61] reshape2_1.4.4              ade4_1.7-23                 gtable_0.3.6               
 [64] tzdb_0.5.0                  hms_1.1.3                   data.table_1.17.0          
 [67] utf8_1.2.4                  XVector_0.48.0              foreach_1.5.2              
 [70] pillar_1.10.2               mitools_2.4                 splines_4.5.0              
 [73] survival_3.8-3              deldir_2.0-4                tidyselect_1.2.1           
 [76] locfit_1.5-9.12             Biostrings_2.76.0           reformulas_0.4.0           
 [79] IRanges_2.42.0              SummarizedExperiment_1.38.0 stats4_4.5.0               
 [82] xfun_0.52                   statmod_1.5.0               matrixStats_1.5.0          
 [85] stringi_1.8.7               UCSC.utils_1.4.0            yaml_2.3.10                
 [88] boot_1.3-31                 evaluate_1.0.3              codetools_0.2-20           
 [91] interp_1.1-6                BiocManager_1.30.25         cli_3.6.4                  
 [94] RcppParallel_5.1.10         Rdpack_2.6.4                munsell_0.5.1              
 [97] GenomeInfoDb_1.44.0         png_0.1-8                   parallel_4.5.0             
[100] latticeExtra_0.6-30         jpeg_0.1-11                 bitops_1.0-9               
[103] pwalign_1.4.0               scales_1.3.0                insight_1.2.0              
[106] crayon_1.5.3                rlang_1.1.6                

# REFERENCES
See full manuscritp for list of references.






