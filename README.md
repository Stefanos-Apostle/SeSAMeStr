
![SeSAMeStr_logo](www/SeSAMeStr_logo.jpg)

![Fig1](www/figure1.jpg)

SeSAMeStr Pipeline workflow. (A) Function and input requirements to run the pipeline. (B) Output directory structure with populated results. (C) Key QC figures generated in the preprocessing step (Detection rate, SNP allele frequency heatmap, Red-Green QQ Plot, Intensity Plot). (D) Figures generated during the PCA analysis based on input formula. (E) Key figures generated from DML analysis step based on input formula (Volcano plot, heatmap of significant CpG’s, GO Enrichment, Test Set Enrichment).

# SeSAMeStr

This pipeline is meant to act as a first pass analysis of DNA CpG Methylation using the SeSAMe R-Package.
This will take care of Quality Control, PCA, and Basic DML analysis. All files and plots will be output to 
output directory for ease of loading data back into R for either rerunning, replotting, or continuation of 
downstream analysis.

## Step 1: Install SeSAMeStr R-Package from Github

```
library(devtools)
devtools::install_github("Stefanos-Apostle/SeSAMeStr")
library(SeSAMeStr)
```

## Step 2: Fill out the SeSAMe_STREET_Sample_Sheet

Download from Github
```
wget https://github.com/Stefanos-Apostle/SeSAMeStr/blob/main/SeSAMe_STREET_Sample_Sheet.xlsx
```

## Step 3: Create the Output directory Architecture

```
.
├── DML
│   ├── DMR_Analysis
│   ├── GO_Enrichment
│   ├── Heatmaps
│   ├── testEnrichments
│   │	└──custom_sets
│   └── Volcano_plots
├── DimRed
└── QC
```

## Step 4: Run the SeSAMe_STREET() Function

```
SeSAMeStr(Idat_dir = "path_to/Idat_dir",
              out_dir = "path_to/output",
              sample_sheet = "path_to/SeSAMe_STREET_Sample_Sheet.xlsx",
              prep = "TQCDPB",
              formula = ~ Condition1 + Condition2 + ...,
              subsample = NA,
              cores = 4)
```

## Future Work

This has been built primarily for the MM285 platform, but is prepared to be extended to other platforms as well. This will happen if a personal/lab project requires it or if there is of interest by others in this extension.

## Cite

Cite SeSAMeStr v1.0.0 by;

Apostle, S., Fagnocchi, L., & Pospisilik, J. A. (2023). SeSAMeStr: An Automated Pipeline for SeSAMe Methylation Array Analysis (Version 1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7510575
