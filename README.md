# Caenorhabditis_RNAseq_Browser 
Web-based Shiny App for browsing and on-demand analysis of *Caenorhabditis* RNA-seq data from the modENCODE project ([Gerstein *et al* 2010](https://pubmed.ncbi.nlm.nih.gov/21177976/)).  
[This app is deployed via Shinyapps.io](https://bryantlabuw.shinyapps.io/Caenorhabditis_RNAseq_Browser/)

This app is based on the *Strongyloides* RNA-seq Browser. For more information, [please see the G3 paper associated with the original project](https://pubmed.ncbi.nlm.nih.gov/33823530/).

## Table of Contents  
1. [General Information](#general-information)
2. [App Setup & Deployment](#app-setup-&-deployment)
3. [App Features](#app-features)
4. [Sources](#sources)
5. [License](#license)
6. [Authors](#authors)

## General Information
This repository contains source code for the web-based *Caenorhabditis* RNA-seq Browser. This app is deployed via Shinyapps.io but can also be run locally. See App Setup and App Features sections below for additional details.  

The sections below describe the contents of the primary subfolders within this repository:

### Data  
This folder contains pre-processed data files, including study design files, gene annotations, and digital gene expression lists (vDGEList) containing variance-stabilized, filtered, TMM-normalized RNA-seq data.

### Server
Server files for the Shiny app.

### UI
User interface files for the Shiny app. Includes custom css and additional README files with methods details.

### www
Static files that can be interactively downloaded within the *Caenorhabditis* RNA-seq Browser environment.

## App Setup & Deployment
To access a stable deployment of the *Caenorhabditis* RNA-seq Browser Web App, please visit:   [hallemlab.shinyapps.io/strongyloides_rnaseq_browser/](https://hallemlab.shinyapps.io/Strongyloides_RNAseq_Browser/)  


To run the latest version locally from GitHub, use the following command in R/RStudio:  
`library(shiny)`  
`shiny::runGitHub(repo = 'Caenorhabditis_RNAseq_Browser', username = 'BryantLabUW')`  

To run a specific release locally use the following commands in R/RStudio:  
  * For PCs --  
    `library(shiny)`  
    `shiny::runUrl('https://github.com/BryantLabUW/Caenorhabditis_RNAseq_Browser/archive/<RELEASE_VERSION>.zip') ` 

  * For Macs --  
    `library(shiny)`  
    `shiny::runUrl('https://github.com/BryantLabUW/Caenorhabditis_RNAseq_Browser/archive/<RELEASE_VERSION>.tar.gz')`  

Please note: the download step for runURL/runGitHub takes a substantial amount of time. We recommend downloading this archive and running the application locally. 

## App Features  
The *Caenorhabditis* RNA-seq Shiny Browser enables users to browse *Caenorhabditis* bulk RNA-seq datasets generated as part of the modENCODE project and perform on-demand analyses including differential expression and gene set enrichment. Data from the following species are currently included: *C. elegans*, *C. briggsae*, *C. brenneri*, *C. japonica*, and *C. remanei*. The app permits browsing RNA-seq data in two modes:

  1. Browse by Life Stage Mode
  2. Browse by Gene Mode  
  
Features of the app include:  

* Search for gene(s) of interest using stable geneIDs or keywords
* Extract gene expression values for genes of interest
  - Display gene expresion across life stages as a heatmap (all genes of interest) or a boxplot (individual genes)
  - Display gene expression across life stages for individual genes and their known *Caenorhabditis* homologs
  - Download log2 counts per million expression for genes of interest as .xslx
* On demand limma-voom-based pairwise differential gene expression analysis
  - Display results as interactive volcano plots and datatables
  - Download results as .pdf (plots) or .xlsx (datatables)
* Download raw/pre-processed data using user-friendly dropdown menu
  - Study design files (.csv)
  - Log2 counts per million expression for all genes and all samples (.csv)
  - Variance-stabilized DGEList object (R object; primary data input for the app)


## Sources
* [Shiny](https://shiny.rstudio.com/) - UI framework
* [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) - Differential gene expression
* *Caenorhabditis* RNA-seq datasets:
  - [Gerstein *et al* 2010](https://pubmed.ncbi.nlm.nih.gov/21177976/)
* [WormBase ParaSite](https://parasite.wormbase.org/index.html) - Gene annotations
* [DIYTranscriptomics](http://diytranscriptomics.com/) - Virtual asynchronous course where the authors learned best practices for RNA-seq data analysis; provided primary pipeline for data pre-processing and analysis

## License  
This project is licensed under the MIT License. 

## Authors  
* [Astra Bryant, PhD](https://github.com/astrasb)

