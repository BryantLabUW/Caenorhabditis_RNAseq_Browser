### Overview  
The *Caenorhabditis* RNA-seq Browser enables users to browse *Caenorhabditis* bulk RNA-seq data sets generated as part of the modENCODE project and perform on-demand differential expression analyses. Data from the following species are currently included in this app: *C. elegans*, *C. briggsae*, *C. brenneri*, *C. japonica*, and *C. remanei*. The app permits
browsing RNA-seq data in two modes:

1.  Browse by Gene Mode
2.  Browse by Life Stage Mode

### App Features  
Features of the app include:

* Search for gene(s) of interest using gene names, stable geneIDs, or keywords
* Extract gene expression values for genes of interest
  - Display gene expression across life stages as a heatmap (all
genes of interest) or a box plot (individual genes)
  - Display gene expression across life stages for individual genes and their known homologs in other *Caenorhabditis* species
-   Download log2 counts per million expression for genes of
interest as .xslx
-   On demand limma-voom-based pairwise differential gene expression
analysis
-   Display results as interactive volcano plots and data tables
-   Download results as .pdf (plots) or .xlsx (datatables)
-   Download raw/pre-processed data using user-friendly drop-down menu
-   Study design files (.csv)
-   Log2 counts per million expression for all genes and all samples
(.csv)
-   DGEList object (R object; primary data input
for the app)


To perform Gene Set Enrichment Analysis of differentially expressed genes, we recommend using one of the many tools already available, including:  
- [WormBase Gene Set Enrichment Analysis app](https://wormbase.org//tools/enrichment/tea/tea.cgi)
- [WormCat](http://www.wormcat.com/) ([Holdorf *et al* 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7017019/))