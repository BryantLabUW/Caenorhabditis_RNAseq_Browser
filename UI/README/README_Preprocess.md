### Read Alignment and Counting
Each genome and gff file were downloaded from WormBase.org version WS290. Reads were aligned to each genome using `STAR` (v2.7.6a, `--alignIntronMax 30000 --alignMatesGapMax 30000`) and the species-specific WS290 GTF file for each genome. PCR duplicates were removed using `seldup` ([Warner et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581053/)).  Read counts were obtained for each gene (CDS region only which is labeled as "CDS"  in *C. briggsae* and *C. elegans* and as "coding_exon" for *C. remanei*, *C. japonica*, and *C. brenneri*) using `featureCounts` (from the software package `subread-2.0.6`) using default settings. Only uniquely mapping reads were counted.  Additionally read counts were obtained for the CDS regions and for the full transcripts using the `featureCounts options -M --fraction` so that multimappers were counted by splitting them equally among all of the locations where they aligned. Alignment and counting was performed by LaDeana Hillier (Waterston Lab, UW). 

### Gene Annotation
Read data for each species was imported into R and annotated with information from WormBase ParaSite BiomaRT. 

### Filtering and Normalization Steps
Raw reads were quantified as log2 counts per million (CPM) using the `EdgeR` package,
then filtered to remove transcripts with low counts.  Non-discarded gene values are
normalized using the trimmed mean of M-values method (TMM, [Robinson and Oshlack](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25))
to permit between-samples comparisons. 

### Variance Stabilization and EList Generation
For *C. briggsae* and *C. elegans*, samples contain significant numbers of biological replicates. Thus these samples were additionally pre-processed for differential gene expression analyses. The mean-variance relationship was modeled using a precision weights approach ([Law *et al* 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)), using the `voom` function in the `limma` package. This function produced a variance-stabilized DGEList object that includes precision weights for each gene to control for heteroscedasity, as well as normalized expression values on the log2 scale (log2 counts-per-million). 

For *C. remanei*, *C. japonica*, and *C. brenneri*, the lack of consistent biological replicates precludes `limma::voom` processing and downstream differential expression calculations. To enable app users to visualize gene count data, a DGEList object was manually constructed that contained filtered, normalized log2CPM values, gene annotations information, and sample information. This file can be loaded into the RNA-seq Browser for plotting. 

For *C. elegans* embryonic timeline data, samples are from four previously published time series. Biological replicates are defined by a previously established inferred timeline that unified samples between these embryo time series ([Boeck *et al*., 2016](https://pubmed.ncbi.nlm.nih.gov/27531719/)). Group names represent the average time for biological replicate pairs, which were selected to represent 80 minute developmental intervals. This timing was chosen to match previously published differential expression analyses.

### Saved Output Files
During pre-processing, the following required data files were saved and
are imported into Shiny application upon species initialization:

1.  For all species, a gene annotation R object (`_geneAnnotations`).  
2.  For *C. elegans* and *C. briggsae*, a DGEList R object containing variance-stabilized, filtered and normalized log2CPM data
    (`_DGEList`), and a matrix of  variance-stabilized, filtered and normalized gene expression data (`log2cpm.csv`).  
3.  For *C. remanei*, *C. japonica*, and *C. brenneri*, a DGEList R object containing filtered and normalized log2CPM data (`_DGEList`), and a matrix of filtered and normalized gene gene expression data (`log2cpm.csv`).  

These data sets can be downloaded from within the Browser App.
