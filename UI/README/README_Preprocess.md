### Data Source
The data sets included in the *Caenorhabditis* RNA-seq Browser are part of the modENCODE project ([Gerstein *et al* 2010](https://pubmed.ncbi.nlm.nih.gov/21177976/), and were provided by LaDeana Hillier (Waterson Lab, University of Washington). 

More information on the samples included in the browser for each species are available by downloading the study design files. 

### Read Alignment and Counting
Each genome and gff file were downloaded from WormBase.org version WS290. Reads were aligned to each genome using `STAR` (v2.7.6a, `--alignIntronMax 30000 --alignMatesGapMax 30000`) and the species-specific WS290 GTF file for each genome. PCR duplicates were removed using `seldup` ([Warner et al., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581053/)).  Read counts were obtained for each gene (CDS region only which is labeled as "CDS"  in *C. briggsae* and *C. elegans* and as "coding_exon" for *C. remanei*, *C. japonica*, and *C. brenneri*) using `featureCounts` (from the software package `subread-2.0.6`) using default settings. Only uniquely mapping reads were counted.  Additionally read counts were obtained for the CDS regions and for the full transcripts using the `featureCounts options -M --fraction` so that multimappers were counted by splitting them equally among all of the locations where they aligned. 

### Gene Annotation
Read data for each species was imported into R and annotated with information from WormBase ParaSite BiomaRT. 

Raw reads were quantified as counts per million using the EdgeR package, then filtered to remove transcripts with low counts (less than 1 count-per-million). A list of discarded genes and their expression values across life stages was saved. Non-discarded gene values were normalized using the trimmed mean of M-values method (TMM, [Robinson and Oshlack](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25) ) to permit between-samples comparisons. The mean-variance relationship was modeled using a precision weights approach [Law *et al* 2014](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29).  

### Filtering and Normalization Steps

Raw reads were quantified as log2 counts per million (CPM) using the `EdgeR` package,
then filtered to remove transcripts with low counts (less than 1
log2CPM in at least 1 sample). Non-discarded gene values are
normalized using the trimmed mean of M-values method [(TMM, Robinson and
Oshlack)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25)
to permit between-samples comparisons. The mean-variance relationship
was modeled using a precision weights approach [(Law *et al*
2014)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29).

These steps produced a variance-stabilized digital gene expression list (vDGEList), which serves as the primary input to the Shiny App.

### Saved Output Files

During pre-processing, the following required data files were saved and
are imported into Shiny application upon species initialization:

1.  a gene annotation R object (`_geneAnnotations`)
2.  the variance-stabilized vDGEList, saved as an R object
    (`_vDGEList`)
3.  a matrix of variance-stabilized gene expression data, extracted from
    the vDGEList (`log2cpm_filtered_norm_voom.csv`) - this data
    is downloadable from within the Browser App
