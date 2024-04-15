Sub folders contain multiple (3x) read count files for each biological sample. These replicates reflect three methods for obtaining read counts. 

For all files, genome and GFF/GTF files were downloaded from WormBase.org version WS290. Reads were aligned to each genome using `STAR` (v2.7.6a, `--alignIntronMax 30000 --alignMatesGapMax 30000`) and the species-specific WS290 GTF file for each genome. PCR duplicates were removed using `seldup` ([Warner *et al*., 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581053/)).

For files ending "CDS.unique.only.ws290.txt": read counts were obtained for each gene (CDS region only which is labeled as "CDS"  in *C. briggsae* and *C. elegans* and as "coding_exon" for *C. remanei*, *C. japonica*, and *C. brenneri*) using `featureCounts` (from the software package `subread-2.0.6`) using default settings. Only uniquely mapping reads were counted. These files are the ones used in the *Caenorhabditis* RNA-seq Browser.

For files ending with "exon.splitMultimappersEvenly.ws290.txt": read counts were obtained for full transcripts using the featureCounts290 options -M --fraction so that multimappers were counted by splitting them equally among all of the locations where they aligned.

For files ending "CDS.splitMultimappersEvennly.ws290.txt": read counts were obtained CDS regions using the featureCounts290 options -M --fraction so that multimappers were counted by splitting them equally among all of the locations where they aligned.

Alignment and counting was performed by LaDeana Hillier (Waterston Lab, UW). 