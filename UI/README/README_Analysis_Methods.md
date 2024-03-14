### On-demand Differential Gene Expression

For *C. elegans* and *C. briggsae*, the limma package ( [Ritchie *et al*
2015](https://pubmed.ncbi.nlm.nih.gov/25605792/), [Phipson *et al*
2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5373812/)) is used to
conduct pairwise differential gene expression (DEG) analyses between
life stages. DEG analyses are conducted on
variance-stabilized, filtered, normalized log2 counts per million (log2CPM) values. We use
empirical Bayes smoothing of gene-wise standard deviations to [provide
increased statistical power](https://www.degruyter.com/doi/10.2202/1544-6115.1027). All
reported p-values are adjusted for multiple comparisons using the
Benjamini-Hochberg method. If users are performing multiple closely
related pairwise comparisons they may select the option to correct for
multiple contrasts using the global method (using the slider input
provided). For each gene, we determine whether expression is up-regulated,
down-regulated, or not significantly different based on a cutoff of
p&lt;0.05.

### Data Visualization

For heat maps of log2CPM gene expression values, columns
(life stages) are ordered using Spearman clustering of the user-defined subset. Rows are ordered using Pearson
clustering of expression of the user-selected gene subset. For *C. elegans* and *C. briggsae*, gene expression plots use 
variance-stabilized, filtered, and normalized log2CPM data. For all other species, these plots show filtered and normalized log2CPM data.
