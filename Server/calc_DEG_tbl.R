## Calculate Differential Gene Expression for a specified pairwise comparison.
## Requires a Bayes-smoothed linear model
## Produces a Top Table of differentially expressed genes. 
## This function is called by limma_ranking.R file
calc_DEG_tbl <- function (ebFit, coef) {
    myTopHits.df <- limma::topTable(ebFit, adjust ="BH", 
                                    coef=coef, number=40000, 
                                    sort.by="logFC") %>%
        as_tibble() %>%
        dplyr::rename(tStatistic = t, LogOdds = B, BH.adj.P.Val = adj.P.Val) %>%
        dplyr::relocate(WormBaseID, Description, InterPro,
                        GS1_homologID, GS1_percent_homology,
                        GS2_homologID, GS2_percent_homology,
                        GS3_homologID, GS3_percent_homology,
                        GS4_homologID, GS4_percent_homology, .after = LogOdds) %>%
        dplyr::relocate(ends_with("geneName"), .after = geneID)

    myTopHits.df
}
