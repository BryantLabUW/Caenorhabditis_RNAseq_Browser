# Load pre-calculated data generated by pre-processing scripts
# 
# Import a file containing life stage legend terms
lifestage_legend <- suppressMessages(read_tsv("./Data/Life_stage_legend.txt",
                             # col_types = "cc",
         na = c("", "NA", "na")))

## Set Expression threshold values for plotting and saving DEGs ----
adj.P.thresh <- 0.05
lfc.thresh <- 1