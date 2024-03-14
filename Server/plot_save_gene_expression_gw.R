## GW: Heatmap/Table of Gene Expression ----
generateHeatmapTable <- reactive({
  req(vals$genelist,input$displayedGene,vals$genelist.Log2CPM,vals$DGEList.filtered.norm)
  
  # Set gene to display
  vals$gene_of_interest <- vals$genelist$geneID

  if (input$displayedGene == "Data Table") {
   
    excluded.genes <- dplyr::anti_join(vals$submitted.genelist, 
                                       vals$genelist,
                                       by = "geneID") %>%
      left_join(vals$annotations, by = "geneID") # Add gene annotations
    setProgress(0.2)
    source('Server/switch_species.R', local = T)
    gene_vals <- vals$genelist.Log2CPM %>%
      dplyr::filter(geneID %in% vals$gene_of_interest) %>%
      dplyr::summarize(mean = mean(log2CPM), .groups = "drop_last") %>%
      pivot_wider(names_from = life_stage,
                  id_cols = geneID,
                  values_from = mean) %>%
      left_join(vals$annotations, by = "geneID") %>%
        dplyr::mutate(WormBaseLink = paste0("<a href='https://wormbase.org/species/C_", species, "/gene/", geneID,"' target = '_blank'>", geneID,"</a>"))%>%
        dplyr::mutate(GS1_homologID = paste0("<a href='https://wormbase.org/species/C_", species.GS1, "/gene/", sub("\\S* \\| " ,"", GS1_homologID),"' target = '_blank'>", GS1_homologID,"</a>"))%>%
        dplyr::mutate(GS2_homologID = paste0("<a href='https://wormbase.org/species/C_", species.GS2, "/gene/", sub("\\S* \\| " ,"", GS2_homologID),"' target = '_blank'>", GS2_homologID,"</a>"))%>%
        dplyr::mutate(GS3_homologID = paste0("<a href='https://wormbase.org/species/C_", species.GS3, "/gene/", sub("\\S* \\| " ,"", GS3_homologID),"' target = '_blank'>", GS3_homologID,"</a>"))%>%
        dplyr::mutate(GS4_homologID = paste0("<a href='https://wormbase.org/species/C_", species.GS4, "/gene/", sub("\\S* \\| " ,"", GS4_homologID),"' target = '_blank'>", GS4_homologID,"</a>"))%>%
        dplyr::relocate(geneName, .after = geneID) %>%
        dplyr::relocate(WormBaseLink, .before = Description) %>%
         {suppressMessages(dplyr::full_join(.,excluded.genes))} 
    
    n_num_cols <- ncol(gene_vals)
    n_num_values <- nlevels(vals$DGEList.filtered.norm$targets$group)
    setProgress(0.4)
    gene_vals.datatable <- gene_vals %>%
      DT::datatable(rownames = FALSE,
                    escape = FALSE,
                    options = list(autoWidth = TRUE,
                                   scrollX = TRUE,
                                   scrollY = '300px',
                                   scrollCollapse = TRUE,
                                   searchHighlight = TRUE, 
                                   pageLength = 10, 
                                   lengthMenu = c("5",
                                                  "10",
                                                  "25",
                                                  "50",
                                                  "100"),
                                   initComplete = htmlwidgets::JS(
                                     "function(settings, json) {",
                                     paste0("$(this.api().table().container()).css({'font-size': '", "10pt", "'});"),
                                     "}"),
                                   columnDefs = list(
                                     list(
                                       targets = c(n_num_values+4, n_num_values+5),
                                       render = JS(
                                         "function(data, type, row, meta) {",
                                         "return type === 'display' && data.length > 20 ?",
                                         "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                         "}")
                                     ),
                                     list(targets = "_all",
                                          class="dt-right")
                                   ),
                                   rowCallback = JS(c(
                                     "function(row, data){",
                                     "  for(var i=0; i<data.length; i++){",
                                     "    if(data[i] === null){",
                                     "      $('td:eq('+i+')', row).html('NA')",
                                     "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
                                     "    }",
                                     "  }",
                                     "}"
                                   )
                                   )
                                   
                    )
      )
    setProgress(0.8)
    gene_vals.datatable <-  gene_vals.datatable %>%
      DT::formatRound(columns=c(3:(n_num_values+2)), 
                      digits=3)
    
    gene_vals.datatable <-  gene_vals.datatable %>%
      DT::formatRound(columns=c(n_num_cols, 
                                n_num_cols-2,
                                n_num_cols-4,
                                n_num_cols-6), 
                      digits=2)
    
  } else {
    # Make a heatmap for all the genes using the Log2CPM values
    myheatcolors <- RdBu(75)
    
    diffGenes <- vals$diffGenes.df %>%
      dplyr::select(!geneID) %>%
      as.matrix()
    rownames(diffGenes) <- rownames(vals$DGEList.filtered.norm$E)
    subset.diffGenes<- diffGenes[vals$gene_of_interest,]

    renaming.rows <- vals$DGEList.filtered.norm$genes[rownames(subset.diffGenes),c("geneID", "geneName")]
    rownames(subset.diffGenes) <- paste0(renaming.rows$geneID, " | ", renaming.rows$geneName)

    setProgress(0.2)
    options(warn = 2)
    spearmancorr <- try(cor(subset.diffGenes, method="spearman"), silent = T)
    
    validate(need(!is(spearmancorr, "try-error"), 'Error: heatmap column clustering cannot be performed by likely because the expression levels of the selected genes are too similar in too many life stages.
                  Please try again with more genes, or switch to single-gene plots.'))
    clustColumns <- hclust(as.dist(1-cor(subset.diffGenes, method="spearman")), method="complete")
    options(warn = 1)

    setProgress(0.4)
    
    clustRows <- hclust(as.dist(1-cor(t(subset.diffGenes), 
                                      method="pearson")), 
                        method="complete") 
    par(cex.main=1.2)
    vals$HeatmapRowOrder <- order.dendrogram(ladderize(as.dendrogram(clustRows)))
    vals$HeatmapColOrder <- order.dendrogram(ladderize(seriate_dendrogram(as.dendrogram(clustColumns),
                                                                          as.dist(1-cor(diffGenes, method="spearman")))))
    setProgress(0.6)
    
    hovertext <- as.data.frame(subset.diffGenes) %>%
      round(digits = 2)
    
    colnames(hovertext) <- vals$DGEList.filtered.norm$targets$samples
    hovertext[] <- lapply(seq_along(hovertext), function(x){
      paste0("GeneID: ", rownames(hovertext), "<br>",
             "Log2CPM: ", hovertext[,x], "<br>",
             "Life Stage: ", vals$DGEList.filtered.norm$targets$group[x],
             "<br>",
             "Sample: ", colnames(hovertext)[x])
    })
    
    showticklabels <- if(length(vals$gene_of_interest)<20){c(TRUE,TRUE)} else {c(TRUE,FALSE)}
    setProgress(0.8)
    p <- heatmaply(subset.diffGenes,
                   colors = rev(myheatcolors),
                   Rowv= ladderize(as.dendrogram(clustRows)),
                   Colv=ladderize(as.dendrogram(clustColumns)),
                   show_dendrogram = c(TRUE, TRUE),
                   showticklabels = showticklabels,
                   scale='row', #rows are scaled to have mean zero and standard deviation one. 
                   plot_method = "plotly",
                   branches_lwd = 0.2,
                   key.title = "Row-scaled Z Score",
                   cexRow=1.2, cexCol=1.2,
                   margins = c(100, 50, 10, 0),
                   colorbar_len = 0.5,
                   colorbar_ypos = 0.5,
                   colorbar_xpos = 1,
                   custom_hovertext = hovertext)
    setProgress(1)
    p
  }
})

## GW: Single Gene Expression Plots ----
generateGenePlot <- reactive({
  req(vals$genelist.Log2CPM, vals$genelist)
  vals$gene_of_interest <- input$displayedGene
  setProgress(0.25)
  
  #Plot Log2CPM values for an individual gene
  gene_vals <- vals$genelist.Log2CPM %>%
    left_join(vals$annotations, by = "geneID") %>%
    dplyr::select(geneID, life_stage, log2CPM, geneName) %>%
    dplyr::mutate(geneID = paste0(geneID, " | ", geneName)) %>%
    dplyr::filter(geneID == vals$gene_of_interest)
  setProgress(0.5)
  # Identify the identity of the primary species
  source('Server/switch_species.R', local = T)
  if (species == 'elegans'){gene_vals$life_stage <- factor(gene_vals$life_stage, levels = c("EE_028", "EE_118", "EE_199","EE_285", "EE_378", "EE_467", "EE_574", "EE_674", "L1", "L2", "L3", "L4", "L4_male", "YA", "Dauer", "DauerEntry", "DauerExit") )
  }
  p <- suppressWarnings (ggplot(gene_vals) + 
                           aes(x = life_stage, y = log2CPM, 
                               fill = life_stage) +
                           stat_boxplot(geom = "errorbar", 
                                        width = 0.25,
                                        size = .5,
                                        show.legend = F) +
                           geom_boxplot(show.legend = F, 
                                        alpha = 0.5) +
                           geom_jitter(size = 3,
                                       stroke = .75,
                                       show.legend = F,
                                       aes(color = life_stage)) +
                           labs(y="log2 CPM expression", x = "Life Stage",
                                title= "Log2 Counts per Million (CPM) Expression",
                                subtitle=paste("Selected gene:",
                                               gene_vals$geneID[1])) +
                           theme_Publication() + 
                           theme(aspect.ratio=2/3))
  
  # Identify the identity of the primary species
  source('Server/switch_species.R', local = T)

  if (species == 'elegans') {
    p <- p +
      theme(axis.text.x = element_text(
              angle = 45,
              hjust = 1),
            aspect.ratio=2/8)
  } 
  
  vals$gene_plot <- p
  p
})


## GW: Gene Homologs Expression Plots ----
##  Take a list of user-provided genes, identify gene homologs in other *Caenorhabditis*
## species, collect expression data from those gene homologs, and plot
fetch_homologs <- reactive({
  req(vals$genelist.Log2CPM, vals$genelist)
 
  Primary.species <- vals$genelist.Log2CPM %>%
    dplyr::filter(geneID %in% vals$genelist$geneID) %>%
    left_join(vals$annotations, by = "geneID") %>%
    dplyr::select(geneID, life_stage,log2CPM,
                  GS1_homologID, GS2_homologID,
                  GS3_homologID, GS4_homologID)

  genelist.allspecies <- Primary.species %>%
    ungroup() %>%
    dplyr::select(geneID, GS1_homologID, GS2_homologID, GS3_homologID, GS4_homologID) %>%
    unique() %>%
    separate_wider_delim(GS1_homologID:GS4_homologID, " | ", names = c(NA, ""),names_sep = "")
  
  vals$homologous_genes <- genelist.allspecies
 
  # Identify the identity of the primary species, and the homologous species
  source('Server/switch_species.R', local = T)
  
  # Load expression data for GS1-4 species
  load(file = paste0("./Data/",species.GS1,"_DGEList"))
  species.GS1.Log2CPM<-DGEList.filtered.norm$E %>%
    as_tibble(rownames = "geneID")%>%
    setNames(nm = c("geneID", 
                    as.character(DGEList.filtered.norm$targets$group))) %>%
    pivot_longer(cols = -geneID,
                 names_to = "life_stage", 
                 values_to = "log2CPM") %>%
    group_by(geneID, life_stage) %>%
    dplyr::filter(geneID %in% genelist.allspecies$GS1_homologID)
  remove(DGEList.filtered.norm)
  
  load(file = paste0("./Data/",species.GS2,"_DGEList"))
  species.GS2.Log2CPM<-DGEList.filtered.norm$E %>%
    as_tibble(rownames = "geneID")%>%
    setNames(nm = c("geneID", 
                    as.character(DGEList.filtered.norm$targets$group))) %>%
    pivot_longer(cols = -geneID,
                 names_to = "life_stage", 
                 values_to = "log2CPM") %>%
    group_by(geneID, life_stage) %>%
    dplyr::filter(geneID %in% genelist.allspecies$GS2_homologID)
  remove(DGEList.filtered.norm)
  
  load(file = paste0("./Data/",species.GS3,"_DGEList"))
  species.GS3.Log2CPM<-DGEList.filtered.norm$E %>%
    as_tibble(rownames = "geneID")%>%
    setNames(nm = c("geneID", 
                    as.character(DGEList.filtered.norm$targets$group))) %>%
    pivot_longer(cols = -geneID,
                 names_to = "life_stage", 
                 values_to = "log2CPM") %>%
    group_by(geneID, life_stage) %>%
    dplyr::filter(geneID %in% genelist.allspecies$GS3_homologID)
  remove(DGEList.filtered.norm)
  
  load(file = paste0("./Data/",species.GS4,"_DGEList"))
  species.GS4.Log2CPM<-DGEList.filtered.norm$E %>%
    as_tibble(rownames = "geneID")%>%
    setNames(nm = c("geneID", 
                    as.character(DGEList.filtered.norm$targets$group))) %>%
    pivot_longer(cols = -geneID,
                 names_to = "life_stage", 
                 values_to = "log2CPM") %>%
    group_by(geneID, life_stage) %>%
    dplyr::filter(grepl(genelist.allspecies$GS4_homologID[1], geneID))
  remove(DGEList.filtered.norm)

  life_stage_types <- lifestage_legend %>%
    dplyr::select(-group) %>%
    colnames()

  plot.tbl <- bind_rows(
    Primary.species =vals$genelist.Log2CPM,
    GS1.species = species.GS1.Log2CPM,
    GS2.species = species.GS2.Log2CPM,
    GS3.species = species.GS3.Log2CPM,
    GS4.species = species.GS4.Log2CPM,
    .id = "id") %>%
    dplyr::mutate(across("id", str_replace_all, c("Primary.species" = paste0("C. ",species), 
                                                  "GS1.species" = paste0("C. ",species.GS1),
                                                  "GS2.species" = paste0("C. ",species.GS2),
                                                  "GS3.species" = paste0("C. ",species.GS3),
                                                  "GS4.species" = paste0("C. ",species.GS4)
                                                  ))) %>%
    dplyr::mutate(id = factor(id, levels = c(paste0("C. ",species),
                                             paste0("C. ",species.GS1),
                                             paste0("C. ",species.GS2),
                                             paste0("C. ",species.GS3),
                                             paste0("C. ",species.GS4))))%>%
    dplyr::mutate(life_stage= factor(life_stage, levels = life_stage_types))
  
})

## GW: Set Up UI Outputs ----
output$CPM.homologs <- renderPlot({
  req(input$displayedGene != "All Genes")
  req(input$displayedGene != "Data Table")
  withProgress({
    plot.tbl <- fetch_homologs()
    
    gene_lookup <- dplyr::filter(vals$genelist, 
                                 geneReference %in% input$displayedGene)
    set_displayed <- dplyr::filter(vals$homologous_genes, 
                                   geneID %in% gene_lookup$geneID)%>%
      as.character()
    
    plot.tbl <- plot.tbl %>%
      dplyr::filter(geneID %in% set_displayed) 
    
    mylevels <- unique(plot.tbl[order(plot.tbl$id), "geneID"])
    plot.tbl <- plot.tbl %>%
      mutate(geneID = factor(geneID, levels = mylevels$geneID)) %>%
      group_by(id, geneID)
    
    p<-suppressWarnings (ggplot(plot.tbl) + 
                           aes(x = life_stage, y = log2CPM, fill = life_stage) +
                           stat_boxplot(geom = "errorbar", 
                                        width = 0.25,
                                        size = .5,
                                        show.legend = F) +
                           geom_boxplot(show.legend = F, 
                                        alpha = 0.5) +
                           geom_jitter(size = 3,
                                       stroke = .75,
                                       show.legend = F,
                                       aes(color = life_stage)) +
                           labs(y="log2 CPM expression", x = "Life Stage",
                                title= "Expression of Homologous Genes Across Species",
                                subtitle=paste("Selected gene:",
                                               vals$gene_of_interest)) +
                           facet_grid(~id + geneID, 
                                      scales = "free_x",
                                      space = "free_x") +
                           theme_Publication() + 
                           theme(axis.text.x = element_text(
                                   angle = 45,
                                   hjust = 1))
    )
   
    vals$homolog_plot <- p
    p
  }, message = "Loading Homolog Plot")
})

output$CPM <- renderPlot({
  req(input$displayedGene != "All Genes")
  req(input$displayedGene != "Data Table")
  withProgress(generateGenePlot(), message = "Loading Plot")
})


output$CPM.datatable <- renderDT({
  req(input$displayedGene == "Data Table")
  withProgress(generateHeatmapTable(), message = "Loading Data Table")
})

output$CPMPlotly <- renderPlotly({
  req(input$displayedGene == "All Genes")
  vals$genelist
  withProgress(generateHeatmapTable(), message = "Loading Heatmap")
})

## GW: Switch what type of GW Plot is produced ----
observe({
  req(vals$genelist)
  req(input$displayedGene)
  vals$genelist
  
  if(input$displayedGene == "All Genes"){
    removeUI(
      selector = "#CPMPlotlydiv"
    )
    removeUI(
      selector = "#CPMdiv"
    )
    
    removeUI(
      selector = "#CPMTablediv"
    )
    
    insertUI(
      selector = '#GenePlotDiv',
      where = "beforeBegin",
      ui = tagList(div(id = "CPMPlotlydiv",
                       h5("Log2 Counts per Million (CPM) Expression Across Life Stages"),
                       plotlyOutput('CPMPlotly')
      ))
    )
    
  }else if(input$displayedGene == "Data Table"){
    removeUI(
      selector = "#CPMTablediv"
    )
    insertUI(
      selector = '#GenePlotDiv',
      where = "beforeBegin",
      ui = tagList(div(id = "CPMTablediv",
                       h5("Log2 Counts per Million (CPM) Expression Across Life Stages"),
                       DTOutput('CPM.datatable')
      ))
    )
    removeUI(
      selector = "#CPMdiv"
    )
    
    removeUI(
      selector = "#CPMPlotlydiv"
    )
    
  }else{
    removeUI(
      selector = "#CPMdiv"
    )
    insertUI(
      selector = '#GenePlotDiv',
      where = "beforeBegin",
      ui = tagList(div(id = "CPMdiv",
                       tabsetPanel(id = "genePlotTab",
                                   type = "pills",
                                   tabPanel("Selected Gene", 
                                            plotOutput('CPM')),
                                   tabPanel(p(em("Caenorhabditis"), "Homologs"), 
                                            plotOutput('CPM.homologs'),
                                            em('Note: homologous genes are identified based on WormBase annotations and may not be accurate.
                                               In addition, C. elegans and C. briggsae data are variance-stabilized log2CPM values, other species\' data is not variance-stabilized.'))
                       )
                       
      ))
    )
    removeUI(
      selector = "#CPMPlotlydiv"
    )
    
    removeUI(
      selector = "#CPMTablediv"
    )
  }
})

## GW: Save Gene Plots ----
output$downloadGenePlot <- downloadHandler(
  
  filename = function(){
    if (input$displayedGene == "All Genes" | input$displayedGene == "Data Table") {
      type <- 'GeneExpression_'
    } else {
    type <- switch(input$genePlotTab,
           `Caenorhabditis Homologs` = "HomologExpression_",
          `Selected Gene` = 'GeneExpression_'
           )}
    paste(type,input$displayedGene, '_',Sys.Date(),'.pdf', sep='')
  },
  content = function(file){
    withProgress ({
      
      if (input$displayedGene == "All Genes" | input$displayedGene == "Data Table") {
        setProgress(0.2)
        # Make a heatmap for all the genes using the Log2CPM values
        myheatcolors <- RdBu(75)
        diffGenes <- vals$diffGenes.df %>%
          dplyr::select(!geneID)
        colnames(diffGenes) <- vals$DGEList.filtered.norm$target$group
        diffGenes <- diffGenes %>%
          as.matrix()
        
        rownames(diffGenes) <- rownames(vals$DGEList.filtered.norm$E)
        subset.diffGenes<- diffGenes[vals$gene_of_interest,]
        
        clustColumns <- hclust(as.dist(1-cor(subset.diffGenes, method="spearman")), method="complete")
        
        colnames(subset.diffGenes) <- paste0(vals$DGEList.filtered.norm$targets$group,
                                             "...",
                                             substr(vals$DGEList.filtered.norm$targets$samples,
                                                    nchar(
                                                      as.character(vals$DGEList.filtered.norm$targets$samples[1]))-2, nchar(
                                                        as.character(vals$DGEList.filtered.norm$targets$samples[1])))
        )
        
        setProgress(0.4)
        
        clustRows <- hclust(as.dist(1-cor(t(subset.diffGenes),
                                          method="pearson")),
                            method="complete")
        
        par(cex.main=1.2)
        setProgress(0.6)
        
        showticklabels <- if(length(vals$gene_of_interest)<20){c(TRUE,TRUE)} else {c(TRUE,FALSE)}
        
        p<-ggheatmap_local(subset.diffGenes,
                           colors = rev(myheatcolors),
                           Rowv= ladderize(as.dendrogram(clustRows)),
                           Colv=ladderize(as.dendrogram(clustColumns)),
                           show_dendrogram = c(TRUE, TRUE),
                           key.title = "Row Z Score",
                           branches_lwd = 0.5,
                           showticklabels = showticklabels,
                           scale='row',
                           cexRow=1.2, cexCol=1.2,
                           main = "Log2 Counts per Million (CPM) Expression Across Life Stages")
        setProgress(0.8)
        ggsave(file,
               plot = p,
               width = 11,
               height = 5,
               device = cairo_pdf)
        setProgress(1)
      } else {
        if (input$genePlotTab == "Selected Gene"){
          p<-vals$gene_plot
          ggsave(file,
                 plot = p,
                 width = 5,
                 height = 4,
                 device = cairo_pdf)
        } else if (input$genePlotTab == "Caenorhabditis Homologs"){
          p<-vals$homolog_plot
          ggsave(file,
                 plot = p,
                 width = 9,
                 height = 4,
                 device = cairo_pdf)
        }
        
      }
    },
    message = "Saving Plots")
  }
  
)

## GW: Save Excel Table with Gene Expression Data ----
output$downloadbuttonsGenes <- renderUI({
  req(vals$genelist,vals$genelist.Log2CPM,input$displayedGene)
  if (nrow(vals$genelist)>1) {req(vals$HeatmapRowOrder)}
  if ((input$displayedGene != "All Genes" | input$displayedGene != "Data Table") & (isTruthy(input$genePlotTab) && input$genePlotTab == "Caenorhabditis Homologs")){
    save.tbl <- fetch_homologs()
    
    #Identify which set of homologs to download the data for
    set_displayed <- dplyr::filter(vals$homologous_genes, 
                                   geneID %in% input$displayedGene)%>%
      as.character()

    genelist.expression <- suppressMessages(save.tbl %>%
      dplyr::filter(geneID %in% set_displayed) %>%
      dplyr::group_by(id, geneID, life_stage) %>%
      dplyr::summarise(avg = round(median(log2CPM),2), 
                       low_hinge = round(fivenum(log2CPM)[2],2), 
                       up_hinge = round(fivenum(log2CPM)[4],2), 
                       .groups = "drop_last") %>%
      tidyr::unite("IQR", low_hinge, up_hinge, sep = " to ") %>%
      tidyr::unite("output", avg, IQR, sep = ", IQR = ") %>%
      dplyr::select(-id) %>%
      pivot_wider(id_cols = life_stage,
                  names_from = c(geneID),
                  names_sep = "-",
                  values_from = output))
    vals$expressionnotes <-"Data values are median with IQR"
    genelist.expression <- list(genelist.expression)
    
  } else {
    vals$genelist.Log2CPM$sampleID <- rep(as.character(vals$DGEList.filtered.norm$targets$samples),
                                          times =  nrow(vals$genelist))
    
    genelist.expression <-  vals$genelist.Log2CPM %>%
      left_join(vals$genelist, .,by = "geneID") %>%
      pivot_wider(id_cols = geneID,
                  names_from = c(life_stage,sampleID),
                  names_sep = "-",
                  values_from = log2CPM)
    
    if (isTruthy(vals$HeatmapRowOrder)){
      row.order <- tibble(OldOrder = vals$HeatmapRowOrder,
                          NewOrder = seq_along(vals$HeatmapRowOrder)) %>%
        dplyr::arrange(OldOrder)
      
      col.order <- tibble(OldOrder = vals$HeatmapColOrder,
                          NewOrder = seq_along(vals$HeatmapColOrder)) %>%
        dplyr::arrange(desc(NewOrder))
      
      genelist.expression <- genelist.expression %>%
        add_column(NewOrder = row.order$NewOrder, .before = "geneID") %>%
        dplyr::arrange(NewOrder) %>%
        dplyr::select(!NewOrder) %>%
        .[c(1, col.order$OldOrder+1)]
    }
    genelist.expression <- genelist.expression %>%
      list("User-selected Genes" = . )
    
    # Add back on genes that were submitted by the user but don't appear in the list of genes for which there is available data.
    excluded.genes <- dplyr::anti_join(vals$submitted.genelist,
                                       vals$genelist,
                                       by = "geneID")
    
    genelist.expression <-lapply(genelist.expression, function (x) {
      dplyr::full_join(x,excluded.genes, by = "geneID")
    })
    
    vals$expressionnotes <- "Data are log2 counts per million"
    # Add gene annotations
    genelist.expression <-lapply(genelist.expression, function (x) {
      vals$annotations %>%
        dplyr::left_join(x,., by = "geneID") %>%
        dplyr::relocate(geneName, .after = geneID)
    })
  }
  output$heatmap_data_download <- generate_excel_report(c("User-selected Genes"),
                                                        genelist.expression,
                                                        name = paste(input$selectSpecies_GW, "RNA-seq Gene Expression"),
                                                        filename_prefix = "Gene_Expression_Data_",
                                                        subtitle_prefix = "Log2CPM Expression:",
                                                        expressionnotes = vals$expressionnotes,
                                                        n_header_rows = 6)
  
  tagList(
    downloadButton("downloadGenePlot",
                   "Download Plot as PDF",
                   class = "btn-primary"),
    
    downloadButton("heatmap_data_download",
                   "Download Data as Excel",
                   class = "btn-primary")
  )
  
})