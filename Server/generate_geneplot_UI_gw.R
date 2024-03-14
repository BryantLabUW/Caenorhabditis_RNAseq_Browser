## GW: Generate Panel for Gene-wise Plots ----
output$genePlotPanel_GW <- renderUI({
    parse_ids()
    req(vals$genelist)
    isolate({
        tagList(div(id = "CPMPlotlydiv_parent",
                    panel(
                        heading = tagList(h5(shiny::icon("fas fa-chart-bar"),
                                             "Gene Expression Across Life Stages")),
                        status = "primary",
                        
                        tagList(div(id = "GenePlotDiv",
                                    uiOutput("downloadbuttonsGenes")
                        ))
                    )
        ))
        
    })
    
})

## GW: Generate Responsive Selection for Gene to Display ----
output$geneDisplaySelection_GW <- renderUI({
  parse_ids()
  req(vals$genelist)
  vals$genelist <- vals$genelist %>%
    left_join(vals$annotations, by = "geneID") %>%
    dplyr::mutate(geneReference = paste0(geneID, " | ", geneName)) %>%
    dplyr::select(geneID, geneReference)
  
    if (length(vals$genelist$geneID)>1) {
        choices <- c("All Genes", "Data Table",vals$genelist$geneReference)
    } else {
        choices <- c(vals$genelist$geneReference,"Data Table")}
    tagList(div(id = "geneDisplaySelectionPanel",
                panel(
                    heading = tagList(h5(shiny::icon("fas fa-filter"),
                                         "Pick Gene to Display")),
                    status = "warning",
                    pickerInput("displayedGene",
                                NULL, 
                                choices,
                                options = list(style = 'btn btn-warning'))
                )
    ))
})

## GW: Generate Legend Explaining the Life Stages ----
output$lifeStageLegend_GW <- renderDT({
    lifestage_legend.df <- lifestage_legend %>%
        dplyr::select(any_of(unique(vals$DGEList.filtered.norm$targets$group))) %>%
        as.data.frame()
    rownames(lifestage_legend.df)<- c("<b>Life Stage</b>")
    lifestage_legend.DT <- lifestage_legend.df %>%
        DT::datatable(rownames = TRUE,
                      escape = FALSE,
                      class = "table-bordered",
                      colnames = c("Abbr." = 1),
                      options = list(scrollX = TRUE,
                                     ordering = FALSE,
                                     dom = 'tS'
                      )
        )
    lifestage_legend.DT
})

output$Legend_GW <- renderUI({
    req(vals$genelist)
    vals$genelist
    panel(
        heading = tagList(h5(shiny::icon("fas fa-book-open"),
                             "Sample ID Legend")),
        status = "info",
        DTOutput("lifeStageLegend_GW")
    )
})
