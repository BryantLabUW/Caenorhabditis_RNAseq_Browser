## GW: Generate/Reset Gene File Upload ----
output$genefile_upload <- renderUI({
    input$resetGenes
    input$speciesGW
    input$tab
    fileInput('loadfile',
              h6('Gene Stable ID List (.csv)'),
              multiple = FALSE)
})


## GW: Clear Genes with Reset Button ----
observeEvent(input$resetGenes, {
    updatePickerInput(session, "displayedComparison_GW", choices = "", selected = "")
    vals$genelist <- NULL
    vals$HeatmapRowOrder <- NULL
    updateTextAreaInput(session,"idtext",value = "")
    updatePickerInput(session, "displayedGene", choices = "", selected = "")
    removeUI(selector = "#CPMPlotlydiv_parent")
    removeUI(selector = "#geneDisplaySelectionPanel")
    updateSelectInput(session, "selectContrast_GW", selected = "")
    updateSelectInput(session, "selectTarget_GW", selected = "")
    updateTextAreaInput(session,"multiContrasts_GW",value = "")
    
    
    vals$comparison_GW <- NULL
    vals$targetStage_GW <- NULL 
    vals$contrastStage_GW <- NULL
    vals$limmacontrast_GW <- NULL
    vals$target.contrast.options <- NULL
    vals$target.contrast.options <- vals$DGEList.filtered.norm$targets$group
})

## Reset Elements on Tab Change ----   
     observeEvent(input$tab, {
         vals$DGEList.filtered.norm <- NULL
         shinyjs::reset("GW")
         shinyjs::reset("LS")
         vals$genelist <- NULL
         vals$HeatmapRowOrder <- NULL
         updateTextAreaInput(session,"idtext",value = "")
         updatePickerInput(session, "displayedGene", choices = "", selected = "")
         removeUI(selector = "#CPMPlotlydiv_parent")
         removeUI(selector = "#geneDisplaySelectionPanel")
         updateSelectInput(session, "selectContrast_GW", selected = "")
         updateSelectInput(session, "selectTarget_GW", selected = "")
         updateTextAreaInput(session,"multiContrasts_GW",value = "")
         vals$comparison_GW <- NULL
         updateSelectInput(session, "selectContrast_LS", selected = "")
         updateSelectInput(session, "selectTarget_LS", selected = "")
         updateTextAreaInput(session,"multiContrasts_LS",value = "")
         vals$target.contrast.options <- NULL
     })
