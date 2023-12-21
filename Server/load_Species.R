# This script loads species-specific data for the Caenorhabditis RNAseq Browser

# GW: Load species data in Gene-wise tab ----
observeEvent(input$speciesGW, {
    input$selectSpecies_GW
    vals$genelist <- NULL
    vals$HeatmapRowOrder <- NULL
    vals$target.contrast.options <- NULL
    updateTextAreaInput(session,"idtext",value = "")
    updatePickerInput(session, "displayedGene", choices = "", selected = "")
    removeUI(selector = "#CPMPlotlydiv_parent")
    removeUI(selector = "#geneDisplaySelectionPanel")
    updateSelectInput(session, "selectContrast_GW", selected = "")
    updateSelectInput(session, "selectTarget_GW", selected = "")
    updateTextAreaInput(session,"multiContrasts_GW",value = "")
    vals$comparison_GW <- NULL
    
    species <- switch(input$selectSpecies_GW,
                      `C. elegans` = 'elegans',
                      `C. briggsae` = 'briggsae',
                      `C. brenneri` = 'brenneri',
                      `C. japonica` = "japonica",
                      `C. remanei` = "remanei")

    
    withProgress({
        
        # Import a variance-stabilized DGEList created by voom transformation command.
        # Outputs: E = normalized CPMexpression values on the log2 scale
        load(file = paste0("./Data/",species,"_vDGEList"))
        vals$v.DGEList.filtered.norm <- v.DGEList.filtered.norm
        
        setProgress(value = .25)
        
        vals$target.contrast.options <- vals$v.DGEList.filtered.norm$targets$group
        # Import a tidy dataframe containing gene annotations for all genes in the genome (including those that are excluded from this database.)
        
        load(file = paste0("./Data/",species,"_geneAnnotations"))
        vals$annotations <- annotations
        
        setProgress(value = .5)
        
        # Parse vDGEList into a tibble containing Log2CPM information
        vals$Log2CPM<-v.DGEList.filtered.norm$E %>%
            as_tibble(rownames = "geneID")%>%
            setNames(nm = c("geneID", 
                            as.character(v.DGEList.filtered.norm$targets$group))) %>%
            pivot_longer(cols = -geneID,
                         names_to = "life_stage", 
                         values_to = "log2CPM") %>%
            group_by(geneID, life_stage)
        
        vals$diffGenes.df <- v.DGEList.filtered.norm$E %>%
            as_tibble(rownames = "geneID", .name_repair = "unique")
        
        setProgress(value = .75)
        
        ## Fit a linear model to the data
        vals$fit <- lmFit(v.DGEList.filtered.norm, v.DGEList.filtered.norm$design)
        
        setProgress(value = 1)
        
    }, message = "Loading Species Database")
    
    
})

# GW: Load experiment information ----
StudyInfo.filename.GW <- reactive({
    req(input$selectSpecies_GW)
    
    species <- switch(input$selectSpecies_GW,
                      `C. elegans` = 'elegans',
                      `C. briggsae` = 'briggsae',
                      `C. brenneri` = 'brenneri',
                      `C. japonica` = "japonica",
                      `C. remanei` = "remanei")
    
    Info.type <- switch(input$which.Experimental.Info.GW,
                        `Study Design` = '_study_design.txt',
                        `Log2CPM Gene Counts` = '_log2cpm_filtered_norm_voom.csv',
                        `vDGEList` = "_vDGEList")
    
    file.location <- switch(input$which.Experimental.Info.GW,
                            `Study Design` = './www/',
                            `Log2CPM Gene Counts` = './www/',
                            `vDGEList` = "./Data/"
    )
    Info.file <- paste0(file.location,species, Info.type)
    Info.file
    
})

output$StudyInfo.panel.GW <- renderUI({
    output$StudyInfo.file.GW <- downloadHandler(
        filename = function() {
            Info.file <- StudyInfo.filename.GW()
            file.location <- switch(input$which.Experimental.Info.GW,
                                    `Study Design` = './www/',
                                    `Log2CPM Gene Counts` = './www/',
                                    `vDGEList` = "./Data/"
            )
            str_remove(Info.file, file.location)
        },
        content = function(file){
            Info.file <- StudyInfo.filename.GW()
            file.copy(Info.file, file)
        }
    )
    
    downloadButton("StudyInfo.file.GW","Download",
                   class = "btn-default")
})



# LS: Load species data in life stage tab ----
observeEvent(input$speciesLS, {
    input$selectSpecies_LS
    
    updateSelectInput(session, "selectContrast_LS", selected = "")
    updateSelectInput(session, "selectTarget_LS", selected = "")
    updateTextAreaInput(session,"multiContrasts_LS",value = "")
    
    species <- switch(input$selectSpecies_LS,
                      `C. elegans` = 'elegans',
                      `C. briggsae` = 'briggsae',
                      `C. brenneri` = 'brenneri',
                      `C. japonica` = "japonica",
                      `C. remanei` = "remanei")
    
    withProgress({
        # Import a variance-stabilized DGEList created by voom transformation command.
        # Outputs: E = normalized CPMexpression values on the log2 scale
        load(file = paste0("./Data/",species,"_vDGEList"))
        vals$v.DGEList.filtered.norm <- v.DGEList.filtered.norm
        
        setProgress(value = .25)
        
        # Import a tidy dataframe containing gene annotations for all genes in the genome (including those that are excluded from this database.)
        load(file = paste0("./Data/",species,"_geneAnnotations"))
        vals$annotations <- annotations
        
        setProgress(value = .5)
       
        # Parse vDGEList into a tibble containing Log2CPM information
        vals$Log2CPM<-v.DGEList.filtered.norm$E %>%
            as_tibble(rownames = "geneID")%>%
            setNames(nm = c("geneID", 
                            as.character(v.DGEList.filtered.norm$targets$group))) %>%
            pivot_longer(cols = -geneID,
                         names_to = "life_stage", 
                         values_to = "log2CPM") %>%
            group_by(geneID, life_stage)
        
        vals$diffGenes.df <- v.DGEList.filtered.norm$E %>%
            as_tibble(rownames = "geneID", .name_repair = "unique")
        
        setProgress(value = .75)
        
        ## Fit a linear model to the data
        vals$fit <- lmFit(v.DGEList.filtered.norm, v.DGEList.filtered.norm$design)
        
        setProgress(value = 1)
    }, message = "Loading Species Database")
    
    
})

# LS: Load experiment information ----
StudyInfo.filename.LS <- reactive({
    req(input$selectSpecies_LS)
    
    species <- switch(input$selectSpecies_LS,
                      `C. elegans` = 'elegans',
                      `C. briggsae` = 'briggsae',
                      `C. brenneri` = 'brenneri',
                      `C. japonica` = "japonica",
                      `C. remanei` = "remanei")
    
    Info.type <- switch(input$which.Experimental.Info.LS,
                        `Study Design` = '_study_design.txt',
                        `Log2CPM Gene Counts` = '_log2cpm_filtered_norm_voom.csv',
                        `vDGEList` = "_vDGEList")
    
    file.location <- switch(input$which.Experimental.Info.LS,
                            `Study Design` = './www/',
                            `Log2CPM Gene Counts` = './www/',
                            `vDGEList` = "./Data/"
    )
    Info.file <- paste0(file.location,species, Info.type)
    Info.file
    
})

output$StudyInfo.panel.LS <- renderUI({
    output$StudyInfo.file.LS <- downloadHandler(
        filename = function() {
            Info.file <- StudyInfo.filename.LS()
            file.location <- switch(input$which.Experimental.Info.LS,
                                    `Study Design` = './www/',
                                    `Log2CPM Gene Counts` = './www/',
                                    `vDGEList` = "./Data/"
            )
            str_remove(Info.file, file.location)
        },
        content = function(file){
            Info.file <- StudyInfo.filename.LS()
            file.copy(Info.file, file)
        }
    )
    
    downloadButton("StudyInfo.file.LS","Download",
                   class = "btn-default")
})

# About Tab: Download experiment information ----
StudyInfo.filename.About <- reactive({
    Info.type <- switch(input$which.Experimental.Info.About,
                        `C. elegans Study Design` = 'elegans_study_design.txt',
                        `C. elegans Log2CPM Gene Counts` = 'elegans_log2cpm_filtered_norm_voom.csv',
                        `C. elegans vDGEList` = "elegans_vDGEList",
                        `C. briggsae Study Design` = 'briggsae_study_design.txt',
                        `C. briggsae Log2CPM Gene Counts` = 'briggsae_log2cpm_filtered_norm_voom.csv',
                        `C. briggsae vDGEList` = "briggsae_vDGEList",
                        `C. brenneri Study Design` = 'brenneri_study_design.txt',
                        `C. brenneri Log2CPM Gene Counts` = 'brenneri_log2cpm_filtered_norm_voom.csv',
                        `C. brenneri vDGEList` = "brenneri_vDGEList",
                        `C. japonica Study Design` = 'japonica_study_design.txt',
                        `C. japonica Log2CPM Gene Counts` = 'japonica_log2cpm_filtered_norm_voom.csv',
                        `C. japonica vDGEList` = "japonica_vDGEList",
                        `C. remanei Study Design` = 'remanei_study_design.txt',
                        `C. remanei Log2CPM Gene Counts` = 'remanei_log2cpm_filtered_norm_voom.csv',
                        `C. remanei vDGEList` = "remanei_vDGEList"
                        )
                        
    
    file.location <- switch(input$which.Experimental.Info.About,
                            `C. elegans Study Design` = './www/',
                            `C. elegans Log2CPM Gene Counts` = './www/',
                            `C. elegans vDGEList` = "./Data/",
                            `C. briggsae Study Design` = './www/',
                            `C. briggsae Log2CPM Gene Counts` = './www/',
                            `C. briggsae vDGEList` = "./Data/",
                            `C. brenneri Study Design` = './www/',
                            `C. brenneri Log2CPM Gene Counts` = './www/',
                            `C. brenneri vDGEList` = "./Data/",
                            `C. japonica Study Design` = './www/',
                            `C. japonica Log2CPM Gene Counts` = './www/',
                            `C. japonica vDGEList` = "./Data/",
                            `C. remanei Study Design` = './www/',
                            `C. remanei Log2CPM Gene Counts` = './www/',
                            `C. remanei vDGEList` = "./Data/" 
    )
    Info.file <- paste0(file.location, Info.type)
    Info.file
    
})

output$StudyInfo.panel.About <- renderUI({
    output$StudyInfo.file.About <- downloadHandler(
        filename = function() {
            Info.file <- StudyInfo.filename.About()
            file.location <- switch(input$which.Experimental.Info.About,
                                    `C. elegans Study Design` = './www/',
                                    `C. elegans Log2CPM Gene Counts` = './www/',
                                    `C. elegans vDGEList` = "./Data/",
                                    `C. briggsae Study Design` = './www/',
                                    `C. briggsae Log2CPM Gene Counts` = './www/',
                                    `C. briggsae vDGEList` = "./Data/",
                                    `C. brenneri Study Design` = './www/',
                                    `C. brenneri Log2CPM Gene Counts` = './www/',
                                    `C. brenneri vDGEList` = "./Data/",
                                    `C. japonica Study Design` = './www/',
                                    `C. japonica Log2CPM Gene Counts` = './www/',
                                    `C. japonica vDGEList` = "./Data/",
                                    `C. remanei Study Design` = './www/',
                                    `C. remanei Log2CPM Gene Counts` = './www/',
                                    `C. remanei vDGEList` = "./Data/"   
            )
            str_remove(Info.file, file.location)
        },
        content = function(file){
            Info.file <- StudyInfo.filename.About()
            file.copy(Info.file, file)
        }
    )
    
    downloadButton("StudyInfo.file.About","Download",
                   class = "btn-primary")
})