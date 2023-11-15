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
                      `S. stercoralis` = 'Ss',
                      `C. elegans` = 'Cele',
                      `C. briggsae` = 'Cbri',
                      `C. brenneri` = 'Cbre',
                      `C. japonica` = "Cjap",
                      `C. remanei` = "Crem")

    
    withProgress({
        
        # Import a variance-stabilized DGEList created by voom transformation command.
        # Outputs: E = normalized CPMexpression values on the log2 scale
        load(file = paste0("./Data/",species,"_vDGEList"))
        vals$v.DEGList.filtered.norm <- v.DEGList.filtered.norm
        
        setProgress(value = .25)
        
        vals$target.contrast.options <- vals$v.DEGList.filtered.norm$targets$group
        # Import a tidy dataframe containing gene annotations for all genes in the genome (including those that are excluded from this database.)
        load(file = paste0("./Data/",species,"_geneAnnotations"))
        vals$annotations <- as_tibble(annotations, rownames = "geneID")
        
        setProgress(value = .5)
        
        # Parse vDGEList into a tibble containing Log2CPM information
        vals$Log2CPM<-v.DEGList.filtered.norm$E %>%
            as_tibble(rownames = "geneID")%>%
            setNames(nm = c("geneID", 
                            as.character(v.DEGList.filtered.norm$targets$group))) %>%
            pivot_longer(cols = -geneID,
                         names_to = "life_stage", 
                         values_to = "log2CPM") %>%
            group_by(geneID, life_stage)
        
        vals$diffGenes.df <- v.DEGList.filtered.norm$E %>%
            as_tibble(rownames = "geneID", .name_repair = "unique")
        
        setProgress(value = .75)
        
        ## Fit a linear model to the data
        vals$fit <- lmFit(v.DEGList.filtered.norm, v.DEGList.filtered.norm$design)
        
        setProgress(value = 1)
        
    }, message = "Loading Species Database")
    
    
})

# GW: Load experiment information ----
StudyInfo.filename.GW <- reactive({
    req(input$selectSpecies_GW)
    
    species <- switch(input$selectSpecies_GW,
                      `S. stercoralis` = 'Ss',
                      `C. elegans` = 'Cele',
                      `C. briggsae` = 'Cbri',
                      `C. brenneri` = 'Cbre',
                      `C. japonica` = "Cjap",
                      `C. remanei` = "Crem")
    
    Info.type <- switch(input$which.Experimental.Info.GW,
                        `Study Design` = '_studyDesign.txt',
                        `Log2CPM Gene Counts` = 'RNAseq_log2cpm_filtered_norm_voom.csv',
                        `vDGEList` = "_vDGEList",
                        `Discarded Gene Counts` = "RNAseq_discardedGene_counts.csv")
    
    file.location <- switch(input$which.Experimental.Info.GW,
                            `Study Design` = './www/',
                            `Log2CPM Gene Counts` = './www/',
                            `vDGEList` = "./Data/",
                            `Discarded Gene Counts` = "./www/"
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
                                    `vDGEList` = "./Data/",
                                    `Discarded Gene Counts` = "./www/"
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
                      `S. stercoralis` = 'Ss',
                      `C. elegans` = 'Cele',
                      `C. briggsae` = 'Cbri',
                      `C. brenneri` = 'Cbre',
                      `C. japonica` = "Cjap",
                      `C. remanei` = "Crem")
    
    withProgress({
        # Import a variance-stabilized DEGList created by voom transformation command.
        # Outputs: E = normalized CPMexpression values on the log2 scale
        load(file = paste0("./Data/",species,"_vDGEList"))
        vals$v.DEGList.filtered.norm <- v.DEGList.filtered.norm
        
        setProgress(value = .25)
        
        # Import a tidy dataframe containing gene annotations for all genes in the genome (including those that are excluded from this database.)
        load(file = paste0("./Data/",species,"_geneAnnotations"))
        vals$annotations <- as_tibble(annotations, rownames = "geneID")
        
        setProgress(value = .5)
        
        # Parse vDGEList into a tibble containing Log2CPM information
        vals$Log2CPM<-v.DEGList.filtered.norm$E %>%
            as_tibble(rownames = "geneID")%>%
            setNames(nm = c("geneID", 
                            as.character(v.DEGList.filtered.norm$targets$group))) %>%
            pivot_longer(cols = -geneID,
                         names_to = "life_stage", 
                         values_to = "log2CPM") %>%
            group_by(geneID, life_stage)
        
        vals$diffGenes.df <- v.DEGList.filtered.norm$E %>%
            as_tibble(rownames = "geneID", .name_repair = "unique")
        
        setProgress(value = .75)
        
        ## Fit a linear model to the data
        vals$fit <- lmFit(v.DEGList.filtered.norm, v.DEGList.filtered.norm$design)
        
        setProgress(value = 1)
    }, message = "Loading Species Database")
    
    
})

# LS: Load experiment information ----
StudyInfo.filename.LS <- reactive({
    req(input$selectSpecies_LS)
    
    species <- switch(input$selectSpecies_LS,
                      `S. stercoralis` = 'Ss',
                      `C. elegans` = 'Cele',
                      `C. briggsae` = 'Cbri',
                      `C. brenneri` = 'Cbre',
                      `C. japonica` = "Cjap",
                      `C. remanei` = "Crem")
    
    Info.type <- switch(input$which.Experimental.Info.LS,
                        `Study Design` = '_studyDesign.txt',
                        `Log2CPM Gene Counts` = 'RNAseq_log2cpm_filtered_norm_voom.csv',
                        `vDGEList` = "_vDGEList",
                        `Discarded Gene Counts` = "RNAseq_discardedGene_counts.csv")
    
    file.location <- switch(input$which.Experimental.Info.LS,
                            `Study Design` = './www/',
                            `Log2CPM Gene Counts` = './www/',
                            `vDGEList` = "./Data/",
                            `Discarded Gene Counts` = "./www/"
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
                                    `vDGEList` = "./Data/",
                                    `Discarded Gene Counts` = "./www/"
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
                        `Ss Study Design` = 'Ss_studyDesign.txt',
                        `Ss Log2CPM Gene Counts` = 'SsRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Ss vDGEList` = "Ss_vDGEList",
                        `Ss Discarded Gene Counts` = "SsRNAseq_discardedGene_counts.csv",
                        `Cele Study Design` = 'Cele_studyDesign.txt',
                        `Cele Log2CPM Gene Counts` = 'CeleRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Cele vDGEList` = "Cele_vDGEList",
                        `Cele Discarded Gene Counts` = "CeleRNAseq_discardedGene_counts.csv",
                        `Cbri Study Design` = 'Cbri_studyDesign.txt',
                        `Cbri Log2CPM Gene Counts` = 'CbriRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Cbri vDGEList` = "Cbri_vDGEList",
                        `Cbri Discarded Gene Counts` = "CbriRNAseq_discardedGene_counts.csv",
                        `Cbre Study Design` = 'Cbre_studyDesign.txt',
                        `Cbre Log2CPM Gene Counts` = 'CbreRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Cbre vDGEList` = "Cbre_vDGEList",
                        `Cbre Discarded Gene Counts` = "CbreRNAseq_discardedGene_counts.csv",
                        `Cjap Study Design` = 'Cjap_studyDesign.txt',
                        `Cjap Log2CPM Gene Counts` = 'CjapRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Cjap vDGEList` = "Cjap_vDGEList",
                        `Cjap Discarded Gene Counts` = "CjapRNAseq_discardedGene_counts.csv",
                        `Crem Study Design` = 'Crem_studyDesign.txt',
                        `Crem Log2CPM Gene Counts` = 'CremRNAseq_log2cpm_filtered_norm_voom.csv',
                        `Crem vDGEList` = "Crem_vDGEList",
                        `Crem Discarded Gene Counts` = "CremRNAseq_discardedGene_counts.csv")
    
    file.location <- switch(input$which.Experimental.Info.About,
                            `Ss Study Design` = './www/',
                            `Ss Log2CPM Gene Counts` = './www/',
                            `Ss vDGEList` = "./Data/",
                            `Ss Discarded Gene Counts` = "./www/",
                            `Cele Study Design` = './www/',
                            `Cele Log2CPM Gene Counts` = './www/',
                            `Cele vDGEList` = "./Data/",
                            `Cele Discarded Gene Counts` = "./www/",
                            `Cbri Study Design` = './www/',
                            `Cbri Log2CPM Gene Counts` = './www/',
                            `Cbri vDGEList` = "./Data/",
                            `Cbri Discarded Gene Counts` = "./www/",
                            `Cbre Study Design` = './www/',
                            `Cbre Log2CPM Gene Counts` = './www/',
                            `Cbre vDGEList` = "./Data/",
                            `Cbre Discarded Gene Counts` = "./www/",
                            `Cjap Study Design` = './www/',
                            `Cjap Log2CPM Gene Counts` = './www/',
                            `Cjap vDGEList` = "./Data/",
                            `Cjap Discarded Gene Counts` = "./www/",
                            `Crem Study Design` = './www/',
                            `Crem Log2CPM Gene Counts` = './www/',
                            `Crem vDGEList` = "./Data/",
                            `Crem Discarded Gene Counts` = "./www/"  
    )
    Info.file <- paste0(file.location, Info.type)
    Info.file
    
})

output$StudyInfo.panel.About <- renderUI({
    output$StudyInfo.file.About <- downloadHandler(
        filename = function() {
            Info.file <- StudyInfo.filename.About()
            file.location <- switch(input$which.Experimental.Info.About,
                                    `Ss Study Design` = './www/',
                                    `Ss Log2CPM Gene Counts` = './www/',
                                    `Ss vDGEList` = "./Data/",
                                    `Ss Discarded Gene Counts` = "./www/",
                                    `Cele Study Design` = './www/',
                                    `Cele Log2CPM Gene Counts` = './www/',
                                    `Cele vDGEList` = "./Data/",
                                    `Cele Discarded Gene Counts` = "./www/",
                                    `Cbri Study Design` = './www/',
                                    `Cbri Log2CPM Gene Counts` = './www/',
                                    `Cbri vDGEList` = "./Data/",
                                    `Cbri Discarded Gene Counts` = "./www/",
                                    `Cbre Study Design` = './www/',
                                    `Cbre Log2CPM Gene Counts` = './www/',
                                    `Cbre vDGEList` = "./Data/",
                                    `Cbre Discarded Gene Counts` = "./www/",
                                    `Cjap Study Design` = './www/',
                                    `Cjap Log2CPM Gene Counts` = './www/',
                                    `Cjap vDGEList` = "./Data/",
                                    `Cjap Discarded Gene Counts` = "./www/",
                                    `Crem Study Design` = './www/',
                                    `Crem Log2CPM Gene Counts` = './www/',
                                    `Crem vDGEList` = "./Data/",
                                    `Crem Discarded Gene Counts` = "./www/"   
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