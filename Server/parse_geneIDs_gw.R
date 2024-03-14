## GW: Generate Main Gene Input Panel ----
output$genePanelinputs <- renderUI({
    tagList(
        panel(
            id = "GeneInputBox",
            heading = tagList(h5(shiny::icon("fas fa-dna"), "Step 1: Input Genes / Keywords")),
            status = "primary",
            ### GeneID (text box)
            h5('Pick Genes', class = 'text-info', style = "margin: 0px 0px 5px 0px"),
            p(tags$em('Users may type gene names, stable IDs, or keywords that will be matched against InterPro gene descriptions. Search terms may be separated using commas, semicolons, or new lines. Users may also upload a .csv file containing search terms.', style = "color: #7b8a8b")),
            p(tags$em("Type 'everything' or 'all genes' to display all genes in the genome. Warning: this will take a long time to process.", style = "color: #7b8a8b")),
            p(tags$em(tags$b('Note: Please hit the Clear button if switching between typing and uploading inputs.', style = "color: #bf9232"))),
            textAreaInput('idtext',
                          h6('Gene IDs or Keywords'),
                          rows = 5, 
                          resize = "vertical"),
            
            ### Upload list of GeneIDs
            uiOutput('genefile_upload'),
            
            ### Action Button
            actionButton('goGW',
                         'Submit',
                         icon = icon("fas fa-share"),
                         class = "btn-primary"),
            
            actionButton('resetGenes', 'Clear',
                         icon = icon("far fa-trash-alt"))
        )
    )
})

## GW: Parse Gene Inputs ----
parse_ids <- eventReactive(input$goGW,{
    vals$genelist <- NULL
    vals$HeatmapRowOrder <- NULL
    validate(
        need(isTruthy(vals$DGEList.filtered.norm), 
             "Please re-select a species for analysis")
    )
    
    validate(
        need({isTruthy(input$loadfile) | isTruthy(input$idtext)}, "Please input gene ids")
    )
    isolate({
        withProgress({
          
        if (isTruthy(input$idtext)){
            terms <- input$idtext %>%
                gsub("\\n",",",.) %>% #replace any new lines with commas
                trimWhiteSpace %>% #remove leading and trailing white space from string
                str_split(pattern = ",|;") %>%
                unlist()
        } else if (isTruthy(input$loadfile)){
            file <- input$loadfile
            ext <- tools::file_ext(file$datapath)
            validate(need(ext == "csv", "Please upload a csv file"))
            suppressWarnings(
                terms <- read.csv(file$datapath, 
                                  header = FALSE, 
                                  colClasses = "character", 
                                  strip.white = T,
                                  na.strings = "") %>%
                    as_tibble() %>%
                    
                    pivot_longer(cols = everything(), 
                                 values_to = "geneID",
                                 values_drop_na = T) 
            )
            terms <- terms$geneID
            
        } 
        setProgress(0.1)
        
            genelist <- vals$annotations %>%
            dplyr::select(geneID) %>%
              dplyr::ungroup()
        
        if (any(grepl('everything|all genes', terms, ignore.case = TRUE))) {
            # Text input matches the strings 'everything' or 'all genes'
            genelist <- genelist
        } else {
            inc <- 0.1/nrow(terms)
            # Search for gene IDs
            terms.cleaned <- gsub("\\.[0-9]$","",terms) #strip any transcript values from the inputted list
            geneindex.geneID<-sapply(terms.cleaned, function(y) {
                incProgress(amount = inc)
                yy <- gsub("^\\s+|\\s+$", "", y) #remove any number of whitespace from start or end
                grepl(paste0("\\<",yy,"\\>"), 
                      vals$annotations$geneID,
                      ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            geneindex.geneNames<-sapply(terms.cleaned, function(y) {
              incProgress(amount = inc)
              yy <- gsub("^\\s+|\\s+$", "", y) #remove any number of whitespace from start or end
              grepl(paste0("\\<",yy,"\\>"),
                    vals$annotations$geneName,
                    ignore.case = TRUE)
            }) %>%
              rowSums() %>%
              as.logical()
            
            if (input$selectSpecies_GW != "C. elegans"){
            geneindex.GS1<-sapply(terms.cleaned, function(y) {
              incProgress(amount = inc)
              yy <- gsub("^\\s+|\\s+$", "", y) #remove any number of whitespace from start or end
              grepl(paste0("\\<",yy,"\\>"),
                    vals$annotations$GS1_homologID,
                    ignore.case = TRUE)
            }) %>%
              rowSums() %>%
              as.logical()
            } else
              geneindex.GS1 <-logical(1)

           # Search WormBase Gene Description Terms
            # geneindex.description<-sapply(terms.cleaned, function(y) {
            #     incProgress(amount = inc)
            #   yy <- gsub("^\\s+|\\s+$", "", y) #remove any number of whitespace from start or end
            #   grepl(paste0("\\<",yy,"\\>"),
            #           vals$annotations$Description,
            #           ignore.case = TRUE)
            # }) %>%
            #     rowSums() %>%
            #     as.logical()
           
             # Search InterPro Terms
            geneindex.interpro<-sapply(terms.cleaned, function(y) {
              incProgress(amount = inc)
              grepl(gsub("^\\s+|\\s+$", "", y), #remove any number of whitespace from start or end
                    vals$annotations$InterPro,
                    ignore.case = TRUE)
            }) %>%
              rowSums() %>%
              as.logical()
            
            # Search WormBaseID Terms
            geneindex.wbid<-sapply(terms.cleaned, function(y) {
              incProgress(amount = inc)
              yy <- gsub("^\\s+|\\s+$", "", y) #remove any number of whitespace from start or end
              grepl(paste0("\\<",yy,"\\>"),
                    vals$annotations$WormBaseID,
                    ignore.case = TRUE)
            }) %>%
              rowSums() %>%
              as.logical()
            
            geneindex <- geneindex.geneID | geneindex.geneNames | geneindex.wbid | geneindex.interpro | geneindex.GS1 
            genelist <- dplyr::filter(genelist,geneindex) 
        }
       
        if (nrow(genelist) == 0){
            disable("goLifeStage_GW")
        } else {enable("goLifeStage_GW")}
        
        # Produces error message if genelist is empty
        validate(
            need(nrow(genelist) != 0, "RNA-seq data unavailable for submitted genes (submitted names may be invalid). Please try a new search.")
        )
       
        # Save record of original genelist before filtering, removing rows that contain the word 'gene'
        vals$submitted.genelist <- genelist %>%
            dplyr::filter(!grepl("gene", geneID, ignore.case = T))
        
        # Remove genes from the list that aren't part of vals$Log2CPM
        genelist <- genelist %>%
            dplyr::filter(geneID %in% vals$Log2CPM$geneID)
        
        setProgress(0.95)
        
        if (nrow(genelist) == 0){
            disable("goLifeStage_GW")
        } else {enable("goLifeStage_GW")}
        
        # Produces error message if genelist is empty after removing genes not included in the RNAseq dataset
        validate(
            need(nrow(genelist) != 0, "RNA-seq data unavailable for submitted genes. Please try a new search.")
        )
        
        vals$genelist <- genelist
        vals$genelist.Log2CPM <- vals$Log2CPM %>%
            dplyr::filter(geneID %in% genelist$geneID)
        },message = "Parsing gene IDs...")
    })
})