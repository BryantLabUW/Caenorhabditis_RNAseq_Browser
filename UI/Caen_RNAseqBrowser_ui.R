

# Header ----
navbarPage(title = h3(em("Caenorhabditis"), "RNA-seq Browser"),
           windowTitle = "Caen-RNA-seq Browser",
           theme = shinytheme("united"), 
           collapsible = F,
           id = "tab",
           
           # Gene Browser Tab ----
           tabPanel(title = h4("Browse By Gene"),
                    value = "GW",
                    useShinyjs(),
                    div(id = "GW",
                        
                        ## Fluid Row 1: General usage notes to orient users ----
                        fluidRow(
                            column(9,
                                   alert(status = "success",
                                         dismissible = TRUE,
                                         id = "userNotes_GW",
                                         "Due to server constraints sessions run on shinyapps.io will time out after several minutes of inactivity. If this happens, please reload the server connection."
                                         )),
                            column(3, 
                                   alert(status = "info",
                                         dismissible = TRUE,
                                         id = "bugReport_GW", 
                                         "To report a bug, please email Dr. Astra Bryant at astrab@uw.edu"))
                                   
                        ),
                        
                        ## Fluid Row 2: Species Selection Panel + Download Study Info Dropdown Menu ----    
                        fluidRow(
                            
                            column(3,
                                   
                                   panel(
                                       heading = tagList(h5(shiny::icon("fas fa-archive"),
                                                            "Select Species")),
                                       status = "primary",
                                       id = "speciesPanelID_GW",
                                       selectInput("selectSpecies_GW",
                                                   label = h6("Pick a species"),
                                                   choices = list("C. elegans",
                                                                  "C. briggsae",
                                                                  "C. brenneri",
                                                                  "C. japonica",
                                                                  "C. remanei")
                                       ),
                                       actionButton('speciesGW',
                                                    'Initialize',
                                                    icon = icon("fas fa-share"),
                                                    class = "btn-primary")
                                   )),
                            column(3,
                                   conditionalPanel(condition = "input.speciesGW != 0",
                                                    panel(heading = tagList(h5(shiny::icon("fas fa-info-circle"),
                                                                               "Download Study Info Files")),
                                                          status = "default",
                                                          pickerInput("which.Experimental.Info.GW",
                                                                      NULL, 
                                                                      c("Study Design",
                                                                        "Log2CPM Gene Counts",
                                                                        "vDGEList"),
                                                                      options = list(style = 'btn btn-default')),
                                                          uiOutput("StudyInfo.panel.GW")
                                                    )
                                   ),
                                   offset = 6
                            )
                        ),
                        ## Fluid Row 3: Input Gene Panel + Life Stage Legend + Gene Dropdown Menu + Gene Expression Plots Panel----
                        fluidRow(
                            
                            column(3,
                                   
                                   conditionalPanel(condition = "input.speciesGW != 0",
                                                    uiOutput("genePanelinputs")
                                   )),
                            column(3,
                                   conditionalPanel(condition = "input.goGW !=0",
                                                    id = "geneSelection_conditionalPanel",
                                                    uiOutput("geneDisplaySelection_GW"))
                                   
                            ),
                            column(6,
                                   conditionalPanel(condition = "input.goGW !=0 && output.geneDisplaySelection_GW",
                                                    id = "lifeStageLegend_GW",
                                                    uiOutput("Legend_GW"))
                                   
                            ),
                            
                            column(9,
                                   conditionalPanel(condition = "input.goGW != 0",
                                                    id = "geneplot_conditionalPanel",
                                                    uiOutput("genePlotPanel_GW")
                                                    
                                   )
                            )
                            
                        ),
                        ## Fluid Row 4: Input Pairwise Comparisons Panel + Volcano Plot Panel + Contrasts Dropdown Menu ----
                        fluidRow(
                            column(3,
                                   conditionalPanel(condition = "(output.downloadbuttonsGenes || output.volcano_GW)",
                                                    id = "lifeStageInputPanel",
                                                    uiOutput("pairwiseSelector_GW")
                                   )
                            ),
                            column(3,
                                   conditionalPanel(condition = "output.downloadbuttonsGenes && input.goLifeStage_GW !=0",
                                                    id = "contrastSelectionPanel_GW",
                                                    uiOutput('contrastDisplaySelection_GW')
                                   )
                                   
                            ),
                            
                            column(9,
                                   conditionalPanel(condition = "(output.downloadbuttonsGenes && input.goLifeStage_GW != 0)",
                                                    uiOutput("volcano_GW")
                                   )
                            )),
                        ## Fluid Row 5: Download Options for Saving DGE Table Results Panel + Differential Gene Expression Table ----
                        fluidRow(
                            column(3,
                                   conditionalPanel(condition = "output.downloadbuttonsGenes && input.goLifeStage_GW != 0 && output.contrastDisplaySelection_GW && output.volcano_GW",
                                                    panel(
                                                        heading = tagList(h5(shiny::icon("fas fa-file-download"),
                                                                             "DGE Table Download Options")),
                                                        status = "primary",
                                                        uiOutput('downloadSelectionMenu_GW'),
                                                        prettyCheckboxGroup("download_DGEdt_direction_GW",
                                                                           h6("Differential Expression Type"),
                                                                           status = "default",
                                                                           icon = icon("check"),
                                                                           choiceNames = c("Upregulated",
                                                                                           "Downregulated",
                                                                                           "No difference"),
                                                                           choiceValues = c("Up",
                                                                                            "Down",
                                                                                            "NotSig"),
                                                                           selected = c("Up",
                                                                                        "Down",
                                                                                        "NotSig")),
                                                       
                                                        textInput("percentDGE_GW",
                                                                  h6("Select Top % of genes, filtered by LogFC value"),
                                                                  "100"),
                                                        h6("Filter across all comparisons?"),
                                                        checkboxInput("download_DGEdt_across_GW",
                                                                      p("Yes, only download genes with selected expression types in all searched pairwise comparisons.")),
                                                        h6("Include names of submitted genes not found in RNA-seq dataset?"),
                                                        checkboxInput("download_missing_genes_GW",
                                                                      p("Yes, include geneIDs for submitted genes not included in RNA-seq dataset.")),
                                                        uiOutput('downloadbuttonGW')
                                                    )
                                                    
                                   )),
                            column(9,
                                   conditionalPanel(condition = "output.downloadbuttonsGenes && input.goLifeStage_GW != 0 && output.contrastDisplaySelection_GW && output.volcano_GW",
                                                    panel(
                                                        heading = tagList(h5(shiny::icon("fas fa-table"),
                                                                             "Pairwise Differential Gene Expression: Table")),
                                                        status = "primary",
                                                        
                                                        DTOutput('highlight.df')        
                                                        
                                                    )
                                   )
                            )
                        )
                    )
           ),
           
           # Life Stage Browser Tab ----
           tabPanel(h4("Browse by Life Stage"),
                    value = "LS",
                    useShinyjs(),
                    div(id = "LS",
                        fluidRow(
                            column(9,
                                   alert(status = "success",
                                         dismissible = TRUE,
                                         id = "userNotes_LS",
                                         "Due to server constraints sessions run on shinyapps.io will time out after several minutes of inactivity. If this happens, please reload the server connection."
                                          )),
                            column(3, 
                                   alert(status = "info",
                                         dismissible = TRUE,
                                         id = "bugReport_LS", 
                                         "To report a bug, please email Dr. Astra Bryant at astrab@uw.edu"))
                        ),
                        
                        ## Fluid Row 2: Species Selection Panel + Life Stage Legend + Download Study Info Dropdown Menu ----
                        fluidRow(
                            column(3,
                                   panel(
                                       heading = tagList(h5(shiny::icon("fas fa-archive"),
                                                            "Select Species")),
                                       status = "primary",
                                       id = "speciesPanelID_LS",
                                       selectInput("selectSpecies_LS",
                                                   h6("Pick a species"),
                                                   choices = list("C. elegans",
                                                                  "C. briggsae",
                                                                  "C. brenneri",
                                                                  "C. japonica",
                                                                  "C. remanei")),
                                       actionButton('speciesLS',
                                                    'Initialize',
                                                    icon = icon("fas fa-share"),
                                                    class = "btn-primary")
                                   )
                            ),
                            column(7,
                                   
                                   conditionalPanel(condition = "input.speciesLS != 0",
                                                    id = "lifeStageLegend_LS",
                                                    uiOutput("Legend_LS"))
                            ),
                            column(2,
                                   conditionalPanel(condition = "input.speciesLS != 0",
                                                    panel(heading = tagList(h5(shiny::icon("fas fa-info-circle"),
                                                                               "Download Study Info Files")),
                                                          status = "default",
                                                          pickerInput("which.Experimental.Info.LS",
                                                                      NULL, 
                                                                      c("Study Design",
                                                                        "Log2CPM Gene Counts",
                                                                        "vDGEList"),
                                                                      options = list(style = 'btn btn-default')),
                                                          uiOutput("StudyInfo.panel.LS")
                                                    )
                                   )
                            )
                        ),
                        
                        
                        ## Fluid Row 3: Input Pairwise Comparisons Panel + Volcano Plot Panel + Contrasts Drop down Menu ----    
                        fluidRow(
                            column(3,
                                   
                                   conditionalPanel(condition = 'input.speciesLS != 0',
                                                    uiOutput("pairwiseSelector_LS")
                                   )),
                            column(3,
                                   conditionalPanel(condition = "output.pairwiseSelector_LS && input.goLS !=0",
                                                    id = "contrastDisplaySelectionPanel_LS",
                                                    uiOutput('contrastDisplaySelection_LS'))
                            ),
                            
                            column(9,
                                   conditionalPanel(condition = "output.pairwiseSelector_LS && input.goLS != 0",
                                                    uiOutput('volcano_LS')               
                                   )
                            )),
                        
                        ## Fluid Row 4: Download Options for Saving DGE Table Results Panel + Differential Gene Expression Table ----
                        
                        fluidRow(
                            column(3,
                                   conditionalPanel(condition = "output.pairwiseSelector_LS && input.goLS != 0 && output.volcano_LS",
                                                    panel(
                                                        heading = tagList(h5(shiny::icon("fas fa-file-download"),
                                                                             "DGE Table Download Options")),
                                                        status = "primary",
                                                        uiOutput('downloadSelectionMenu_LS'),
                                                        prettyCheckboxGroup("download_DGEdt_direction_LS",
                                                                           h6("Differential Expression Type"),
                                                                           status = "default",
                                                                           icon = icon("check"),
                                                                           choiceNames = c("Upregulated",
                                                                                       "Downregulated",
                                                                                       "No difference"),
                                                                           choiceValues = c("Up",
                                                                                            "Down",
                                                                                            "NotSig"),
                                                                           selected = c("Up",
                                                                                        "Down",
                                                                                        "NotSig")),
                                                        textInput("percentDGE_LS",
                                                                  h6("Select top % of genes, filtered by logFC value"),
                                                                  "100"),
                                                        h6("Filter across all comparisons?"),
                                                        checkboxInput("download_DGEdt_across_LS",
                                                                      p("Yes, only download genes with selected expression types in all searched pairwise comparisons.")),
                                                        
                                                    uiOutput('downloadbuttonLS')
                                                    )
                                       
                                   )),
                            column(9,
                                   conditionalPanel(condition = "output.pairwiseSelector_LS && input.goLS != 0 && output.volcano_LS",
                                                    panel(
                                                        heading = tagList(h5(shiny::icon("fas fa-table"),
                                                                             "Pairwise Differential Gene Expression: Table")),
                                                        status = "primary",
                                                        DTOutput('tbl_LS')
                                                    )
                                   )
                            )
                            
                        )
                    )
           ),
           
           
           # About Tab ----
           tabPanel(h4("About (v1.0.0)"),
                    value = "about",
                    fluidRow(
                        column(8,
                               panel(heading =  tagList(h5(shiny::icon("fas fa-question-circle"),
                                                           "App Overview and Features")),
                                     status = "primary",
                                     id = "About_Overview",
                                     includeMarkdown('UI/README/README_Features.md')
                               )
                        ),
                        
                        
                        column(4,
                               panel(heading =  tagList(h5(shiny::icon("fas fa-drafting-compass"),
                                                           "Authors and Release Notes")),
                                     status = "primary",
                                     id = "About_Updates",
                                     includeMarkdown('UI/README/README_Updates.md')
                               )
                        )
                    ),
                        
                    fluidRow(
                        column(8,
                               panel(heading =  tagList(h5(shiny::icon("fas fa-archive"),
                                                           "Data, Pre-Processing, and Analysis Methods")),
                                     status = "primary",
                                     id = "About_Methods",
                                     tabsetPanel(
                                       type = "pills",
                                       tabPanel(
                                         title = tags$em("Data"),
                                         includeMarkdown('UI/README/README_Data.md')
                                       ),
                                       tabPanel(
                                         title = tags$em("Pre-processing"),
                                         includeMarkdown('UI/README/README_Preprocess.md')
                                       ),
                                       tabPanel(
                                         title = tags$em("Analysis Methods"),
                                         includeMarkdown('UI/README/README_Analysis_Methods.md')
                                       ))
                               )
                               ),
                        column(4,
                               panel(heading =  tagList(h5(shiny::icon("fas fa-cloud-download-alt"),
                                                           "Data Availability")),
                                     status = "primary",
                                     p('For each species, the following datasets used can be
        downloaded using the dropdown menu and download button below:',
                                       tags$ol(
                                           tags$li('Study design file (.csv)'),
                                           tags$li('Filtered, normalized log2CPM gene counts (.csv)'),
                                           tags$li('Variance-stabilized Digital Gene Expression List (vDGEList; R object)')
                                       )),
                                     
                                     pickerInput("which.Experimental.Info.About",
                                                 NULL, 
                                                 choices = list(
                                                     `C. elegans` = c("C. elegans Study Design",
                                                                    "C. elegans Log2CPM Gene Counts",
                                                                    "C. elegans vDGEList"),
                                                     `C. briggsae` = c("C. briggsae Study Design",
                                                                         "C. briggsae Log2CPM Gene Counts",
                                                                         "C. briggsae vDGEList"),
                                                     `C. brenneri` = c("C. brenneri Study Design",
                                                                            "C. brenneri Log2CPM Gene Counts",
                                                                            "C. brenneri vDGEList"),
                                                     `C. japonica` = c("C. japonica Study Design",
                                                                      "C. japonica Log2CPM Gene Counts",
                                                                      "C. japonica vDGEList"),
                                                     `C. remanei` = c("C. remanei Study Design",
                                                                      "C. remanei Log2CPM Gene Counts",
                                                                      "C. remanei vDGEList")
                                                 ),
                                                 options = list(style = 'btn btn-primary',
                                                                title = "Select a file to download")),
                                     uiOutput("StudyInfo.panel.About")
                                     
                               )
                        )
                        )
           ),
           
           # Tutorials Tab ----
           tabPanel(h4("Tutorials"),
                    value = "tutorials",
                    
                    fluidRow(
                        column(12,
                               panel(
                                     status = "primary",
                                     id = "Tutorial_Workflow",
                                     img(src='Caen_Tutorial_workflow.png', style = "width: 70vw; display: block; margin-left: auto; margin-right: auto")
                               )
                        )
                    )
                    
           )
           
)
