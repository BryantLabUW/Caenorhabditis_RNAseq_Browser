# Identify the identity of the primary species, and the homologous species
species <- switch(input$selectSpecies_GW,
                  `C. elegans` = 'elegans',
                  `C. briggsae` = 'briggsae',
                  `C. brenneri` = 'brenneri',
                  `C. japonica` = "japonica",
                  `C. remanei` = "remanei")
species.GS1 <- switch(species,
                      'elegans' = 'briggsae',
                      'briggsae' = 'elegans',
                      'brenneri' = 'elegans',
                      'japonica' = 'elegans',
                      'remanei' = 'elegans')
species.GS2 <- switch(species,
                      'elegans' = 'brenneri',
                      'briggsae' = 'brenneri',
                      'brenneri' = 'briggsae',
                      'japonica' = 'briggsae',
                      'remanei' = 'briggsae')
species.GS3 <- switch(species,
                      'elegans' = 'remanei',
                      'briggsae' = 'remanei',
                      'brenneri' = 'remanei',
                      'japonica' = 'remanei',
                      'remanei' = 'japonica')
species.GS4 <- switch(species,
                      'elegans' = 'japonica',
                      'briggsae' = 'japonica',
                      'brenneri' = 'japonica',
                      'japonica' = 'brenneri',
                      'remanei' = 'brenneri')