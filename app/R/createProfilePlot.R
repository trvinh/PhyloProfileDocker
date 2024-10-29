#' Profile plot
#' @param data data for heatmap plot (from reactive fn "dataHeat")
#' @param clusteredDataHeat clustered data (from reactive fn "clusteredDataHeat"
#' @param applyCluster choose clustered data or not (from input$applyCluster)
#' @param parameters plot parameters (colors, size, variable names, ...)
#' @param inSeq subset sequences for customized profile (input$inSeq)
#' @param inTaxa subset taxa for customized profile (input$inTaxa)
#' @param rankSelect selected taxonomy rank (input$rankSelect)
#' @param inSelect selected taxon name (input$inSelect)
#' @param taxonHighlight highlighted taxon (input$taxonHighlight)
#' @param geneHighlight highlighted gene (input$geneHighlight)
#' @param typeProfile either "mainProfile" or "customizedProfile"
#' @param taxDB path to taxonomy database (taxonomyMatrix.txt file required!)
#' @param superRank taxonomy rank for division lines (e.g. superkingdom)
#' @param allTaxa list of all input (super)taxa and their ncbi IDs
#' @param mode plot mode, either "fast" or "normal"
#' @taxDB Path to the taxonomy DB files
#' @return info for selected point on the profile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createProfilePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("plot.ui")),
        br(),
        downloadButton(ns("profileDownload"),"Download profile",
                       class = "butDL"),
        br(),br(),hr(),
        uiOutput(ns("brushedInfo.ui")),
        radioButtons(
            ns("showSettings"), "Show subplot settings",
            choices = c("show", "hide"), selected = "hide", inline = TRUE
        ),
        conditionalPanel(
            condition = "input.showSettings == 'show'",
            column(
                2,
                radioButtons(
                    inputId = ns("xAxisSub"),
                    label = "X-axis",
                    choices = list("taxa", "genes"),
                    selected = "taxa",
                    inline = TRUE
                )
            ),
            column(
                2,
                createPlotSize(ns("widthSub"), "Width (px)", 1200)
            ),
            column(
                2,
                createPlotSize(ns("heightSub"), "Height (px)", 1200)
            ),
            column(
                2,
                createTextSize(ns("fontSizeSub"), "Axis size", 14,100)
            ),
            column(
                2,
                createTextSize(ns("legendSizeSub"), "Legend size", 8,100)
            ),
            column(
                2,
                selectInput(
                    ns("legendPosSub"), label = "Legend position:",
                    choices = list(
                        "Right" = "right", "Left" = "left", "Top" = "top", 
                        "Bottom" = "bottom", "Hide" = "none"
                    ),
                    selected = "right"
                )
            ),
            br(),
            downloadButton(
                ns("profileSubDownload"),"Download subplot", class = "butDL"
            ),
            strong(""),
            downloadButton(
                ns("dataSubDownload"),"Download selected data", class = "butDL"
            ),
            ns = NS(id)
        ),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        )
    )
}

createProfilePlot <- function(
    input, output, session, data, clusteredDataHeat, applyCluster, parameters,
    inSeq, inTaxa, rankSelect, inSelect, taxonHighlight, geneHighlight,
    typeProfile, taxDB, superRank, allTaxa, mode = "normal"
) {
    # data for heatmap ---------------------------------------------------------
    dataHeat <- reactive({
        if (is.null(data())) stop("Profile data is NULL!")
        if (typeProfile() == "customizedProfile") {
            if (is.null(inTaxa()) | is.null(inSeq())) return()
            dataHeat <- dataCustomizedPlot(data(), inTaxa(), inSeq())
            if (applyCluster() == TRUE) {
                dataHeat <- dataCustomizedPlot(
                    clusteredDataHeat(), inTaxa(), inSeq()
                )
            }
        } else {
            dataHeat <- dataMainPlot(data())
            if (applyCluster() == TRUE) {
                dataHeat <- dataMainPlot(clusteredDataHeat())
            }
        }
        # add all input taxa (if available)
        if (!(is.null(allTaxa()))) {
            mergedDf <- merge(
                dataHeat, allTaxa(), by = c("supertaxonID","supertaxon"),
                all = TRUE
            )
            mergedDf$geneID[is.na(mergedDf$geneID)] <- dataHeat$geneID[1]
            mergedDf$supertaxon <- factor(
                mergedDf$supertaxon, levels = levels(allTaxa()$supertaxon)
            )
            return(mergedDf)
        }
        return(dataHeat)
    })

    # basic profile plot -------------------------------------------------------
    basicProfile <- reactive({
        if (is.null(dataHeat())) stop("Profile data is NULL!")
        if (typeProfile() == "customizedProfile") {
            if (length(inSeq()) == 0 || length(inTaxa()) == 0) return()
            if ("all" %in% inSeq() & "all" %in% inTaxa()) return()
        }
        if (mode() == "fast")
            return(heatmapPlottingFast(dataHeat(), parameters()))
        return(heatmapPlotting(dataHeat(), parameters()))
    })

    # get superRank ------------------------------------------------------------
    getSuperRank <- reactive({
        if (!is.null(superRank()) && superRank() == "") return(NULL)
        else return(superRank())
    })

    # render heatmap profile ---------------------------------------------------
    output$plot <- renderPlot({
        if (is.null(basicProfile())) return()
        withProgress(message = 'PLOTTING...', value = 0.5, {
            p <- highlightProfilePlot(
                basicProfile(), dataHeat(), taxonHighlight(), rankSelect(),
                geneHighlight(), taxDB(), parameters()$xAxis
            )
            refLine <- TRUE
            if (mode() == "fast") refLine <- FALSE
            addRankDivisionPlot(
                p, dataHeat(), taxDB(), rankSelect(), getSuperRank(),
                parameters()$xAxis, parameters()$font, 
                parameters()$groupLabelSize, parameters()$groupLabelDist, 
                parameters()$groupLabelAngle, refLine
            )
        })
    })

    output$plot.ui <- renderUI({
        ns <- session$ns

        if (typeProfile() == "customizedProfile") {
            if (is.null(inSeq()[1]) | is.null(inTaxa()[1])) return()
            else if (inSeq()[1] == "all" & inTaxa()[1] == "all")  return()
        }

        shinycssloaders::withSpinner(
            if (mode() == "fast") {
                plotOutput(
                    ns("plot"),
                    width = parameters()$width,
                    height = parameters()$height,
                    brush = ns("plotBrush")
                )
            } else 
                plotOutput(
                    ns("plot"),
                    width = parameters()$width,
                    height = parameters()$height,
                    click = ns("plotClick")
                )
        )
    })

    output$profileDownload <- downloadHandler(
        filename = function() {
            c("profile.svg")
        },
        content = function(file) {
            ggsave(
                file,
                plot = addRankDivisionPlot(
                    basicProfile(), dataHeat(), taxDB(), rankSelect(),
                    getSuperRank(), parameters()$xAxis, parameters()$font,
                    parameters()$groupLabelSize, parameters()$groupLabelDist,
                    parameters()$groupLabelAngle, FALSE
                ),
                width = parameters()$width * 0.035,
                height = parameters()$height * 0.035,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            )
        }
    )

    # get info of brushed area -------------------------------------------------
    brushedData <- reactive({
        if (is.null(input$plotBrush)) return(data.frame())
        dt <- brushedPoints(
            dataHeat(), input$plotBrush, xvar = "supertaxon", yvar = "geneID"
        )
        dt$geneID <- droplevels(dt$geneID)
        dt$supertaxon <- droplevels(dt$supertaxon)
        return(dt)
    })
    
    # create subplot -----------------------------------------------------------
    # ** show/hide subplot settings --------------------------------------------
    observe({
        if (mode() == "fast") {
            updateRadioButtons(
                session, "showSettings", "Show subplot settings",
                choices = c("show", "hide"), selected = "show", inline = TRUE
            )
        } else {
            updateRadioButtons(
                session, "showSettings", "Show subplot settings",
                choices = c("show", "hide"), selected = "hide", inline = TRUE
            )
        }
    })
    
    # ** auto update subplot size ----------------------------------------------
    observe({
        req(nrow(brushedData()) > 0)
        df <- brushedData()
        nrTaxa <- length(unique(as.character(df$supertaxon)))
        nrGene <- length(unique(as.character(df$geneID)))
        adaptedSize <- adaptPlotSize(
            nrTaxa, nrGene, parameters()$xAxis, parameters()$dotZoom
        )
        if (length(adaptedSize) > 0) {
            h <- adaptedSize[1]
            hv <- adaptedSize[2]
            wv <- adaptedSize[3]
            if (hv > 25000) hv <- 25000
            if (wv > 25000) wv <- 25000
            if (h <= 20) {
                updateSelectInput(
                    session, "legendPosSub", label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "top"
                )
                updateNumericInput(session, "widthSub", value = wv  + 50)
            } else if (h <= 30) {
                updateNumericInput(session, "widthSub", value = wv + 50)
            } else {
                updateNumericInput(session, "widthSub", value = wv)
            }
            updateNumericInput(session, "heightSub", value = hv)
        }
    })
    
    # ** get settings for subplot; replace xAxis, font size and legend pos -----
    subPlotParameters <- reactive({
        plotParameters <- parameters()
        plotParameters$xAxis <- input$xAxisSub
        plotParameters$mainLegend <- input$legendPosSub
        plotParameters$xSize <- input$fontSizeSub
        plotParameters$ySize <- input$fontSizeSub
        plotParameters$legendSize <- input$legendSizeSub
        return(plotParameters)
    })
    
    # ** render subplot --------------------------------------------------------
    output$brushedInfo.ui <- renderUI({
        ns <- session$ns
        req(nrow(brushedData()) > 0)
        tagList(
            h4("SELECTED DATA"),
            br(),
            plotOutput(
                ns("subPlot"),
                width = input$widthSub,
                height = input$heightSub,
                click = ns("plotClick")
            )
        )
    })
    
    output$subPlot <- renderPlot({
        req(nrow(brushedData()) > 0)
        heatmapPlotting(brushedData(), subPlotParameters())
    })
    
    # ** download subplot ------------------------------------------------------
    output$profileSubDownload <- downloadHandler(
        filename = function() {
            c("subProfile.svg")
        },
        content = function(file) {
            ggsave(
                file,
                plot = heatmapPlotting(brushedData(), subPlotParameters()),
                width = input$widthSub * 0.035, height = input$heightSub *0.035,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            )
        }
    )
    
    # ** download data of subplot ----------------------------------------------
    output$dataSubDownload <- downloadHandler(
        filename = function() {
            c("subData.out")
        },
        content = function(file){
            dataOut <- brushedData() %>% select(
                geneID, supertaxonID, supertaxon, orthoID, var1, var2, geneName
            )
            if (!(identical(dataOut$geneID, dataOut$geneName)))
                dataOut <- subset(dataOut, select = -c(geneName))
            write.table(
                dataOut, file, sep = "\t", row.names = FALSE, quote = FALSE
            )
        }
    )
    
    # get info of clicked point on heatmap plot --------------------------------
    selectedpointInfo <- reactive({
        # get selected supertaxon name
        taxaList <- getNameList(taxDB())
        rankName <- rankSelect()
        inSelect <- taxaList$ncbiID[taxaList$fullName == inSelect()]
        
        dataHeat <- dataHeat()
        if (mode() == "fast")
            dataHeat <- brushedData()
            
        if (is.null(dataHeat)) {
            message("WARNING: Data for heatmap is NULL!")
            return()
        }
        
        if (typeProfile() == "customizedProfile") {
            # get sub-dataframe of selected taxa and sequences
            if (is.null(inSeq()[1]) | is.null(inTaxa()[1])) {
                message("WARNING: Subset taxa or genes is NULL!")
                return()
            }
            if (inTaxa()[1] == "all" & inSeq()[1] != "all") {
                # select data from dataHeat for selected sequences only
                dataHeat <- subset(dataHeat, geneID %in% inSeq())
            } else if (inSeq()[1] == "all" & inTaxa()[1] != "all") {
                # select data from dataHeat for selected taxa only
                dataHeat <- subset(dataHeat, supertaxon %in% inTaxa())
            } else {
                # select data from dataHeat for selected sequences and taxa
                dataHeat <- subset(dataHeat, geneID %in% inSeq()
                                   & supertaxon %in% inTaxa())
            }
            
            # drop all other supertaxon that are not in sub-dataframe
            dataHeat$supertaxon <- factor(dataHeat$supertaxon)
            dataHeat$geneID <- factor(dataHeat$geneID)
        }
        
        # get values
        if (is.null(input$plotClick$x)) return()
        else {
            # get cooridiate point
            if (parameters()$xAxis == "genes") {
                corX <- round(input$plotClick$y);
                corY <- round(input$plotClick$x)
            } else {
                corX <- round(input$plotClick$x);
                corY <- round(input$plotClick$y)
            }
            # get geneID
            genes <- levels(dataHeat$geneID)
            geneID <- toString(genes[corY])
            # get supertaxon (spec)
            supertaxa <- levels(dataHeat$supertaxon)
            spec <- toString(supertaxa[corX])
            # get var1 and var2 score
            var1 <- NA
            if (!is.na(dataHeat$var1[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var1 <- max(
                    na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec])
                )
            }
            var2 <- NA
            if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var2 <- {
                    max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]))
                }
            }
            # get percentage of present species and total of taxa
            Percent <- NA
            if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                         & dataHeat$supertaxon == spec][1])) {
                Percent <- {
                    max(
                        na.omit(
                            dataHeat$presSpec[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]
                        )
                    )
                }
            }
            presentTaxa <- NA
            if (!is.na(dataHeat$presentTaxa[dataHeat$geneID == geneID
                                            & dataHeat$supertaxon == spec][1])){
                presentTaxa <- max(
                    na.omit(dataHeat$presentTaxa[dataHeat$geneID == geneID
                                                 & dataHeat$supertaxon == spec])
                )
            }
            totalTaxa <- NA
            if (
                !is.na(dataHeat$totalTaxa[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec][1])
            ) {
                totalTaxa <- max(
                    na.omit(
                        dataHeat$totalTaxa[dataHeat$geneID == geneID
                                           & dataHeat$supertaxon == spec]
                    )
                )
            }
            
            
            # get ortholog ID
            orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]
            totalOrtho <- length(orthoID)
            
            # get working taxonomy level
            strain <- "Y"
            if (nlevels(as.factor(dataHeat$totalTaxa)) > 1) strain <- "N"
            
            if (is.na(Percent)) return()
            else {
                info <- c(
                    geneID,
                    list(orthoID),
                    totalOrtho,
                    spec,
                    round(var1, 2),
                    round(Percent, 2),
                    round(var2, 2),
                    presentTaxa,
                    totalTaxa,
                    strain
                )
                return(info)
            }
        }
    })
    
    return(selectedpointInfo)
}
