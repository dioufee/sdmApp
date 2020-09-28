#@utor: Aboubacar HEMA,aboubacarhema94@gmail.com
#@utor: Aboubacar HEMA,aboubacarhema94@gmail.com
#library needed
library(shiny)
library(sf) # classes and functions for vector data
library(raster)# classes and functions for raster data
library(ggplot2)
library(biomod2)
library(grid)
library(rhandsontable)
library(haven)
library(DT)
library(shinyBS)
library(data.table)
library(readxl)
library(shinyFiles)
library(shinydashboard)
library(SSDM)
library(automap)
library(blockCV)
library(tidyverse)
library(ggpubr)
library(CENFA)
library(dismo)
library(randomForest)
library(kernlab)
# No data message
msg_nodata <- function(tab_import=FALSE) {
  if ( tab_import ) {
    txt <- "Please select a dataset to upload"
  } else {
    txt <- "Please upload data in the Microdata tab!"
  }
  fluidRow(
    column(12, h2("No data available", align="center")),
    column(12, strong(txt)))
}

noInputData <- function(prefix="btn_a_micro_", uri) {
  btn <- myActionButton(paste0(prefix, uri), label=("Load dataset"), "primary")
  fluidRow(
    column(12, h3("No input data available!"), class="wb-header"),
    column(12, p("Go to the Data Upload tab to upload a dataset"), class="wb-header-hint"),
    column(12, p("Go back to the Data Upload tab by clicking the button below and load a dataset."), align="center"),
    column(12, div(btn, align="center")))
}
#functions
genObserver_menus <-
  function(pat="btn_results_", n=1, updateVal) {
    res <- paste0('observeEvent(input$',pat,n,', {
                  curid <- "',pat,n,'"
                  nn <- names(input)
                  nn <- nn[grep("',pat,'",nn)]
                  nn <- setdiff(nn, curid)
                  for (btnid in nn) {
                  updateButton(session, btnid, style="default")
                  }
                  obj$',updateVal,' <- "',pat,n,'"
                  updateButton(session, curid, style="primary")
  });
                  ')
    res
  }

# which data.frames ware available in the global environment?
ex <- ls()
if ( length(ex) > 0 ) {
  available_dfs <- ex[sapply(ex, function(x) {
    is.data.frame(get(paste(x)))
  })]
  if ( length(available_dfs) == 0 ) {
    available_dfs <- NULL
  }
}

TimesRasters<-function(x,y){
  z<-x * y
  names(z)<-names(x)
  return(z)
}

myActionButton <- function(inputId, label, btn.style="", css.class="") {
  if ( btn.style %in% c("primary","info","success","warning","danger","inverse","link")) {
    btn.css.class <- paste("btn", btn.style, sep="-")
  } else {
    btn.css.class <- ""
  }
  tags$button(id=inputId, type="button", class=paste("btn action-button", btn.css.class, css.class, collapse=" "), label)
}


summarise_fold<-function(sb){
  records<-sb$records
  records$fold<-1:nrow(records)
  records <- records[,c(5,1,2,3,4)] %>% 
    dplyr::mutate(Pourcentage= round((test_0 + test_1)*100/(test_0 + test_1+train_1 + train_0),digits = 0))
  return(records)
}
ggR_P<-function(rasterLayer){
  samp <- raster::sampleRegular(rasterLayer, 5e+05, asRaster = TRUE)
  map_df <- raster::as.data.frame(samp, xy = TRUE, centroids = TRUE, 
                                  na.rm = TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  basePlot1 <- ggplot2::ggplot() + ggplot2::geom_raster(data = map_df, 
                                                        ggplot2::aes_string(y = "Northing", x = "Easting", fill = "MAP"))
  basePlot1<-basePlot1 + ggplot2::theme_bw() + ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggtitle(label = names(rasterLayer)) + 
    theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
    ggplot2::scale_fill_gradientn(name = "Probability \n of occurence", colours = rev(terrain.colors(10)))
  return(basePlot1)
}

sdmApp_RasterPlot<-function(rasterLayer){
  samp <- raster::sampleRegular(rasterLayer, 5e+05, asRaster = TRUE)
  map_df <- raster::as.data.frame(samp, xy = TRUE, centroids = TRUE, 
                                  na.rm = TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  basePlot1 <- ggplot2::ggplot() + ggplot2::geom_raster(data = map_df, 
                                                        ggplot2::aes_string(y = "Northing", x = "Easting", fill = "MAP"))
  # basePlot1<-basePlot1 + ggplot2::theme_bw() + ggplot2::labs(x = "Longitude", y = "Latitude") +
  #   ggtitle(label = names(rasterLayer)) + 
  #basePlot1<-basePlot1 + theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
  basePlot1<-basePlot1 + ggtitle(label = names(rasterLayer))  #+ ggplot2::scale_fill_gradientn(name = " ", colours = rev(terrain.colors(10)))
  return(basePlot1)
}



PASpecies<-function(rasterLayer){
  samp <- raster::sampleRegular(rasterLayer, 5e+05, asRaster = TRUE)
  map_df <- raster::as.data.frame(samp, xy = TRUE, centroids = TRUE, 
                                  na.rm = TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  map_df$MAP<-factor(map_df$MAP)
  basePlot <- ggplot2::ggplot() + ggplot2::geom_raster(data = map_df, 
                                                       ggplot2::aes_string(y = "Northing", x = "Easting", fill = "MAP"))
  basePlot<-basePlot + ggplot2::theme_bw() + ggplot2::labs(x = "Longitude", y = "Latitude") + ggtitle(label = names(rasterLayer)) + theme(plot.title = element_text(hjust = 0.5, size = 10))+scale_fill_manual(values=c("white","green"),name="Specie",labels=c("Absence","Presence"))
  
  return(basePlot)
}
Explorer<-function (blocks, rasterLayer, speciesData, num) {
  # records<-blocks$records
  # records$fold<-1:nrow(records)
  # records <- records[,c(5,1,2,3,4)] %>% 
  #   mutate(calcul= round((test_0 + test_1)*100/(test_0 + test_1+train_1 + train_0),digits = 0))
  # records.p <- ggtexttable(records,rows = NULL, theme = ttheme("mGreen"))
  # 
  polyObj <- blocks$blocks
  folds <- blocks$folds
  kmax <- length(folds)
  species <- blocks$species
  speciesData <- sf::st_as_sf(speciesData)
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster = TRUE)
  map_df <- raster::as.data.frame(samp, xy = TRUE, centroids = TRUE, 
                                  na.rm = TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  mid <- stats::median(map_df$MAP)
  basePlot <- ggplot2::ggplot() + ggplot2::geom_raster(data = map_df, 
                                                       ggplot2::aes_string(y = "Northing", x = "Easting", fill = "MAP")) + 
    ggplot2::scale_fill_gradient2(low = "darkred", mid = "yellow", 
                                  high = "darkgreen", midpoint = mid) + ggplot2::guides(fill = FALSE) + 
    ggplot2::theme_bw() + ggplot2::labs(x = "", y = "")
  trainSet <- unlist(folds[[num]][1])
  testSet <- unlist(folds[[num]][2])
  training <- speciesData[trainSet, ]
  testing <- speciesData[testSet, ]
  plotPoly <- polyObj[polyObj$folds ==num,]               
  plotPoly <- sf::st_as_sf(plotPoly)
  if (is.null(species)) {
    if (class(blocks) == "SpatialBlock") {
      ptr <- basePlot + ggplot2::geom_sf(data = plotPoly, 
                                         color = "red", fill = "orangered4", alpha = 0.04, 
                                         size = 0.2) + ggplot2::geom_sf(data = training, 
                                                                        alpha = 0.7, color = "blue", size = 2) + 
        ggplot2::ggtitle("Training set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
      pts <- basePlot + ggplot2::geom_sf(data = plotPoly, 
                                         color = "red", fill = "orangered4", alpha = 0.04, 
                                         size = 0.2) + ggplot2::geom_sf(data = testing, 
                                                                        alpha = 0.7, color = "blue", size = 2) + 
        ggplot2::ggtitle("Testing set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
    else {
      ptr <- basePlot + ggplot2::geom_sf(data = training, 
                                         alpha = 0.7, color = "blue", size = 2) + 
        ggplot2::ggtitle("Training set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
      pts <- basePlot + ggplot2::geom_sf(data = testing, 
                                         alpha = 0.7, color = "blue", size = 2) + 
        ggplot2::ggtitle("Testing set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
  }
  else {
    if (class(blocks) == "SpatialBlock") {
      ptr <- basePlot + ggplot2::geom_sf(data = plotPoly, 
                                         color = "red", fill = "orangered4", alpha = 0.04, 
                                         size = 0.2) + ggplot2::geom_sf(data = training, 
                                                                        ggplot2::aes(color = get(species)), show.legend = "point", 
                                                                        alpha = 0.7, size = 2) + ggplot2::labs(color = species) + 
        ggplot2::ggtitle("Training set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
      pts <- basePlot + ggplot2::geom_sf(data = plotPoly, 
                                         color = "red", fill = "orangered4", alpha = 0.04, 
                                         size = 0.2) + ggplot2::geom_sf(data = testing, 
                                                                        ggplot2::aes(color = get(species)), show.legend = "point", 
                                                                        alpha = 0.7, size = 2) + ggplot2::labs(color = species) + 
        ggplot2::ggtitle("Testing set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
    else {
      ptr <- basePlot + ggplot2::geom_sf(data = training, 
                                         ggplot2::aes(color = get(species)), show.legend = "point", 
                                         alpha = 0.7, size = 2) + ggplot2::labs(color = species) + 
        ggplot2::ggtitle("Training set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
      pts <- basePlot + ggplot2::geom_sf(data = testing, 
                                         ggplot2::aes(color = get(species)), show.legend = "point", 
                                         alpha = 0.7, size = 2) + ggplot2::labs(color = species) + 
        ggplot2::ggtitle("Testing set") + theme(plot.title = element_text(hjust = 0.5, size = 10))
    }
  }
  plot(ggpubr::ggarrange(ptr, pts,common.legend = TRUE))
}

cirad_hema <-function()
{
  #ui
  ui<-navbarPage(id="cirad","SDMs GUI",theme = "readtable",
                 tabPanel("Help/About",
                          uiOutput("ui_about")),
                 tabPanel("Data Upload",uiOutput("ui_import_data")),
                 tabPanel("Data Preparation",uiOutput("ui_preparation")),
                 # tabPanel("ggplot2",uiOutput("ui_show_microdata")),        
                 # tabPanel("Interactive map"),
                 tabPanel("Spatial blocking",uiOutput("ui_blockCV")),
                 tabPanel("Modeling",uiOutput("ui_Models")),
                 tabPanel("R-Code")
                 
                 
  )
  
  
  
  
  
  
  
  #server
  server<-function(input, output, session){
    data <- reactiveValues(Env = stack(), Occ = data.frame(), dir = getwd(), ESDM = NULL, esdms = list(), Stack = NULL)
    obj<-reactiveValues()
    obj$inputdata <- NULL
    obj$sdcObj <- NULL
    #obj$sdcObj <- createSdcObj(testdata,
    #  keyVars=c('roof','walls','water'),
    #  numVars=c('expend','income','savings'), strataVar="sex", w='sampling_weight')
    obj$code_read_and_modify <- c()
    obj$code_setup <- c()
    obj$code_anonymize <- c()
    obj$code <- "require(sdcMicro)"
    obj$transmat <- NULL
    obj$last_warning <- NULL
    obj$last_error <- NULL
    obj$comptime <- 0
    obj$microfilename <- NULL # name of uploaded file
    obj$lastaction <- NULL
    obj$anon_performed <- NULL # what has been applied?
    obj$rbs <- obj$sls <- NULL
    obj$setupval_inc <- 0
    obj$inp_sel_viewvar1 <- NULL
    obj$inp_sel_anonvar1 <- NULL
    obj$lastreport <- NULL # required to show the last saved report
    obj$lastdataexport <- NULL # required to show the last saved exported data
    obj$lastproblemexport <- NULL # required to show the last exported sdcproblem
    obj$lastproblemexport1 <- NULL # required to show the last exported sdcproblem (undo-page)
    obj$lastscriptexport <- NULL # required to show the last saved script
    obj$ldiv_result <- NULL # required for l-diversity risk-measure
    obj$suda2_result <- NULL # required for suda2 risk-measure
    obj$hhdata <- NULL # household-file data required for merging
    obj$hhdata_applied <- FALSE # TRUE, if mergeHouseholdData() has been applied
    obj$hhdata_selected <- FALSE # TRUE, if selectHouseholdData() has been applied
    
    # stores the current selection of the relevant navigation menus
    obj$cur_selection_results <- "btn_results_1" # navigation for Results/Risks
    obj$cur_selection_exports <- "btn_export_results_1" # navigation for export
    obj$cur_selection_script <- "btn_export_script_1" # navigation for reproducibility/script
    obj$cur_selection_microdata <- "btn_menu_microdata_1" # navigation for microdata
    obj$cur_selection_import <- "btn_import_data_1" # navigation for import
    obj$cur_selection_anon <- "btn_sel_anon_1" # navigation for anonymization
    #obj$cur_selection_imputation<-"btn_impute_results_1" # navigation for imputation
    # for stata-labelling
    obj$stata_labs <- NULL
    obj$stata_varnames <- NULL
    
    ####################data
    
    ###############################################################
    
    ######################## About/Help content ##########################################
    #
    
    output$ui_about <- renderUI({
      out <- fluidRow(
        column(width = 8, offset = 2, h2(("sdmApp"))),
        column(width = 8, offset = 2, p("This graphical user interface of",code("sdmgui")," offers the possibility to run 10 state-of-the-art modeling techniques to describe and model the relationships between a given species and its environment. It is an attempt to define the ecological niche of a particular species using environmental variables (temperature, precipitation, ...) with the potential use of making, for instance, future projections under climate and land use change scenarios. Although it has been mostly developed for ecologists that aim to predict species distribution, biomod2 can also be used to model any binomial data (for instance, gene, markers, ecosystem...) in function of any explanatory variables. 
                                        Even if you are not an expert in the",code("R"),"programming language. Detailed information on how to use this graphical user-interface (GUI) can be found in a tutorial (a so-called vignette) that is included in the",code("ensae"),"package.
                                        The vignette is available on",tags$a("GitHub pages", href="https://github.com/Abson-dev", target="_blank") 
                                        
                                        
        )))
      
      out <- list(out, fluidRow(
        column(width = 8, offset = 2, h4(("Contact and Feedback"))),
        column(width = 8, offset = 2, p("In case you have any suggestions or bug reports, please file an issue at the",
                                        tags$a("issue tracker", href="https://github.com/Abson-dev", target="_blank"),"in our",
                                        tags$a("GitHub repo", href="https://github.com/Abson-dev", target="_blank"),".")),
        column(width = 8, offset = 2, p("Before reporting any bugs, please make sure that you are working with an up-to-date",tags$b("R"),"installation and
                                        that all packages have been updated. You can do so by entering",code("update.packages(ask=FALSE)"),"into your",code("R"),"prompt."))
        ))
      out
    })
    ######################## End About/Help content ##################################################################################
    ##################################################################################################################################
    ########################  Data Upload content ##################################################################################
    ##################################################################################################################################
    #Import data
    # specific (gui)-options for csv-import
    
    # Environmental variable loading
    load.var <- reactiveValues(factors = c(), formats = c(), norm = TRUE,  vars = list())
    working.directory <- "C:\\Users\\DELLDRAMOMO\\Desktop\\Package\\data\\"
    example = system.file("extdata", package = "SSDM")
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'envfiles', session=session,
                      roots=c(wd = working.directory,
                              example = example,
                              home = '/home',
                              root = '/'),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'envfiles', session=session,
                      roots=c(wd = working.directory,
                              example = example,
                              disks),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    } else {
      shinyFileChoose(input, 'envfiles', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/user',
                                root = '/'),
                      filetypes=c('',"grd", "tif", "asc","sdat", "rst", "nc", "tif", "envi", "bil", "img"))
    }
    observeEvent(input$envfiles,{
      if(!is.integer(input$envfiles)){
        load.var$vars = lapply(input$envfiles$files, function(x) x[[length(x)]])
        names(load.var$vars) <- unlist(load.var$vars)
      }
    })
    
    output$factors <- renderUI({
      selectInput('factors', 'Categorical', load.var$vars, multiple = TRUE, selectize = TRUE)
    })
    observeEvent(input$load, {
      validate(
        need(length(load.var$vars) > 0, 'Choose environment variable files first !')
      )
      if(Sys.info()[['sysname']] == 'Linux') {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      'home' = '/home',
                      'root' = '/')
      } else if (Sys.info()[['sysname']] == 'Windows') {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      input$envfiles$root)
      } else {
        path = switch(input$envfiles$root,
                      'wd' = working.directory,
                      'example' = example,
                      'home' = '/home',
                      'root' = '/')
      }
      for(i in 2:(length(input$envfiles$files[[1]]))-1){
        path = paste0(path, '/', input$envfiles$files[[1]][i])
      }
      load.var$formats = c()
      for (i in seq_len(length(load.var$vars))) {
        format = paste0('.',strsplit(load.var$vars[[i]], '.', fixed = TRUE)[[1]][2])
        if (!(format %in% load.var$formats)) {load.var$formats = c(load.var$formats, format)}
      }
     
      a = try(withProgress(message = 'Variables loading',
                           load_var(path,
                                    files = unlist(load.var$vars),
                                    format = load.var$formats,
                                    Norm = FALSE,
                                    tmp = FALSE,
                                    categorical = load.var$factors,
                                    verbose = FALSE,
                                    GUI = TRUE)))
      if(inherits(a, 'try-error')){
        output$Envbug <- renderUI(p('Environmental variables loading failed, please check your inputs and try again'))
      } else {
        output$Envbug <- renderUI(p())
        data$Env = a
        for (i in seq_len(length(load.var$vars))) {
          names(data$Env)[i] = strsplit(load.var$vars[[i]], '.', fixed = TRUE)[[1]][1]
        }
        output$layerchoice <- renderUI({
          selectInput('layer', 'Variable', as.list(names(data$Env)), multiple = FALSE, selectize = TRUE)
          
        })
        width <- reactive({
          input$fig_width
        })
        height <- reactive({
          input$fig_height
        })
        string_code <- reactive({
          p <- paste("sdmApp_RasterPlot(map)")
          p <- paste(p, "+ scale_fill_","gradientn", "(name = 'Value',  colours = rev(terrain.colors(10)))",
                     sep = "")
          #p <- paste("+ theme(plot.title = element_text(hjust = 0.5, size = 10))")
          if (input$label_axes) 
            p <- paste(p, "+ labs(x = 'input$lab_x', y = 'input$lab_y')")
          if (input$add_title) 
            p <- paste(p, "+ ggtitle('input$title')")
          if (input$adj_leg == "Change legend")
            p <- paste(p, "+ scale_fill_","gradientn", "(name = 'input$leg_ttl',  colours = rev(terrain.colors(10)))",
                       sep = "")
          # if (input$adj_col)
          #   p <- paste(p, "+ scale_fill_","gradientn", "(name = 'input$leg_ttl',  colours = rev(terrain.colors(10)))",
          #     sep = "")
          p <- paste(p, "+", input$theme)
          if (input$adj_fnt_sz || input$adj_fnt || input$rot_txt || 
              input$adj_leg != "Keep legend as it is" || 
              input$adj_grd) {
            p <- paste(p, paste(" + theme(\n    ", 
                                "plot.title = element_text(hjust = 0.5, size = 10),\n    ",
                                if (input$adj_fnt_sz) 
                                  "axis.title = element_text(size = input$fnt_sz_ttl),\n    ", 
                                if (input$adj_fnt_sz) 
                                  "axis.text = element_text(size = input$fnt_sz_ax),\n    ", 
                                if (input$adj_fnt) 
                                  "text = element_text(family = 'input$font'),\n    ", 
                                if (input$rot_txt) 
                                  "axis.text.x = element_text(angle = 45, hjust = 1),\n    ", 
                                if (input$adj_leg == "Remove legend") 
                                  "legend.position = 'none',\n    ", 
                                if (input$adj_leg == "Change legend") 
                                  "legend.position = 'input$pos_leg',\n    ", 
                                if (input$grd_maj) 
                                  "panel.grid.major = element_blank(),\n    ", 
                                if (input$grd_min) 
                                  "panel.grid.minor = element_blank(),\n    ", 
                                ")", sep = ""), sep = "")
          }
          p <- str_replace_all(p, c(`input\\$lab_x` = as.character(input$lab_x), 
                                    `input\\$lab_y` = as.character(input$lab_y), 
                                    `input\\$title` = as.character(input$title), 
                                    `input\\$palet` = as.character(input$palet), 
                                    `input\\$fnt_sz_ttl` = as.character(input$fnt_sz_ttl), 
                                    `input\\$fnt_sz_ax` = as.character(input$fnt_sz_ax), 
                                    `input\\$font` = as.character(input$font), 
                                    `input\\$leg_ttl` = as.character(input$leg_ttl), 
                                    `input\\$pos_leg` = as.character(input$pos_leg))
          )
          p <- str_replace_all(p, ",\n    \\)", "\n  \\)")
          p
        })
        output$env <- renderPlot(width = width, height = height,{
          if(!is.null(input$layer)){
            i = as.numeric(which(as.list(names(data$Env)) == input$layer))
            if(data$Env[[i]]@data@isfactor) {
              map = !as.factor(data$Env[[i]])
            } else {
              map = data$Env[[i]]
            }
            a =try(eval(parse(text = string_code())))
            if(inherits(a, 'try-error')){
              output$Envbugplot <- renderUI(p('Can not plot this raster! Please verify it and try again.'))
            }
            else{
              output$Envbugplot <- renderUI(p())
              a
            }
          }
        })
        # observeEvent(input$export_raster_plot,{
        #   ggsave(paste0(working.directory,input$layer,".png"),a)
        # })
      }
      updateTabItems(session, "actions", selected = "newdata")
    })
    
    # Occurrences loading
    #load.occ <- reactiveValues(columns = c())
    load.occ <- reactiveValues()
   
# type_file <-reactive({
#         if(input$file_type=="text"){
#               type_file=c('',"csv", "txt")}
#             else {
#               if(input$file_type=="Excel"){
#                 type_file=c('',"xlsx", "xls")
#                 }
#               else{
#                 if(input$file_type=="SPSS"){
#                   type_file=c('',"sav", "zsav","por")}
#                 else{
#                   if(input$file_type=="Stata"){
#                     type_file=c('',"dta")}
#                   else{if(input$file_type == "SAS"){type_file=c('',"sas7bdat")}}
#                   }
#                  }
#               
#             }
#               
#               
#               type_file
#             })
######################################################################"
observeEvent(input$file_type,{
  if(input$file_type=="text"){
    load.occ$type_file=c('',"csv", "txt")}
  else {
    if(input$file_type=="Excel"){
      load.occ$type_file=c('',"xlsx", "xls")
    }
    else{
      if(input$file_type=="SPSS"){
        load.occ$type_file=c('',"sav", "zsav","por")}
      else{
        if(input$file_type=="Stata"){
          load.occ$type_file=c('',"dta")}
        else{if(input$file_type == "SAS"){load.occ$type_file=c('',"sas7bdat")}}
      }
    }
    
  }
    if(Sys.info()[['sysname']] == 'Linux') {
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/home',
                                root = '/'),
                      filetypes=load.occ$type_file)
    } else if (Sys.info()[['sysname']] == 'Windows') {
      d = system('wmic logicaldisk get caption', intern = TRUE)
      disks = c()
      for(i in 2:(length(d)-1)){
        disks = c(disks, substr(d[i],1,2))
      }
      names(disks) = disks
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                disks),
                      filetypes=load.occ$type_file)
    } else {
      shinyFileChoose(input, 'Occ', session=session,
                      roots = c(wd = working.directory,
                                example = example,
                                home = '/user',
                                root = '/'),
                      filetypes=load.occ$type_file)
    }
})
  ###################################
    observeEvent(input$Occ, {
      if(!is.integer(input$Occ)) {
        file = paste0(switch(input$Occ$root,
                             'wd' = working.directory,
                             'example' = example,
                             'home' = '/home',
                             'root' = '/',
                             input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
        if(input$file_type=="text"){
          load.occ$columns = names(read.csv2(file))
          load.occ$df_occ<-read.csv2(file)
          observeEvent(input$sep, {
            if(!is.integer(input$Occ)) {
              file = paste0(switch(input$Occ$root,
                                   'wd' = working.directory,
                                   'example' = example,
                                   'home' = '/home',
                                   'root' = '/',
                                   input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
              load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
              load.occ$df_occ<-read.csv2(file, sep = input$sep, nrows = 0)
            }
          })
          observeEvent(input$Occ, {
            if(!is.integer(input$Occ)) {
              file = paste0(switch(input$Occ$root,
                                   'wd' = working.directory,
                                   'example' = example,
                                   'home' = '/home',
                                   'root' = '/',
                                   input$Occ$root), '/', paste0(unlist(input$Occ$files[[1]])[-1], collapse = '/'))
              load.occ$columns = names(read.csv2(file, sep = input$sep, nrows = 0))
              load.occ$df_occ<-read.csv2(file, sep = input$sep, nrows = 0)
              
            }
          })
        }
        else if (input$file_type == "Excel") {
          load.occ$columns <- names(read_excel(file))
          load.occ$df_occ<-read_excel(file)
        }
        else if (input$file_type == "SPSS") {
          load.occ$columns <- names(read_sav(file))
          load.occ$df_occ<-read_sav(file)
        }
        else if (input$file_type == "Stata") {
          load.occ$columns <- names(read_dta(file))
          load.occ$df_occ<-read_dta(file)
        }
        else if (input$file_type == "SAS") {
          load.occ$columns <- names(read_sas(file))
          load.occ$df_occ<-read_sas(file)
        }
      }
    })
    
    ##############################

    ####################################################""
     output$Xcol <- renderUI({selectInput('Xcol', 'Longitude (X)', load.occ$columns, multiple = FALSE)})
     observeEvent(input$Xcol,{
       load.occ$Ycolumns<-setdiff(load.occ$columns,input$Xcol)
     output$Ycol <- renderUI({selectInput('Ycol', 'Latitude (Y)', load.occ$Ycolumns, multiple = FALSE)})
     observeEvent(input$Ycol,{
     load.occ$Pcol<-setdiff(load.occ$Ycolumns,input$Ycol)
     output$Pcol <- renderUI({selectInput('Pcol', 'Specie column', load.occ$Pcol, multiple = FALSE)})
     })
     })
    observeEvent(input$load2, {
      validate(
        need(length(data$Env@layers) > 0, 'You need to load environmental variable before !'),
        need(length(input$Occ) > 0, 'Choose occurrences file first !')
      )
    load.occ$select<-load.occ$df_occ[,c(input$Xcol,input$Ycol,input$Pcol)]
    load.occ$lon<-input$Xcol
    load.occ$lat<-input$Ycol
    load.occ$spec_select<-input$Pcol
    
     })
    
    ################
     occ_data_df = reactive({
       datatable(load.occ$df_occ,
                 rownames = FALSE,
                 selection="none",
                 options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=20)
       )
     })
     #, options = list(scrollX=TRUE, lengthMenu=list(c(10, 25, 100, -1), c('10', '20', '100', 'All')), pageLength=25), filter="top", rownames=FALSE
     output$occ <- DT::renderDataTable({
       occ_data_df()
     })
    
    output$ui_import_data <- renderUI({
      txt_setup <- "First load rasters file, after species occurence file"
      out<-NULL
      out <- fluidRow(
        column(width = 12, offset = 0, h3("Uploading environmental variables and occurrence table"), class="wb-header"),
        column(width = 12, offset = 0, p("Load the dataset."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      txt_rasters_info<-paste0("You have" ,code(raster::nlayers(data$Env)),"layers.The extent is xmin=",code(raster::extent(data$Env)@xmin),",xmax=",code(raster::extent(data$Env)@xmax),",ymin=",code(raster::extent(data$Env)@ymin),",ymax=",code(raster::extent(data$Env)@ymax))
      out<-list(out,
                sidebarPanel(width = 3,
                  p('Load environmental rasters for model building or model forecasting'),
                  uiOutput('Envbug'),
                  shinyFilesButton('envfiles', 'Raster selection', 'Please select rasters', FALSE, multiple = TRUE),
                  selectInput("categorical_var","Categorical variable",
                              c(No = "No categorical",Yes = "Categorical")),
                  conditionalPanel(
                    condition = "input.categorical_var == 'Categorical'",
                    p('Which variable should be considered as a categorical variable?'),
                    uiOutput('factors')
                  ),
                  myActionButton("load",label=("Load data"), "primary")),
                #actionButton('load', 'Load')),
                mainPanel(width = 6, tabsetPanel(type = "tabs",
                                                 tabPanel("About",
                                                          p(HTML(txt_rasters_info))
                                                 ),
                                                 tabPanel("Preview",
                                                          
                                                          uiOutput('layerchoice'),
                                                          myActionButton("export_raster_plot",label=("Export"), "primary"),
                                                          uiOutput('Envbugplot'),
                                                          plotOutput('env')
                                                 )
                                                 
                                                 
                ),
                id = "tabs"),
                sidebarPanel(width = 3, h4("Change aesthetics"), 
                             tabsetPanel(tabPanel("Text", checkboxInput(inputId = "label_axes", 
                                                                        label = strong("Change labels axes"), 
                                                                        value = FALSE), conditionalPanel(condition = "input.label_axes == true", 
                                                                                                         textInput("lab_x", "X-axis:", value = "label x-axis")), 
                                                  conditionalPanel(condition = "input.label_axes == true", 
                                                                   textInput("lab_y", "Y-axis:", 
                                                                             value = "label y-axis")), checkboxInput(inputId = "add_title", 
                                                                                                                     label = strong("Add title"), value = FALSE), 
                                                  conditionalPanel(condition = "input.add_title == true", 
                                                                   textInput("title", "Title:", 
                                                                             value = "Title")), checkboxInput(inputId = "adj_fnt_sz", 
                                                                                                              label = strong("Change font size"), 
                                                                                                              value = FALSE), conditionalPanel(condition = "input.adj_fnt_sz == true", 
                                                                                                                                               numericInput("fnt_sz_ttl", "Size axis titles:", 
                                                                                                                                                            value = 12), numericInput("fnt_sz_ax", 
                                                                                                                                                                                      "Size axis labels:", value = 10)), 
                                                  checkboxInput(inputId = "rot_txt", label = strong("Rotate text x-axis"), 
                                                                value = FALSE), checkboxInput(inputId = "adj_fnt", 
                                                                                              label = strong("Change font"), value = FALSE), 
                                                  conditionalPanel(condition = "input.adj_fnt == true", 
                                                                   selectInput("font", "Font", choices = c("Courier", 
                                                                                                           "Helvetica", "Times"), selected = "Helvetica"))), 
                                         tabPanel("Theme", conditionalPanel(condition = "input.group != '.'", 
                                                                            checkboxInput(inputId = "adj_col", 
                                                                                          label = strong("Change colours"), 
                                                                                          value = FALSE), conditionalPanel(condition = "input.adj_col", 
                                                                                                                           selectInput(inputId = "palet", label = strong("Select palette"), 
                                                                                                                                       choices = list(Qualitative = c("Accent", 
                                                                                                                                                                      "Dark2", "Paired", "Pastel1", 
                                                                                                                                                                      "Pastel2", "Set1", "Set2", 
                                                                                                                                                                      "Set3"), Diverging = c("BrBG", 
                                                                                                                                                                                             "PiYG", "PRGn", "PuOr", 
                                                                                                                                                                                             "RdBu", "RdGy", "RdYlBu", 
                                                                                                                                                                                             "RdYlGn", "Spectral"), 
                                                                                                                                                      Sequential = c("Blues", "BuGn", 
                                                                                                                                                                     "BuPu", "GnBu", "Greens", 
                                                                                                                                                                     "Greys", "Oranges", "OrRd", 
                                                                                                                                                                     "PuBu", "PuBuGn", "PuRd", 
                                                                                                                                                                     "Purples", "RdPu", "Reds", 
                                                                                                                                                                     "YlGn", "YlGnBu", "YlOrBr", 
                                                                                                                                                                     "YlOrRd")), selected = "set1"))), 
                                                  conditionalPanel(condition = "input.jitter", 
                                                                   checkboxInput("adj_jitter", strong("Change look jitter"), 
                                                                                 FALSE), conditionalPanel(condition = "input.adj_jitter", 
                                                                                                          textInput("col_jitter", "Colour (name or RGB):", 
                                                                                                                    value = "black"), numericInput("size_jitter", 
                                                                                                                                                   "Size:", value = 1), sliderInput("opac_jitter", 
                                                                                                                                                                                    "Opacity:", min = 0, max = 1, 
                                                                                                                                                                                    value = 0.5, step = 0.01), sliderInput("width_jitter", 
                                                                                                                                                                                                                           "Width jitter:", min = 0, max = 0.5, 
                                                                                                                                                                                                                           value = 0.25, step = 0.01))), checkboxInput("adj_grd", 
                                                                                                                                                                                                                                                                       strong("Remove gridlines"), FALSE), 
                                                  conditionalPanel(condition = "input.adj_grd", 
                                                                   checkboxInput("grd_maj", strong("Remove major gridlines"), 
                                                                                 FALSE), checkboxInput("grd_min", 
                                                                                                       strong("Remove minor gridlines"), 
                                                                                                       FALSE)), selectInput("theme", "Theme", 
                                                                                                                            choices = c(bw = "theme_bw()", classic = "theme_classic()", 
                                                                                                                                        dark = "theme_dark()", grey = "theme_grey()", 
                                                                                                                                        light = "theme_light()", line_draw = "theme_linedraw()", 
                                                                                                                                        minimal = "theme_minimal()"), selected = "theme_bw()")), 
                                         tabPanel("Legend", conditionalPanel(condition = "input.group != '.'", 
                                                                             radioButtons(inputId = "adj_leg", label = NULL, 
                                                                                          choices = c("Keep legend as it is", 
                                                                                                      "Remove legend", "Change legend"), 
                                                                                          selected = "Keep legend as it is"), 
                                                                             conditionalPanel(condition = "input.adj_leg=='Change legend'", 
                                                                                              textInput("leg_ttl", "Title legend:", 
                                                                                                        value = "title legend"), selectInput("pos_leg", 
                                                                                                                                             "Position legend", choices = c("right", 
                                                                                                                                                                            "left", "top", "bottom"))))), 
                                         tabPanel("Size", checkboxInput("fig_size", 
                                                                        strong("Adjust plot size on screen"), 
                                                                        FALSE), conditionalPanel(condition = "input.fig_size", 
                                                                                                 numericInput("fig_height", "Plot height (# pixels): ", 
                                                                                                              value = 400), numericInput("fig_width", 
                                                                                                                                         "Plot width (# pixels):", value = 480)), 
                                                  checkboxInput("fig_size_download", 
                                                                strong("Adjust plot size for download"), 
                                                                FALSE), conditionalPanel(condition = "input.fig_size_download", 
                                                                                         numericInput("fig_height_download", 
                                                                                                      "Plot height (in cm):", value = 14), 
                                                                                         numericInput("fig_width_download", 
                                                                                                      "Plot width (in cm):", value = 14)))))
      )
      out
      #####################################################
      out<-list(out,
                column(12,
                       sidebarPanel(width = 3,
                         p('Occurrence table'),
                         selectInput("file_type","Type of file:", list(`text (csv)` = "text", 
                                                                       Excel = "Excel", SPSS = "SPSS", 
                                                                       Stata = "Stata", SAS = "SAS"), selected = "text"),
                         shinyFilesButton('Occ', 'Occurrence selection', 'Please select occurrence file', FALSE),
                         conditionalPanel(condition = "input.file_type=='text'",
                                          radioButtons('sep', 'Separator',
                                                       c(Comma = ',',
                                                         Semicolon = ';',
                                                         Tab = '\t',
                                                         'White space' = ' '),
                                                       ',', inline = TRUE),
                                          radioButtons('dec', 'Decimal',
                                                       c(Point ='.',
                                                         Comma = ','),
                                                       '.', inline = TRUE)),
                         uiOutput('Xcol'),
                         uiOutput('Ycol'),
                         uiOutput('Pcol'),
                         myActionButton("load2",label=("Load"), "primary")
                         #actionButton('load2', 'Select studied specie')
                       ),
                       mainPanel(width = 7, tabsetPanel(type = "tabs",
                                                        tabPanel("Preview",
                                                                 uiOutput('Occbug'),
                                                                 dataTableOutput('occ')))
                                 ,
                                 id = "tabs"))
      )
    }
    )
    ########################################### End Data upload ##############
    ##########################################################################
    
    ###########################################"" Data Preparation#############
    ###############################################################"###########
    
    ######
    output$ui_view_species_data <- renderUI({
      ###############
      ###species selected
      occ_data_select_df = reactive({
        datatable(load.occ$select,
                  rownames = FALSE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=10)
        )
      })
      
      output$occ_data_select <- DT::renderDataTable({
        occ_data_select_df()
      })
      
      ###############  function ######################
      # all variables available in the input data set
      allVars <- reactive({
        inp <- load.occ$select
        if (is.null(inp)) {
          return(NULL)
        }
        cn <- colnames(inp)
        cl <- sapply(1:ncol(inp), function(x) {
          class(inp[[x]])
        })
        names(cn) <- paste0(cn," (",cl,")")
        cn
      })
      ###############  end function ######################
      ###############  function ######################
      dataTypes <- reactive({
        inputdata <- load.occ$select
        if (is.null(inputdata)) {
          return(NULL)
        }
        cn <- colnames(inputdata)
        cl <- sapply(1:ncol(inputdata), function(x) {
          class(inputdata[[x]])
        })
        cl
      })
      
      output$SpeciesTable <- DT::renderDataTable({
        ############# function ########
        sdmData <- reactive({
          inputdata <- load.occ$select
          if (is.null(inputdata)) {
            return(NULL)
          }
          vars <- allVars()
          df <- data.frame(
            "Variable Name"=vars
          )
          df$nrCodes <- sapply(inputdata, function(x) { length(unique(x))} )
          df$nrNA <- sapply(inputdata, function(x) { sum(is.na(x))} )
          df$min<-sapply(inputdata, function(x) { min(x,na.rm = TRUE)} )
          df$max<-sapply(inputdata, function(x) { max(x,na.rm = TRUE)} )
          colnames(df) <- c("Variable name",  "Number of levels", "Number of missing","minimum","maximum")
          rownames(df) <- NULL
          df
        })
        datatable(sdmData(),
                  rownames = FALSE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=10)
        )
      })
      
      output$sumlayers<-DT::renderDataTable({
        summrize_rasters<-reactive({
          df<-data.frame("Layers name"=names(data$Env),"Minimum"=minValue(data$Env),"Maximum"=maxValue(data$Env))
          rownames(df) <- NULL
          df
          })
        datatable(summrize_rasters(),
                  rownames = FALSE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=10)
        )
      })
      
      txt_species_data <- paste0("The loaded dataset  consists of",code(nrow(load.occ$select)),"observations and ",code(ncol(load.occ$select)),"variables. ")
      txt_rasters_info<-paste0("You have" ,code(raster::nlayers(data$Env)),"layers.The extent is xmin=",code(raster::extent(data$Env)@xmin),",xmax=",code(raster::extent(data$Env)@xmax),",ymin=",code(raster::extent(data$Env)@ymin),",ymax=",code(raster::extent(data$Env)@ymax))
      fluidRow(
        mainPanel(width = 8, tabsetPanel(type = "tabs",
                                         tabPanel("Occurence data",
                                                  p(HTML(txt_species_data)),
                                                  dataTableOutput("occ_data_select")
                                         ),
                                         tabPanel("Summarise Occurence data",
                                                  p(HTML(txt_rasters_info)),
                                                  dataTableOutput("SpeciesTable")
                                         ),
                                         tabPanel("Layers summerize",
                                                  dataTableOutput("sumlayers"))
                                         
                                         
        ),
        id = "tabs")
        )
    
        
    })
    
   
    ###############"################################
    
    
    ###############  end function ######################
    ############# function ########
    # dynamically generate inputs (currently used in setup(biomod2)
    # shinyInput <- function(FUN, len, id, ...) {
    #   inputs = character(len)
    #   for (i in seq_len(len)) {
    #     inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...))
    #   }
    #   inputs
    # }
    ###############  end function ######################
    
    
    
    
    
    
    # output$ui_choice_species <- renderUI({
    #   #input$btn_reset_sdc # dependency so that variable-types will get updated!
    #   #out <- NULL
    #   #if (!is.null(obj$last_error)) {
    #   #  out <- list(out, fluidRow(column(12, verbatimTextOutput("ui_lasterror")), class = "wb-error-toast"))
    #   #}
    #   
    #   txt_setup <- "Select the following variables for setting up the SDC problem instance: categorical key variables, continuous key variables (optional), variables selected for PRAM (optional), sample weight (optional), hierarchical identifier (optional), variables to be deleted (optional). Also, specify the parameter alpha and set a seed at the bottom of this page."
    #   txt_setup <- paste(txt_setup, tags$br(), tags$br(), "Tip - Before you start, make sure that variable types are appropriate. If not, go to the Microdata tab and convert variables to numeric or factor.")
    #   out <- NULL
    #   out <- list(out,
    #               fluidRow(column(12, h4("Select variables", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center")),
    #               fluidRow(column(12, DT::dataTableOutput("SpeciesTable", height="100%"))))
    #   out
    # })
    
    ###########################################"
    
    ####ui correlation
    
    #glc <- GLcenfa(x = ENFA_var) 
    
    
    
    #load.occ$spec_select<-input$Pcol
    #coor$pa_dataF <- sf::st_as_sf(coor$Specdata, coords = c("lon","lat"), crs = crs(data$Env))
    #coor$Cor <- raster::extract(data$Env, coor$pa_dataF, df = TRUE)
    #coor$Cor<-coor$Cor[,-1]
    Specdata<-reactive({
      dsf<-load.occ$select
      dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
      dsf
    })
    # correlation matrix
    Z<-reactive({
      CENFA::parScale(data$Env)
    })
    
    
    # Efficient calculation of covariance matrices for Raster* objects
    mat<-reactive({
      CENFA::parCov(Z())
    })
    
    pa_data<-reactive({
      sf::st_as_sf(Specdata(), coords = c("lon","lat"), crs = crs(data$Env))
      
    })
    Cor<-reactive({
      Corr<-raster::extract(data$Env, pa_data(), df = TRUE)
      Corr<-Corr[,-1]
      Corr
    })
    
    p.mat <-reactive({
      p_mat<-ggcorrplot::cor_pmat(Cor())
      p_mat
    })
    output$ui_correlation <- renderUI({
      
      
      output$coor_mat <- DT::renderDataTable({
        datatable(mat(),
                  rownames = TRUE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=20)
        
      )})
      output$coor_plot <- renderPlot({
        ggcorrplot::ggcorrplot(mat(),ggtheme = ggplot2::theme_gray,
                               hc.order = TRUE,
                               type = "lower",
                               p.mat = p.mat(),
                               colors = c("#6D9EC1", "white", "#E46726"))
      })
      fluidRow(column(12, h4("Correlation between rasters"), align="center"),
        mainPanel(width = 8, tabsetPanel(type = "tabs",
                                          tabPanel("Correlation matrix",
                                                   DT::dataTableOutput("coor_mat")
                                          ),
                                          tabPanel("Correlation Plot",
                                                   plotOutput("coor_plot")
                                          )
                                          
                                          
        ),
        id = "tabs")
      )
    })
    
    Specdata<-reactive({
      dsf<-load.occ$select
      dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
      dsf
    })
    
    Specdata_Presence<-reactive({
      dsf<-Specdata()
      dsf<-dsf[dsf[,ncol(dsf)] == 1,]
      sp::coordinates(dsf) <-~lon+lat
      sp::proj4string(dsf) <-raster::crs(data$Env)
      dsf
    })
    
    glc<-reactive({
      GLcenfa(x = data$Env)
    })
    
    mod.enfa<-reactive({
      pr<-Specdata_Presence()
      pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
      CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
    })
 
    ############ end ui correlation
    
    ####ui enfa
    
    output$ui_enfa<-renderUI({
      glc <- glc()
      
      mod.enfa <- mod.enfa()
      if(brStick(s.factor(mod.enfa()))==1){
        
        output$enfa_scatter<-renderPlot({
          CENFA::scatter(x = mod.enfa,y = glc,n=nlayers(data$Env),p=1)
        })
      }
      else{
      observeEvent(input$number_spec,{
       
        output$enfa_scatter<-renderPlot({
          CENFA::scatter(x = mod.enfa,yax=as.numeric(input$number_spec),y = glc,n=nlayers(data$Env),p=1)
          
        })
      })
      }
      
      
      
      marg_spec<-reactive({
        mod.enfa <- mod.enfa()
        data.frame(mod.enfa@co)
      })
      
      output$marg<- DT::renderDataTable({
        datatable(marg_spec(),
                  rownames = TRUE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=20)
                  
        )})
      txt_enfa_info<-paste0('The number of significant factors is',code(brStick(s.factor(mod.enfa()))))
      fluidRow(column(12, h4("Ecological Niche Factor Analysis"), align="center"),
        mainPanel(width = 8, tabsetPanel(type = "tabs",
                                          tabPanel("ENFA",
                                                   conditionalPanel(
                                                     condition = brStick(s.factor(mod.enfa()))>1,
                                                     selectInput('number_spec', 'Please select the number between 2 and the number of significant factors.', 2:brStick(s.factor(mod.enfa())), multiple = FALSE, selectize = TRUE)
                                                   ),
                                                   plotOutput("enfa_scatter"))
                                          ,
                                          tabPanel("Marginality and specialization",
                                                   p(HTML(txt_enfa_info)),
                                                   DT::dataTableOutput("marg")
                                          )),
                  id = "tabs")
        
      )
    })
    
    
    output$ui_preparation_main <- renderUI({
      out <- NULL
      val <- obj$cur_selection_results
      if (val=="btn_preparation_results_1") {
        return(uiOutput("ui_view_species_data"))
      }
      if (val=="btn_preparation_results_2") {
        return(uiOutput("ui_correlation"))
      }
      if (val=="btn_preparation_results_3") {
        return(uiOutput("ui_enfa"))
      }
      
    })	
    
    output$ui_preparation_sidebar_left <- renderUI({
      output$ui_sel_preparation_btns <- renderUI({
        cc1 <- c("Summarise")
        cc2 <- c("Correlation", "ENFA")
        df <- data.frame(lab=c(cc1,cc2), header=NA)
        df$header[1] <- "View"
        df$header[2] <- "Spatial Analysis"
        out <- NULL
        for (i in 1:nrow(df)) {
          id <- paste0("btn_preparation_results_",i)
          if (obj$cur_selection_results==id) {
            style <- "primary"
          } else {
            style <- "default"
          }
          if (!is.na(df$header[i])) {
            out <- list(out, fluidRow(column(12, h4(df$header[i]), align="center")))
          }
          out <- list(out, fluidRow(
            column(12, bsButton(id, label=df$lab[i], block=TRUE, size="extra-small", style=style))
          ))
        }
        out
      })
      # required observers that update the color of the active button!
      eval(parse(text=genObserver_menus(pat="btn_preparation_results_", n=1:3, updateVal="cur_selection_results")))
      return(uiOutput("ui_sel_preparation_btns"))
    })
    output$ui_anonymize_noproblem <- renderUI({
      return(list(
        noInputData(uri="ui_preparation"),
        fluidRow(column(12, tags$br(), p(""), align="center"))
        #fluidRow(column(12, myActionButton("nodata_anonymize_uploadproblem", label="Upload a previously saved problem", btn.style="primary"), align="center"))
      ))
    })
    output$ui_preparation <- renderUI({
      if(length(input$Occ)==0){
        return(uiOutput("ui_anonymize_noproblem"))}
      else{
      fluidRow(
        column(2, uiOutput("ui_preparation_sidebar_left"), class="wb_sidebar"),
        column(10, uiOutput("ui_preparation_main"), class="wb-maincolumn"))
      }
    }
    )
    
    #########################################
    ####################################"
    #BlockCV contents
    # loading the package
    library(blockCV)	
    library(raster)
    library(sf)	
    
    sac<-reactive({
      a = try(withProgress(message = 'Variables loading',
                           spatialAutoRange(rasterLayer = data$Env, 
                                            doParallel = T, 
                                            plotVariograms = TRUE,
                                            showPlots = FALSE)))
      a
    })
    
    range<-reactive({
      sac<-sac()
      round(sac$range,0)
      
    })
    output$ui_spatial_auto_range<-renderUI({
      
      
      tableRange<-reactive({
        sac<-sac()
        sac$rangeTable
      })
      
      
      
      
      output$tableRange <- DT::renderDataTable({
        datatable(tableRange(),
                  rownames = FALSE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=20)

        )})
      observeEvent(input$vario_var,{
        output$variogram<-renderPlot({
          sac<-sac()
          vect<-names(data$Env)
          plot(sac$variograms[[which(vect==input$vario_var)]])
        }) 
      })
      
      output$barchart <- renderPlot({
        sac<-sac()
        sac$plots$barchart
      })
      
      output$mapplot <- renderPlot({
        sac<-sac()
        sac$plots$mapplot
      })
      
      fluidRow(column(12, h4("Spatial"), align="center"),
               mainPanel(width = 8, tabsetPanel(type = "tabs",
                                                tabPanel("barchart",
                                                         p('Spatial autocorrelation ranges in input covariates'),
                                                         plotOutput("barchart")),
                                                tabPanel("mapplot",
                                                         p('Corresponding spatial blocks (the selected block size is based on median spatial autocorrelation range across all input data)'),
                                                         plotOutput("mapplot")),
                                                tabPanel("Autocorrelation range table",
                                                         p('Spatial autocorrelation ranges table in input covariates'),
                                                         DT::dataTableOutput("tableRange")),
                                                tabPanel("Variogramme",
                                                         selectInput('vario_var', 'Please select the predictor to see variogram corresponding', names(data$Env), multiple = FALSE, selectize = TRUE),
                                                         plotOutput("variogram"))
                                                         
                                                ),
                         id = "tabs")
               
      )
    })
    output$ui_spatial_blocks<-renderUI({
      observeEvent(input$number_fold,{
        load.occ$k<-input$number_fold
      })
      
      observeEvent(input$allocation_fold,{
        
        load.occ$allocation_fold<-input$allocation_fold
      })
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      pa_data<-reactive({
        load.occ$pa_data<-sf::st_as_sf(Specdata(), coords = c("lon","lat"), crs = crs(data$Env))
        load.occ$pa_data
      })
      spatialblock<-reactive({
        a = try(withProgress(message = 'Variables loading',
                             spatialBlock(speciesData = pa_data(),
                                          species = load.occ$spec_select,
                                          rasterLayer = data$Env,
                                          theRange = range(), #load.occ$range, # size of the blocks
                                          k = load.occ$k,
                                          showBlocks = TRUE,
                                          selection = load.occ$allocation_fold,
                                          iteration = 100, # find evenly dispersed folds
                                          biomod2Format = FALSE,
                                          xOffset = 0, # shift the blocks horizontally
                                          yOffset = 0)))
        a
        
      })
      
      output$sp_block<-renderPlot({
        spatialblock<-spatialblock()
        spatialblock$plots + geom_sf(data = pa_data(), alpha = 0.5)
      }) 
      
  
      output$sum_fold <- DT::renderDataTable({
        spatialblock<-spatialblock()
        sumfold<-summarise_fold(spatialblock)
        datatable(sumfold,
                  rownames = FALSE,
                  selection="none",
                  options = list(scrollX=TRUE, scrollY=250, lengthMenu=list(c(20, 50, 100, -1), c('20', '50', '100', 'All')), pageLength=20)
                  
        )})
      observeEvent(input$test_fold,{
        load.occ$fold<-input$test_fold
      })
      output$test_train_plot<-renderPlot({
        spatialblock<-spatialblock()
        Explorer(spatialblock, data$Env, pa_data(),load.occ$fold)
      }) 
      
      fluidRow(column(12, h4("Spatial blocking"), align="center"),
               mainPanel(width = 8, tabsetPanel(type = "tabs",
                                                tabPanel("Spatial blocking",
                                                         p('Set spatial bloking parameters'),
                                                         sliderInput("number_fold", "folds", min=1, max=100, value=5),
                                                         selectInput("allocation_fold","allocation of blocks to folds",choices = c("random","systematic"),selected="random"),
                                                         sliderInput("test_fold","Select the number of fold to assign as test dataset",min = 1,max=100,value = 1),
                                                         plotOutput("sp_block"),
                                                         plotOutput("test_train_plot")),
                                                tabPanel("Summarize fold",
                                                         p('Fold summarizing. Purcentage means the purcentage of test dataset'),
                                                         DT::dataTableOutput("sum_fold"))
                                                
               ),
               id = "tabs")
               
      )
      # out <- NULL
      # out <- list(out,
      #             fluidRow(box(title = "Set spatial bloking parameters",
      #                      sliderInput("number_fold", "folds", min=1, max=100, value=5),
      #                      selectInput("allocation_fold","allocation of blocks to folds",choices = c("random","systematic"),selected="random")),
      #                      box(title = "Number fold of test dataset",
      #                          sliderInput("test_fold","Select the number of fold to assign as test dataset",min = 1,max=5,value = 1))) #max=input$number_fold
      # )
      # out<-list(out,
      #           fluidRow(box(title = 'Spatial blocking',
      #                        plotOutput("sp_block")),
      #                    box(title = 'Summarise fold',
      #                        DT::dataTableOutput("sum_fold"))))
      # out<-list(out,
      #           fluidRow(column(12, plotOutput("test_train_plot"))))
      # 
      # out
    })
    
    observeEvent(input$number_no_block_fold,{
      load.occ$number_no_block_fold<-input$number_no_block_fold
    })
    
    Specdata<-reactive({
      dsf<-load.occ$select
      dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
      dsf
    })
    
    kfold<-reactive({
      Specdata<-Specdata()
      set.seed(1994)
      kfold<-dismo::kfold(Specdata,load.occ$number_no_block_fold)
      kfold
    })

        output$ui_blockCV_main <- renderUI({
      out <- NULL
      val <- obj$cur_selection_results
      ## Categorical (defined in controller/ui_results_imputation.R)
      if (val=="btn_impute_results_1") {
        return(uiOutput("ui_spatial_auto_range"))
      }
      if (val=="btn_impute_results_2") {
        return(uiOutput("ui_spatial_blocks"))
      }
    })	
    
    output$ui_blockCV_sidebar_left <- renderUI({
      output$ui_sel_impute_btns <- renderUI({
        cc1 <- c("spatial autocorrelation", "spatial blocks")
        df <- data.frame(lab=c(cc1), header=NA)
        df$header[1] <- "View"
         out <- NULL
        for (i in 1:nrow(df)) {
          id <- paste0("btn_impute_results_",i)
          if (obj$cur_selection_results==id) {
            style <- "primary"
          } else {
            style <- "default"
          }
          if (!is.na(df$header[i])) {
            out <- list(out, fluidRow(column(12, h4(df$header[i]), align="center")))
          }
          out <- list(out, fluidRow(
            column(12, bsButton(id, label=df$lab[i], block=TRUE, size="extra-small", style=style))
          ))
        }
        out
      })
      # required observers that update the color of the active button!
      eval(parse(text=genObserver_menus(pat="btn_impute_results_", n=1:2, updateVal="cur_selection_results")))
      return(uiOutput("ui_sel_impute_btns"))
    })
    
    output$ui_blockCV_noproblem <- renderUI({
      return(list(
        noInputData(uri="ui_blockCV"),
        fluidRow(column(12, tags$br(), p(""), align="center"))
        #fluidRow(column(12, myActionButton("nodata_anonymize_uploadproblem", label="Upload a previously saved problem", btn.style="primary"), align="center"))
      ))
    })
    output$ui_blockCV <- renderUI({
      if(length(input$Occ)==0){
        return(uiOutput("ui_blockCV_noproblem"))}
      else{
      fluidRow(
        column(2, uiOutput("ui_blockCV_sidebar_left"), class="wb_sidebar"),
        column(10, uiOutput("ui_blockCV_main"), class="wb-maincolumn"))
      }
    }
    )
    ########################################### End BlockCV contents ##############
    ##########################################################################
    #Biomod2 contents
    ##########################################################################
    
    ##### Maxent content
    output$ui_MaxEnt<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$MaxEnt,{
        validate(
          need(length(input$var_expl) > 0, 'Choose specie predictors first !')
        )
        
          data$enfa<-raster::subset(data$Env,input$var_expl)
          Specdata<-Specdata()
          set.seed(1994)
          fold<-dismo::kfold(Specdata,input$number_no_block_fold)
          #fold<-kfold()
          model<-list()
          evaluate_model<-list()
          for (i in 1:5) {
           p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
            a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
            #test<-Specdata[fold == i, ] 
            
            occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
            
            bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
            model[[i]] <- dismo::maxent(data$enfa, p, a) #, factors='Sol'
            evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
            
          }
          model_pred<-list()
          auc <- sapply(evaluate_model, function(x){x@auc})
          model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
          model_pred[["AUC"]]<-auc[which.max(auc)]
          model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
          model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
          model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
          observeEvent(input$probaplot,{
            if(input$probaplot=='Probability of occurence(absence/presence)'){
              title_probaplot<-'Probability of occurence(absence/presence)'
              map<-model_pred[["espece"]]}
            if(input$probaplot=='Presence/Absence'){
              title_probaplot<-'Presence/Absence'
              map<-model_pred[["PresenceAbsence"]]
              
            }
            if(input$probaplot=='Probability of occurence(presence)'){
              title_probaplot<-'Probability of occurence(presence)'
              map<-model_pred[["ProbaPresence"]]
            }
            output$proba_occ<-renderPlot({
              if(title_probaplot=='Presence/Absence'){PASpecies(map)}
              else{
                ggR_P(map) 
              }
              
            })
          })
          observeEvent(input$model_ev,{
            if(input$model_ev == 'ROC') {ev<-'ROC'}
            if(input$model_ev == 'density') {ev<-'density'}
            if(input$model_ev == 'boxplot') {ev<-'boxplot'}
            if(input$model_ev == 'kappa') {ev<-'kappa'}
            if(input$model_ev == 'FPR') {ev<-'FPR'}
            if(input$model_ev == 'prevalence') {ev<-'prevalence'}
            output$eval<-renderPlot({
              if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
              else{
                if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
                else{
                  plot(evaluate_model[[which.max(auc)]],ev)
                }
              }
              
              
            })
          })
          observeEvent(input$response_var,{
            output$response_eco<-renderPlot({
              dismo::response(model[[which.max(auc)]],var=input$response_var,main=load.occ$spec_select)
            })
          })
          output$var_importance<-renderPlot({
            plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
          })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("MaxEnt"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors according to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      
      out<-list(out,
                sidebarPanel(
                  selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                              c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                  ),
                  conditionalPanel(
                    condition = "input.choice_block == 'Modelling without spatial blocking'",
                    sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
                  )),
                
                mainPanel(width = 6, tabsetPanel(type = "tabs",
                                                 tabPanel("Specie predictors",
                                                          selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                          myActionButton("MaxEnt",label=("Apply MaxEnt"), "primary"),
                                                          plotOutput("enfa_var")
                                                 ),
                                                 tabPanel("Map",
                                                          
                                                          selectInput('probaplot', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("proba_occ")
                                                          
                                                 ),
                                                 tabPanel("Model Evaluation",
                                                          selectInput('model_ev', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("eval")
                                                 ),
                                                 tabPanel("Variable response",
                                                          selectInput('response_var', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("response_eco")
                                                 ),
                                                 tabPanel("Variable Importance",
                                                          plotOutput("var_importance")
                                                 )
                                                 
                                                 
                ),
                id = "tabs")
      )
      out
      
    })
 ###########"end Maxent ######
    
 #### bioclim contents
  
    output$ui_bioclim<-renderUI({
  
  Specdata<-reactive({
    dsf<-load.occ$select
    dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
    dsf
  })
  
  Specdata_Presence<-reactive({
    dsf<-Specdata()
    dsf<-dsf[dsf[,ncol(dsf)] == 1,]
    sp::coordinates(dsf) <-~lon+lat
    sp::proj4string(dsf) <-raster::crs(data$Env)
    dsf
  })
  
  glc<-reactive({
    GLcenfa(x = data$Env)
  })
  
  mod.enfa<-reactive({
    pr<-Specdata_Presence()
    pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
    CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
  })
  enfa_plot<-reactive({
    glc <- glc()
    
    mod.enfa <- mod.enfa()
    CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
  })
  output$enfa_var<-renderPlot({
    enfa_plot()
  })
  
  observeEvent(input$Bioclim,{
    validate(
      need(length(input$var_expl_Bioclim) > 0, 'Choose specie predictors first !')
    )
    
    data$enfa<-raster::subset(data$Env,input$var_expl_Bioclim)
    Specdata<-Specdata()
    set.seed(1994)
    fold<-dismo::kfold(Specdata,input$number_no_block_fold_bioclim)
    #fold<-kfold()
    model<-list()
    evaluate_model<-list()
    for (i in 1:5) {
      
      # testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
      # testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
      #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
      
      p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      #test<-Specdata[fold == i, ] 
      
      occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      
      bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      model[[i]] <- dismo::bioclim(data$enfa, p) #, factors='Sol'
      evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
      
    }
    model_pred<-list()
    auc <- sapply(evaluate_model, function(x){x@auc})
    model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
    model_pred[["AUC"]]<-auc[which.max(auc)]
    model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
    model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
    model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
    observeEvent(input$probaplot_Bioclim,{
      if(input$probaplot_Bioclim=='Probability of occurence(absence/presence)'){
        title_probaplot_Bioclim<-'Probability of occurence(absence/presence)'
        map<-model_pred[["espece"]]}
      if(input$probaplot_Bioclim=='Presence/Absence'){
        title_probaplot_Bioclim<-'Presence/Absence'
        map<-model_pred[["PresenceAbsence"]]
        
      }
      if(input$probaplot_Bioclim=='Probability of occurence(presence)'){
        title_probaplot_Bioclim<-'Probability of occurence(presence)'
        map<-model_pred[["ProbaPresence"]]
      }
      output$proba_occ_Bioclim<-renderPlot({
        if(title_probaplot_Bioclim=='Presence/Absence'){PASpecies(map)}
        else{
          ggR_P(map) 
        }
        
      })
    })
    observeEvent(input$model_ev_Bioclim,{
      if(input$model_ev_Bioclim == 'ROC') {ev<-'ROC'}
      if(input$model_ev_Bioclim == 'density') {ev<-'density'}
      if(input$model_ev_Bioclim == 'boxplot') {ev<-'boxplot'}
      if(input$model_ev_Bioclim == 'kappa') {ev<-'kappa'}
      if(input$model_ev_Bioclim == 'FPR') {ev<-'FPR'}
      if(input$model_ev_Bioclim == 'prevalence') {ev<-'prevalence'}
      output$eval_Bioclim<-renderPlot({
        if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
        else{
          if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
          else{
            plot(evaluate_model[[which.max(auc)]],ev)
          }
        }
        
        
      })
    })
    observeEvent(input$response_var_Bioclim,{
      output$response_eco<-renderPlot({
        dismo::response(model[[which.max(auc)]],var=input$response_var_Bioclim,main=load.occ$spec_select)
      })
    })
    output$var_importance_Bioclim<-renderPlot({
      plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
    })
  })
  
  out <- NULL
  txt_setup<-'The Bioclim software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Bioclim is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
  out <- fluidRow(
    column(width = 12, offset = 0, h3("Bioclim"), class="wb-header"),
    column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Bioclim method."), class="wb-header-hint"),
    fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
  )
  out<-list(out,
            sidebarPanel(
              selectInput("choice_block_bioclim", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                          c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
              ),
              conditionalPanel(
                condition = "input.choice_block_bioclim == 'Modelling without spatial blocking'",
                sliderInput("number_no_block_fold_bioclim", "Please set the number of fold", min = 1, max = 100, value = 5)
              )),
              
            mainPanel(width = 6, tabsetPanel(type = "tabs",
                                             tabPanel("Specie predictors",
                                                      selectInput('var_expl_Bioclim', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                      myActionButton("Bioclim",label=("Apply Bioclim"), "primary"),
                                                      plotOutput("enfa_var")
                                                      ),
                                             tabPanel("Map",
                                                      
                                                      selectInput('probaplot_Bioclim', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("proba_occ_Bioclim")
                                                      
                                             ),
                                             tabPanel("Model Evaluation",
                                                      selectInput('model_ev_Bioclim', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("eval_Bioclim")
                                                      ),
                                             tabPanel("Variable response",
                                                      selectInput('response_var_Bioclim', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("response_eco")
                                                      ),
                                             tabPanel("Variable Importance",
                                                      plotOutput("var_importance_Bioclim")
                                                      )
                                             
                                             
            ),
            id = "tabs")
  )
  out
  
  
})
 ##### end bioclim
    
    #### domain contents ###
    
    output$ui_domain<-renderUI({
  
  Specdata<-reactive({
    dsf<-load.occ$select
    dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
    dsf
  })
  
  Specdata_Presence<-reactive({
    dsf<-Specdata()
    dsf<-dsf[dsf[,ncol(dsf)] == 1,]
    sp::coordinates(dsf) <-~lon+lat
    sp::proj4string(dsf) <-raster::crs(data$Env)
    dsf
  })
  
  glc<-reactive({
    GLcenfa(x = data$Env)
  })
  
  mod.enfa<-reactive({
    pr<-Specdata_Presence()
    pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
    CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
  })
  enfa_plot<-reactive({
    glc <- glc()
    
    mod.enfa <- mod.enfa()
    CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
  })
  output$enfa_var<-renderPlot({
    enfa_plot()
  })
  
  observeEvent(input$Domain,{
    validate(
      need(length(input$var_expl_Domain) > 0, 'Choose specie predictors first !')
    )
    
    data$enfa<-raster::subset(data$Env,input$var_expl_Domain)
    Specdata<-Specdata()
    set.seed(1994)
    fold<-dismo::kfold(Specdata,input$number_no_block_fold_Domain)
    #fold<-kfold()
    model<-list()
    evaluate_model<-list()
    for (i in 1:5) {
      p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      
      bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      model[[i]] <- dismo::domain(data$enfa, p) #, factors='Sol'
      evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
      
    }
    model_pred<-list()
    auc <- sapply(evaluate_model, function(x){x@auc})
    model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
    model_pred[["AUC"]]<-auc[which.max(auc)]
    model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
    model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
    model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
    
    observeEvent(input$probaplot_Domain,{
      if(input$probaplot_Domain=='Probability of occurence(absence/presence)'){
        title_probaplot_Domain<-'Probability of occurence(absence/presence)'
        map<-model_pred[["espece"]]}
      if(input$probaplot_Domain=='Presence/Absence'){
        title_probaplot_Domain<-'Presence/Absence'
        map<-model_pred[["PresenceAbsence"]]
        
      }
      if(input$probaplot_Domain=='Probability of occurence(presence)'){
        title_probaplot_Domain<-'Probability of occurence(presence)'
        map<-model_pred[["ProbaPresence"]]
      }
      output$proba_occ_Domain<-renderPlot({
        if(title_probaplot_Domain=='Presence/Absence'){PASpecies(map)}
        else{
          ggR_P(map) 
        }
        
      })
    })
    observeEvent(input$model_ev_Domain,{
      if(input$model_ev_Domain == 'ROC') {ev<-'ROC'}
      if(input$model_ev_Domain == 'density') {ev<-'density'}
      if(input$model_ev_Domain == 'boxplot') {ev<-'boxplot'}
      if(input$model_ev_Domain == 'kappa') {ev<-'kappa'}
      if(input$model_ev_Domain == 'FPR') {ev<-'FPR'}
      if(input$model_ev_Domain == 'prevalence') {ev<-'prevalence'}
      output$eval_Domain<-renderPlot({
        if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
        else{
          if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
          else{
            plot(evaluate_model[[which.max(auc)]],ev)
          }
        }
        
        
      })
    })
    observeEvent(input$response_var_Domain,{
      output$response_eco_Domain<-renderPlot({
        dismo::response(model[[which.max(auc)]],var=input$response_var_Domain,main=load.occ$spec_select)
      })
    })
    output$var_importance_Domain<-renderPlot({
      plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
    })
  })
  
  
  out <- NULL
  txt_setup<-'The Domain algorithm computes the Gower distance between environmental variables at any location and those at any of the known locations of occurrence (training sites).'
  out <- fluidRow(
    column(width = 12, offset = 0, h3("Domain"), class="wb-header"),
    column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Domain method."), class="wb-header-hint"),
    fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
  )
  out<-list(out,
            sidebarPanel(
              selectInput("choice_block_Domain", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                          c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
              ),
              conditionalPanel(
                condition = "input.choice_block_Domain == 'Modelling without spatial blocking'",
                sliderInput("number_no_block_fold_Domain", "Please set the number of fold", min = 1, max = 100, value = 5)
              )),
            
            mainPanel(width = 6, tabsetPanel(type = "tabs",
                                             tabPanel("Specie predictors",
                                                      selectInput('var_expl_Domain', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                      myActionButton("Domain",label=("Apply Domain"), "primary"),
                                                      plotOutput("enfa_var")
                                             ),
                                             tabPanel("Map",
                                                      
                                                      selectInput('probaplot_Domain', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("proba_occ_Domain")
                                                      
                                             ),
                                             tabPanel("Model Evaluation",
                                                      selectInput('model_ev_Domain', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("eval_Domain")
                                             ),
                                             tabPanel("Variable response",
                                                      selectInput('response_var_Domain', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("response_eco_Domain")
                                             ),
                                             tabPanel("Variable Importance",
                                                      plotOutput("var_importance_Domain")
                                             )
                                             
                                             
            ),
            id = "tabs")
  )
  out
})
    #### end domain ###
    
    # ###Maxent contents
    # output$ui_Maxent<-renderUI({
    #   
    #   Specdata<-reactive({
    #     dsf<-load.occ$select
    #     dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
    #     dsf
    #   })
    #   
    #   Specdata_Presence<-reactive({
    #     dsf<-Specdata()
    #     dsf<-dsf[dsf[,ncol(dsf)] == 1,]
    #     sp::coordinates(dsf) <-~lon+lat
    #     sp::proj4string(dsf) <-raster::crs(data$Env)
    #     dsf
    #   })
    #   
    #   glc<-reactive({
    #     GLcenfa(x = data$Env)
    #   })
    #   
    #   mod.enfa<-reactive({
    #     pr<-Specdata_Presence()
    #     pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
    #     CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
    #   })
    #   enfa_plot<-reactive({
    #     glc <- glc()
    #     
    #     mod.enfa <- mod.enfa()
    #     CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
    #   })
    #   output$enfa_var<-renderPlot({
    #     enfa_plot()
    #   })
    #   
    #   observeEvent(input$MaxEnt,{
    #     validate(
    #       need(length(input$var_expl) > 0, 'Choose specie predictors first !')
    #     )
    #     
    #     data$enfa<-raster::subset(data$Env,input$var_expl)
    #     Specdata<-Specdata()
    #     set.seed(1994)
    #     fold<-dismo::kfold(Specdata,input$number_no_block_fold)
    #     #fold<-kfold()
    #     model<-list()
    #     evaluate_model<-list()
    #     for (i in 1:5) {
    #       p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
    #       a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
    #       #test<-Specdata[fold == i, ] 
    #       
    #       occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
    #       
    #       bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
    #       model[[i]] <- dismo::mahal(data$enfa, p) #, factors='Sol'
    #       evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
    #       
    #     }
    #     model_pred<-list()
    #     auc <- sapply(evaluate_model, function(x){x@auc})
    #     model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
    #     model_pred[["AUC"]]<-auc[which.max(auc)]
    #     model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
    #     model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
    #     model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
    #     observeEvent(input$probaplot,{
    #       if(input$probaplot=='Probability of occurence(absence/presence)'){
    #         title_probaplot<-'Probability of occurence(absence/presence)'
    #         map<-model_pred[["espece"]]}
    #       if(input$probaplot=='Presence/Absence'){
    #         title_probaplot<-'Presence/Absence'
    #         map<-model_pred[["PresenceAbsence"]]
    #         
    #       }
    #       if(input$probaplot=='Probability of occurence(presence)'){
    #         title_probaplot<-'Probability of occurence(presence)'
    #         map<-model_pred[["ProbaPresence"]]
    #       }
    #       output$proba_occ<-renderPlot({
    #         if(title_probaplot=='Presence/Absence'){PASpecies(map)}
    #         else{
    #           ggR_P(map) 
    #         }
    #         
    #       })
    #     })
    #     observeEvent(input$model_ev,{
    #       if(input$model_ev == 'ROC') {ev<-'ROC'}
    #       if(input$model_ev == 'density') {ev<-'density'}
    #       if(input$model_ev == 'boxplot') {ev<-'boxplot'}
    #       if(input$model_ev == 'kappa') {ev<-'kappa'}
    #       if(input$model_ev == 'FPR') {ev<-'FPR'}
    #       if(input$model_ev == 'prevalence') {ev<-'prevalence'}
    #       output$eval<-renderPlot({
    #         if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
    #         else{
    #           if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
    #           else{
    #             plot(evaluate_model[[which.max(auc)]],ev)
    #           }
    #         }
    #         
    #         
    #       })
    #     })
    #     observeEvent(input$response_var,{
    #       output$response_eco<-renderPlot({
    #         dismo::response(model[[which.max(auc)]],var=input$response_var,main=load.occ$spec_select)
    #       })
    #     })
    #     output$var_importance<-renderPlot({
    #       plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
    #     })
    #   })
    #   
    #   
    #   
    #   
    #   out <- NULL
    #   txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
    #   out <- fluidRow(
    #     column(width = 12, offset = 0, h3("Maxent"), class="wb-header"),
    #     column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
    #     fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
    #   )
    #   out <- list(out,
    #               column(12,offset=0,sidebarPanel(
    #                 selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
    #                             c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
    #                 ),
    #                 conditionalPanel(
    #                   condition = "input.choice_block == 'Modelling without spatial blocking'",
    #                   sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
    #                 )),align="center"))
    #   out <- list(out,
    #               fluidRow(
    #                 box(title = "Specie predictors",
    #                     selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE)),
    #                 box(title = "Ecological Niche Factor Analysis",
    #                     plotOutput("enfa_var"))))
    #   
    #   out <- list(out,
    #               fluidRow(actionButton('MaxEnt', 'Apply')))
    #   out <- list(out,
    #               fluidRow(
    #                 box(title='Probability of occurence',
    #                     selectInput('probaplot', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
    #                     plotOutput("proba_occ")),
    #                 box(title = "Model Evaluation",
    #                     selectInput('model_ev', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
    #                     plotOutput("eval")),
    #                 box(title = 'Variable response',
    #                     selectInput('response_var', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
    #                     plotOutput("response_eco")),
    #                 box(title = 'Variable importance',
    #                     plotOutput("var_importance"))))
    #   out
    #   
    # })
    
    output$ui_mahal<-renderUI({
  
  Specdata<-reactive({
    dsf<-load.occ$select
    dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
    dsf
  })
  
  Specdata_Presence<-reactive({
    dsf<-Specdata()
    dsf<-dsf[dsf[,ncol(dsf)] == 1,]
    sp::coordinates(dsf) <-~lon+lat
    sp::proj4string(dsf) <-raster::crs(data$Env)
    dsf
  })
  
  glc<-reactive({
    GLcenfa(x = data$Env)
  })
  
  mod.enfa<-reactive({
    pr<-Specdata_Presence()
    pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
    CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
  })
  enfa_plot<-reactive({
    glc <- glc()
    
    mod.enfa <- mod.enfa()
    CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
  })
  output$enfa_var<-renderPlot({
    enfa_plot()
  })
  
  observeEvent(input$Mahal,{
    validate(
      need(length(input$var_expl_Mahal) > 0, 'Choose specie predictors first !')
    )
    
    data$enfa<-raster::subset(data$Env,input$var_expl_Mahal)
    Specdata<-Specdata()
    set.seed(1994)
    fold<-dismo::kfold(Specdata,input$number_no_block_fold_Mahal)
    #fold<-kfold()
    model<-list()
    evaluate_model<-list()
    for (i in 1:5) {
      
      # testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
      # testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
      #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
      
      p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      #test<-Specdata[fold == i, ] 
      
      occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
      
      bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
      model[[i]] <- dismo::mahal(data$enfa, p) #, factors='Sol'
      evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
      
    }
    model_pred<-list()
    auc <- sapply(evaluate_model, function(x){x@auc})
    model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
    model_pred[["AUC"]]<-auc[which.max(auc)]
    model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
    model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
    model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
    observeEvent(input$probaplot_Mahal,{
      if(input$probaplot_Mahal=='Probability of occurence(absence/presence)'){
        title_probaplot_Mahal<-'Probability of occurence(absence/presence)'
        map<-model_pred[["espece"]]}
      if(input$probaplot_Mahal=='Presence/Absence'){
        title_probaplot_Mahal<-'Presence/Absence'
        map<-model_pred[["PresenceAbsence"]]
        
      }
      if(input$probaplot_Mahal=='Probability of occurence(presence)'){
        title_probaplot_Mahal<-'Probability of occurence(presence)'
        map<-model_pred[["ProbaPresence"]]
      }
      output$proba_occ_Mahal<-renderPlot({
        if(title_probaplot_Mahal=='Presence/Absence'){PASpecies(map)}
        else{
          ggR_P(map) 
        }
        
      })
    })
    observeEvent(input$model_ev_Mahal,{
      if(input$model_ev_Mahal == 'ROC') {ev<-'ROC'}
      if(input$model_ev_Mahal == 'density') {ev<-'density'}
      if(input$model_ev_Mahal == 'boxplot') {ev<-'boxplot'}
      if(input$model_ev_Mahal == 'kappa') {ev<-'kappa'}
      if(input$model_ev_Mahal == 'FPR') {ev<-'FPR'}
      if(input$model_ev_Mahal == 'prevalence') {ev<-'prevalence'}
      output$eval_Mahal<-renderPlot({
        if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
        else{
          if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
          else{
            plot(evaluate_model[[which.max(auc)]],ev)
          }
        }
        
        
      })
    })
    observeEvent(input$response_var_Mahal,{
      output$response_eco_Mahal<-renderPlot({
        dismo::response(model[[which.max(auc)]],var=input$response_var_Mahal,main=load.occ$spec_select)
      })
    })
    output$var_importance_Mahal<-renderPlot({
      plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
    })
  })
  
  
  
  
  out <- NULL
  txt_setup<-'The Mahal software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Mahal is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
  out <- fluidRow(
    column(width = 12, offset = 0, h3("Mahal"), class="wb-header"),
    column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Mahal method."), class="wb-header-hint"),
    fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
  )
  out<-list(out,
            sidebarPanel(
              selectInput("choice_block_Mahal", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                          c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
              ),
              conditionalPanel(
                condition = "input.choice_block_Mahal == 'Modelling without spatial blocking'",
                sliderInput("number_no_block_fold_Mahal", "Please set the number of fold", min = 1, max = 100, value = 5)
              )),
            
            mainPanel(width = 6, tabsetPanel(type = "tabs",
                                             tabPanel("Specie predictors",
                                                      selectInput('var_expl_Mahal', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                      myActionButton("Mahal",label=("Apply Mahal"), "primary"),
                                                      plotOutput("enfa_var")
                                             ),
                                             tabPanel("Map",
                                                      
                                                      selectInput('probaplot_Mahal', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("proba_occ_Mahal")
                                                      
                                             ),
                                             tabPanel("Model Evaluation",
                                                      selectInput('model_ev_Mahal', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("eval_Mahal")
                                             ),
                                             tabPanel("Variable response",
                                                      selectInput('response_var_Mahal', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                      plotOutput("response_eco_Mahal")
                                             ),
                                             tabPanel("Variable Importance",
                                                      plotOutput("var_importance_Mahal")
                                             )
                                             
                                             
            ),
            id = "tabs")
  )
  out
  # out <- list(out,
  #             column(12,offset=0,sidebarPanel(
  #               selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
  #                           c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
  #               ),
  #               conditionalPanel(
  #                 condition = "input.choice_block == 'Modelling without spatial blocking'",
  #                 sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
  #               )),align="center"))
  # out <- list(out,
  #             fluidRow(
  #               box(title = "Specie predictors",
  #                   selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE)),
  #               box(title = "Ecological Niche Factor Analysis",
  #                   plotOutput("enfa_var"))))
  # 
  # out <- list(out,
  #             fluidRow(actionButton('Mahal', 'Apply')))
  # out <- list(out,
  #             fluidRow(
  #               box(title='Probability of occurence',
  #                   selectInput('probaplot_Mahal', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
  #                   plotOutput("proba_occ_Mahal")),
  #               box(title = "Model Evaluation",
  #                   selectInput('model_ev_Mahal', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
  #                   plotOutput("eval")),
  #               box(title = 'Variable response',
  #                   selectInput('response_var_Mahal', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
  #                   plotOutput("response_eco_Mahal")),
  #               box(title = 'Variable importance',
  #                   plotOutput("var_importance_Mahal"))))
  # out
  
})
    ### end mahal
    
    ### GLM contents
    output$ui_GLM<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$GLM,{
        validate(
          need(length(input$var_expl_GLM) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl_GLM)
        pa_data<-reactive({
          pa_data<-sf::st_as_sf(Specdata(), coords = c("lon","lat"), crs = crs(data$enfa))
          pa_data
        })
        # extract the raster values for the species points as a dataframe
        mydataF<-reactive({
          mydataF <- raster::extract(data$enfa, pa_data(), df = TRUE)
          mydataF <- mydataF[,-1]
          Specdata<-Specdata()
          spec<-Specdata[,load.occ$spec_select]
          mydataF<-cbind(mydataF,spec)
          mydataF<-mydataF %>% dplyr::rename(Species=spec)
          mydataF
        })
        mydataF<-mydataF()
        set.seed(1994)
        fold<-dismo::kfold(Specdata(),input$number_no_block_fold_GLM)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {

          testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          #dsf<-dsf[dsf[,ncol(dsf)] == 1,]


          model[[i]] <- glm(as.factor(Species)~., mydataF[fold != i, ], na.action=na.omit,
                            family = binomial(link = "logit"))
          evaluate_model[[i]] <- dismo::evaluate(testpres, testbackg, model[[i]])

        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        observeEvent(input$probaplot_GLM,{
          if(input$probaplot_GLM=='Probability of occurence(absence/presence)'){
            title_probaplot_GLM<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot_GLM=='Presence/Absence'){
            title_probaplot_GLM<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]

          }
        if(input$probaplot_GLM=='Probability of occurence(presence)'){
          title_probaplot_GLM<-'Probability of occurence(presence)'
          map<-model_pred[["ProbaPresence"]]
        }
        output$proba_occ_GLM<-renderPlot({
          if(title_probaplot_GLM=='Presence/Absence'){PASpecies(map)}
          else{
            ggR_P(map)
          }

        })
        
        })
        observeEvent(input$model_ev_GLM,{
          if(input$model_ev_GLM == 'ROC') {ev_GLM<-'ROC'}
          if(input$model_ev_GLM == 'density') {ev_GLM<-'density'}
          if(input$model_ev_GLM == 'boxplot') {ev_GLM<-'boxplot'}
          if(input$model_ev_GLM == 'kappa') {ev_GLM<-'kappa'}
          if(input$model_ev_GLM == 'FPR') {ev_GLM<-'FPR'}
          if(input$model_ev_GLM == 'prevalence') {ev_GLM<-'prevalence'}
          output$eval_GLM<-renderPlot({
            if(ev_GLM=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev_GLM=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev_GLM)
              }
            }


          })
        })
        observeEvent(input$response_var_GLM,{
          output$response_eco_GLM<-renderPlot({
            dismo::response(model[[which.max(auc)]])#,var=input$response_var,main=load.occ$spec_select)
          })
        })
        output$var_importance_GLM<-renderPlot({
          plot(model[[which.max(auc)]])#, main=load.occ$spec_select,xlab="Purcentage(%)")
        })
      })


      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("GLM"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      
      out<-list(out,
                sidebarPanel(
                  selectInput("choice_block_GLM", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                              c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                  ),
                  conditionalPanel(
                    condition = "input.choice_block_GLM == 'Modelling without spatial blocking'",
                    sliderInput("number_no_block_fold_GLM", "Please set the number of fold", min = 1, max = 100, value = 5)
                  )),
                
                mainPanel(width = 6, tabsetPanel(type = "tabs",
                                                 tabPanel("Specie predictors",
                                                          selectInput('var_expl_GLM', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                          myActionButton("GLM",label=("Apply GLM"), "primary"),
                                                          plotOutput("enfa_var")
                                                 ),
                                                 tabPanel("Map",
                                                          
                                                          selectInput('probaplot_GLM', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("proba_occ_GLM")
                                                          
                                                 ),
                                                 tabPanel("Model Evaluation",
                                                          selectInput('model_ev_GLM', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("eval_GLM")
                                                 ),
                                                 tabPanel("Variable response",
                                                          selectInput('response_var_GLM', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("response_eco_GLM")
                                                 ),
                                                 tabPanel("Variable Importance",
                                                          plotOutput("var_importance_GLM")
                                                 )
                                                 
                                                 
                ),
                id = "tabs")
      )
      out
  
    })
    ### end GLM
    
    ### GAM contents ###
    output$ui_GAM<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$MaxEnt,{
        validate(
          need(length(input$var_expl) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl)
        Specdata<-Specdata()
        set.seed(1994)
        fold<-dismo::kfold(Specdata,input$number_no_block_fold)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {
          
          # testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          # testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
          
          p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          #test<-Specdata[fold == i, ] 
          
          occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          
          bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          model[[i]] <- dismo::maxent(data$enfa, p, a) #, factors='Sol'
          evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
          
        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        # data$proba_occ<-model_pred[["espece"]]
        # data$pre_abs<-model_pred[["PresenceAbsence"]]
        # output$proba_occ<-renderPlot({
        #   ggR_P(model_pred[["espece"]])
        # })
        # output$pre_ab_map<-renderPlot({
        #   PASpecies(model_pred[["PresenceAbsence"]])
        # })
        # output$timesraster<-renderPlot({
        #   ggR_P(model_pred[["ProbaPresence"]])
        # })
        observeEvent(input$probaplot,{
          if(input$probaplot=='Probability of occurence(absence/presence)'){
            title_probaplot<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot=='Presence/Absence'){
            title_probaplot<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]
            
          }
          if(input$probaplot=='Probability of occurence(presence)'){
            title_probaplot<-'Probability of occurence(presence)'
            map<-model_pred[["ProbaPresence"]]
          }
          output$proba_occ<-renderPlot({
            if(title_probaplot=='Presence/Absence'){PASpecies(map)}
            else{
              ggR_P(map) 
            }
            
          })
        })
        observeEvent(input$model_ev,{
          if(input$model_ev == 'ROC') {ev<-'ROC'}
          if(input$model_ev == 'density') {ev<-'density'}
          if(input$model_ev == 'boxplot') {ev<-'boxplot'}
          if(input$model_ev == 'kappa') {ev<-'kappa'}
          if(input$model_ev == 'FPR') {ev<-'FPR'}
          if(input$model_ev == 'prevalence') {ev<-'prevalence'}
          output$eval<-renderPlot({
            if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev)
              }
            }
            
            
          })
        })
        observeEvent(input$response_var,{
          output$response_eco<-renderPlot({
            dismo::response(model[[which.max(auc)]],var=input$response_var,main=load.occ$spec_select)
          })
        })
        output$var_importance<-renderPlot({
          plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
        })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("Maxent"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      out <- list(out,
                  column(12,offset=0,sidebarPanel(
                    selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                                c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                    ),
                    conditionalPanel(
                      condition = "input.choice_block == 'Modelling without spatial blocking'",
                      sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
                    )),align="center"))
      out <- list(out,
                  fluidRow(
                    box(title = "Specie predictors",
                        selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE)),
                    box(title = "Ecological Niche Factor Analysis",
                        plotOutput("enfa_var"))))
      
      out <- list(out,
                  fluidRow(actionButton('MaxEnt', 'Apply')))
      out <- list(out,
                  fluidRow(
                    box(title='Probability of occurence',
                        selectInput('probaplot', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                        plotOutput("proba_occ")),
                    box(title = "Model Evaluation",
                        selectInput('model_ev', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                        plotOutput("eval")),
                    box(title = 'Variable response',
                        selectInput('response_var', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                        plotOutput("response_eco")),
                    box(title = 'Variable importance',
                        plotOutput("var_importance"))))
      out
      
    })
    ### end GAM
    ### BRT contents ##
    output$ui_BRT<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$MaxEnt,{
        validate(
          need(length(input$var_expl) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl)
        Specdata<-Specdata()
        set.seed(1994)
        fold<-dismo::kfold(Specdata,input$number_no_block_fold)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {
          
          # testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          # testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
          
          p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          #test<-Specdata[fold == i, ] 
          
          occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          
          bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          model[[i]] <- dismo::maxent(data$enfa, p, a) #, factors='Sol'
          evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
          
        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        # data$proba_occ<-model_pred[["espece"]]
        # data$pre_abs<-model_pred[["PresenceAbsence"]]
        # output$proba_occ<-renderPlot({
        #   ggR_P(model_pred[["espece"]])
        # })
        # output$pre_ab_map<-renderPlot({
        #   PASpecies(model_pred[["PresenceAbsence"]])
        # })
        # output$timesraster<-renderPlot({
        #   ggR_P(model_pred[["ProbaPresence"]])
        # })
        observeEvent(input$probaplot,{
          if(input$probaplot=='Probability of occurence(absence/presence)'){
            title_probaplot<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot=='Presence/Absence'){
            title_probaplot<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]
            
          }
          if(input$probaplot=='Probability of occurence(presence)'){
            title_probaplot<-'Probability of occurence(presence)'
            map<-model_pred[["ProbaPresence"]]
          }
          output$proba_occ<-renderPlot({
            if(title_probaplot=='Presence/Absence'){PASpecies(map)}
            else{
              ggR_P(map) 
            }
            
          })
        })
        observeEvent(input$model_ev,{
          if(input$model_ev == 'ROC') {ev<-'ROC'}
          if(input$model_ev == 'density') {ev<-'density'}
          if(input$model_ev == 'boxplot') {ev<-'boxplot'}
          if(input$model_ev == 'kappa') {ev<-'kappa'}
          if(input$model_ev == 'FPR') {ev<-'FPR'}
          if(input$model_ev == 'prevalence') {ev<-'prevalence'}
          output$eval<-renderPlot({
            if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev)
              }
            }
            
            
          })
        })
        observeEvent(input$response_var,{
          output$response_eco<-renderPlot({
            dismo::response(model[[which.max(auc)]],var=input$response_var,main=load.occ$spec_select)
          })
        })
        output$var_importance<-renderPlot({
          plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
        })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("Maxent"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      out <- list(out,
                  column(12,offset=0,sidebarPanel(
                    selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                                c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                    ),
                    conditionalPanel(
                      condition = "input.choice_block == 'Modelling without spatial blocking'",
                      sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
                    )),align="center"))
      out <- list(out,
                  fluidRow(
                    box(title = "Specie predictors",
                        selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE)),
                    box(title = "Ecological Niche Factor Analysis",
                        plotOutput("enfa_var"))))
      
      out <- list(out,
                  fluidRow(actionButton('MaxEnt', 'Apply')))
      out <- list(out,
                  fluidRow(
                    box(title='Probability of occurence',
                        selectInput('probaplot', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                        plotOutput("proba_occ")),
                    box(title = "Model Evaluation",
                        selectInput('model_ev', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                        plotOutput("eval")),
                    box(title = 'Variable response',
                        selectInput('response_var', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                        plotOutput("response_eco")),
                    box(title = 'Variable importance',
                        plotOutput("var_importance"))))
      out
      
    })
    ### end BRT
    
    ### RandomForest contents
    output$ui_RF<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$RandomForest,{
        validate(
          need(length(input$var_expl_RF) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl_RF)
        pa_data<-reactive({
          pa_data<-sf::st_as_sf(Specdata(), coords = c("lon","lat"), crs = crs(data$enfa))
          pa_data
        })
        # extract the raster values for the species points as a dataframe
        mydataF<-reactive({
          mydataF <- raster::extract(data$enfa, pa_data(), df = TRUE)
          mydataF <- mydataF[,-1]
          Specdata<-Specdata()
          spec<-Specdata[,load.occ$spec_select]
          mydataF<-cbind(mydataF,spec)
          mydataF<-mydataF %>% dplyr::rename(Species=spec)
          mydataF
        })
        mydataF<-mydataF()
        set.seed(1994)
        fold<-dismo::kfold(Specdata(),input$number_no_block_fold_RF)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {
          
          testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          
          model[[i]] <- randomForest(as.numeric(Species)~., mydataF[fold != i, ], na.action=na.omit,importance=TRUE)
          evaluate_model[[i]] <- dismo::evaluate(testpres, testbackg, model[[i]])
          
        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        observeEvent(input$probaplot_RF,{
          if(input$probaplot_RF=='Probability of occurence(absence/presence)'){
            title_probaplot_RF<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot_RF=='Presence/Absence'){
            title_probaplot_RF<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]
            
          }
          if(input$probaplot_RF=='Probability of occurence(presence)'){
            title_probaplot_RF<-'Probability of occurence(presence)'
            map<-model_pred[["ProbaPresence"]]
          }
          output$proba_occ_RF<-renderPlot({
            if(title_probaplot_RF=='Presence/Absence'){PASpecies(map)}
            else{
              ggR_P(map) 
            }
            
          })
        })
        observeEvent(input$model_ev_RF,{
          if(input$model_ev_RF == 'ROC') {ev_RF<-'ROC'}
          if(input$model_ev_RF == 'density') {ev_RF<-'density'}
          if(input$model_ev_RF == 'boxplot') {ev_RF<-'boxplot'}
          if(input$model_ev_RF == 'kappa') {ev_RF<-'kappa'}
          if(input$model_ev_RF == 'FPR') {ev_RF<-'FPR'}
          if(input$model_ev_RF == 'prevalence') {ev_RF<-'prevalence'}
          output$eval_RF<-renderPlot({
            if(ev_RF=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev_RF=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev_RF)
              }
            }
            
            
          })
        })
        observeEvent(input$response_var_RF,{
          output$response_eco_RF<-renderPlot({
            dismo::response(model[[which.max(auc)]],var=input$response_var_RF,main=load.occ$spec_select)
          })
        })
        output$var_importance_RF<-renderPlot({
          #plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
          randomForest::varImpPlot(model[[which.max(auc)]],main=load.occ$spec_select)
        })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("Random Forest"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      out<-list(out,
                sidebarPanel(
                  selectInput("choice_block_RF", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                              c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                  ),
                  conditionalPanel(
                    condition = "input.choice_block_RF == 'Modelling without spatial blocking'",
                    sliderInput("number_no_block_fold_RF", "Please set the number of fold", min = 1, max = 100, value = 5)
                  )),
                
                mainPanel(width = 6, tabsetPanel(type = "tabs",
                                                 tabPanel("Specie predictors",
                                                          selectInput('var_expl_RF', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                          myActionButton("RandomForest",label=("Apply Random Forest"), "primary"),
                                                          plotOutput("enfa_var")
                                                 ),
                                                 tabPanel("Map",
                                                          
                                                          selectInput('probaplot_RF', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("proba_occ_RF")
                                                          
                                                 ),
                                                 tabPanel("Model Evaluation",
                                                          selectInput('model_ev_RF', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("eval_RF")
                                                 ),
                                                 tabPanel("Variable response",
                                                          selectInput('response_var_RF', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("response_eco_RF")
                                                 ),
                                                 tabPanel("Variable Importance",
                                                          plotOutput("var_importance_RF")
                                                 )
                                                 
                                                 
                ),
                id = "tabs")
      )
      out
      
    })
    ### en RandomForest
    
    ### SVM contents
    output$ui_SVM<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$SVM,{
        validate(
          need(length(input$var_expl_SVM) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl_SVM)
        pa_data<-reactive({
          pa_data<-sf::st_as_sf(Specdata(), coords = c("lon","lat"), crs = crs(data$enfa))
          pa_data
        })
        # extract the raster values for the species points as a dataframe
        mydataF<-reactive({
          mydataF <- raster::extract(data$enfa, pa_data(), df = TRUE)
          mydataF <- mydataF[,-1]
          Specdata<-Specdata()
          spec<-Specdata[,load.occ$spec_select]
          mydataF<-cbind(mydataF,spec)
          mydataF<-mydataF %>% dplyr::rename(Species=spec)
          mydataF
        })
        mydataF<-mydataF()
        Specdata<-Specdata()
        set.seed(1994)
        fold<-dismo::kfold(Specdata,input$number_no_block_fold_SVM)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {
          
          testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
          
          testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          model[[i]] <- ksvm(as.numeric(Species)~., mydataF[fold != i, ], na.action=na.omit)
          evaluate_model[[i]] <- dismo::evaluate(testpres, testbackg, model[[i]])
          
          
        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        observeEvent(input$probaplot_SVM,{
          if(input$probaplot_SVM=='Probability of occurence(absence/presence)'){
            title_probaplot_SVM<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot_SVM=='Presence/Absence'){
            title_probaplot_SVM<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]
            
          }
          if(input$probaplot_SVM=='Probability of occurence(presence)'){
            title_probaplot_SVM<-'Probability of occurence(presence)'
            map<-model_pred[["ProbaPresence"]]
          }
          output$proba_occ_SVM<-renderPlot({
            if(title_probaplot_SVM=='Presence/Absence'){PASpecies(map)}
            else{
              ggR_P(map) 
            }
            
          })
        })
        observeEvent(input$model_ev_SVM,{
          if(input$model_ev_SVM == 'ROC') {ev_SVM<-'ROC'}
          if(input$model_ev_SVM == 'density') {ev_SVM<-'density'}
          if(input$model_ev_SVM == 'boxplot') {ev_SVM<-'boxplot'}
          if(input$model_ev_SVM == 'kappa') {ev_SVM<-'kappa'}
          if(input$model_ev_SVM == 'FPR') {ev_SVM<-'FPR'}
          if(input$model_ev_SVM == 'prevalence') {ev_SVM<-'prevalence'}
          output$eval_SVM<-renderPlot({
            if(ev_SVM=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev_SVM=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev_SVM)
              }
            }
            
            
          })
        })
        observeEvent(input$response_var_SVM,{
          output$response_eco_SVM<-renderPlot({
            dismo::response(model[[which.max(auc)]],var=input$response_var_SVM,main=load.occ$spec_select)
          })
        })
        output$var_importance<-renderPlot({
          plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
        })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("SVM"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      out<-list(out,
                sidebarPanel(
                  selectInput("choice_block_SVM", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                              c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                  ),
                  conditionalPanel(
                    condition = "input.choice_block_SVM == 'Modelling without spatial blocking'",
                    sliderInput("number_no_block_fold_SVM", "Please set the number of fold", min = 1, max = 100, value = 5)
                  )),
                
                mainPanel(width = 6, tabsetPanel(type = "tabs",
                                                 tabPanel("Specie predictors",
                                                          selectInput('var_expl_SVM', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE),
                                                          myActionButton("SVM",label=("Apply SVM"), "primary"),
                                                          plotOutput("enfa_var")
                                                 ),
                                                 tabPanel("Map",
                                                          
                                                          selectInput('probaplot_SVM', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("proba_occ_SVM")
                                                          
                                                 ),
                                                 tabPanel("Model Evaluation",
                                                          selectInput('model_ev_SVM', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("eval_SVM")
                                                 ),
                                                 tabPanel("Variable response",
                                                          selectInput('response_var_SVM', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                                                          plotOutput("response_eco_SVM")
                                                 ),
                                                 tabPanel("Variable Importance",
                                                          plotOutput("var_importance_SVM")
                                                 )
                                                 
                                                 
                ),
                id = "tabs")
      )
      out
      
    })
    ### end SVM contents
    
    ### Combining model contents
    output$ui_Combining<-renderUI({
      
      Specdata<-reactive({
        dsf<-load.occ$select
        dsf<-dsf %>% dplyr::rename(lon=load.occ$lon,lat=load.occ$lat)
        dsf
      })
      
      Specdata_Presence<-reactive({
        dsf<-Specdata()
        dsf<-dsf[dsf[,ncol(dsf)] == 1,]
        sp::coordinates(dsf) <-~lon+lat
        sp::proj4string(dsf) <-raster::crs(data$Env)
        dsf
      })
      
      glc<-reactive({
        GLcenfa(x = data$Env)
      })
      
      mod.enfa<-reactive({
        pr<-Specdata_Presence()
        pr@data$load.occ$spec_select<-as.numeric(pr@data$load.occ$spec_select)
        CENFA::enfa(x = data$Env, s.dat = pr, field = load.occ$spec_select)
      })
      enfa_plot<-reactive({
        glc <- glc()
        
        mod.enfa <- mod.enfa()
        CENFA::scatter(x = mod.enfa, y = glc,n=nlayers(data$Env),p=1)
      })
      output$enfa_var<-renderPlot({
        enfa_plot()
      })
      
      observeEvent(input$MaxEnt,{
        validate(
          need(length(input$var_expl) > 0, 'Choose specie predictors first !')
        )
        
        data$enfa<-raster::subset(data$Env,input$var_expl)
        Specdata<-Specdata()
        set.seed(1994)
        fold<-dismo::kfold(Specdata,input$number_no_block_fold)
        #fold<-kfold()
        model<-list()
        evaluate_model<-list()
        for (i in 1:5) {
          
          # testpres <- mydataF[mydataF[fold == i,ncol(mydataF)] == 1, 1:(ncol(mydataF)-1)]
          # testbackg <- mydataF[mydataF[fold == i,ncol(mydataF)] == 0, 1:(ncol(mydataF)-1)]
          #dsf<-dsf[dsf[,ncol(dsf)] == 1,]
          
          p<-Specdata[Specdata[fold != i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          a<-Specdata[Specdata[fold != i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          #test<-Specdata[fold == i, ] 
          
          occtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 1, 1:(ncol(Specdata)-1)]
          
          bgtest<-Specdata[Specdata[fold == i,ncol(Specdata)] == 0, 1:(ncol(Specdata)-1)]
          model[[i]] <- dismo::maxent(data$enfa, p, a) #, factors='Sol'
          evaluate_model[[i]] <- evaluate(occtest, bgtest, model[[i]], data$enfa)
          
        }
        model_pred<-list()
        auc <- sapply(evaluate_model, function(x){x@auc})
        model_pred[["espece"]]<-predict(data$enfa, model[[which.max(auc)]])
        model_pred[["AUC"]]<-auc[which.max(auc)]
        model_pred[["threshold"]]<- threshold(evaluate_model[[which.max(auc)]], 'spec_sens')
        model_pred[["PresenceAbsence"]]<-model_pred[["espece"]]>model_pred[["threshold"]]
        model_pred[["ProbaPresence"]]<-TimesRasters(model_pred[["espece"]],model_pred[["PresenceAbsence"]])
        # data$proba_occ<-model_pred[["espece"]]
        # data$pre_abs<-model_pred[["PresenceAbsence"]]
        # output$proba_occ<-renderPlot({
        #   ggR_P(model_pred[["espece"]])
        # })
        # output$pre_ab_map<-renderPlot({
        #   PASpecies(model_pred[["PresenceAbsence"]])
        # })
        # output$timesraster<-renderPlot({
        #   ggR_P(model_pred[["ProbaPresence"]])
        # })
        observeEvent(input$probaplot,{
          if(input$probaplot=='Probability of occurence(absence/presence)'){
            title_probaplot<-'Probability of occurence(absence/presence)'
            map<-model_pred[["espece"]]}
          if(input$probaplot=='Presence/Absence'){
            title_probaplot<-'Presence/Absence'
            map<-model_pred[["PresenceAbsence"]]
            
          }
          if(input$probaplot=='Probability of occurence(presence)'){
            title_probaplot<-'Probability of occurence(presence)'
            map<-model_pred[["ProbaPresence"]]
          }
          output$proba_occ<-renderPlot({
            if(title_probaplot=='Presence/Absence'){PASpecies(map)}
            else{
              ggR_P(map) 
            }
            
          })
        })
        observeEvent(input$model_ev,{
          if(input$model_ev == 'ROC') {ev<-'ROC'}
          if(input$model_ev == 'density') {ev<-'density'}
          if(input$model_ev == 'boxplot') {ev<-'boxplot'}
          if(input$model_ev == 'kappa') {ev<-'kappa'}
          if(input$model_ev == 'FPR') {ev<-'FPR'}
          if(input$model_ev == 'prevalence') {ev<-'prevalence'}
          output$eval<-renderPlot({
            if(ev=='density'){density(evaluate_model[[which.max(auc)]])}
            else{
              if(ev=='boxplot'){boxplot(evaluate_model[[which.max(auc)]], col=c('red', 'green'),xlab=load.occ$spec_select)}
              else{
                plot(evaluate_model[[which.max(auc)]],ev)
              }
            }
            
            
          })
        })
        observeEvent(input$response_var,{
          output$response_eco<-renderPlot({
            dismo::response(model[[which.max(auc)]],var=input$response_var,main=load.occ$spec_select)
          })
        })
        output$var_importance<-renderPlot({
          plot(model[[which.max(auc)]], main=load.occ$spec_select,xlab="Purcentage(%)")
        })
      })
      
      
      
      
      out <- NULL
      txt_setup<-'The Maxent software is based on the maximum-entropy approach for modeling species niches and distributions. From a set of environmental (e.g., climatic) grids and georeferenced occurrence localities (e.g. mediated by GBIF), the model expresses a probability distribution where each grid cell has a predicted suitability of conditions for the species. Maxent is a stand-alone Java application and can be used on any computer running Java version 1.5 or later.'
      out <- fluidRow(
        column(width = 12, offset = 0, h3("Maxent"), class="wb-header"),
        column(width = 12, offset = 0, p("The first step is to choose specie predictors accordint to ENFA or other source, afther apply Maxent method."), class="wb-header-hint"),
        fluidRow(column(12, h4("Read Me", tipify(icon("info-circle"), title=txt_setup, placement="bottom"), class="wb-block-title"), align="center"))
      )
      out <- list(out,
                  column(12,offset=0,sidebarPanel(
                    selectInput("choice_block", "Please Choose your model technic (without spatial blocking or with spatial blocking)",
                                c(without="Modelling without spatial blocking",with="Modelling with spatial blocking")
                    ),
                    conditionalPanel(
                      condition = "input.choice_block == 'Modelling without spatial blocking'",
                      sliderInput("number_no_block_fold", "Please set the number of fold", min = 1, max = 100, value = 5)
                    )),align="center"))
      out <- list(out,
                  fluidRow(
                    box(title = "Specie predictors",
                        selectInput('var_expl', 'Please select the specie predictors', names(data$Env), multiple = TRUE, selectize = TRUE)),
                    box(title = "Ecological Niche Factor Analysis",
                        plotOutput("enfa_var"))))
      
      out <- list(out,
                  fluidRow(actionButton('MaxEnt', 'Apply')))
      out <- list(out,
                  fluidRow(
                    box(title='Probability of occurence',
                        selectInput('probaplot', '', c("Probability of occurence(absence/presence)","Presence/Absence","Probability of occurence(presence)"), multiple = FALSE, selectize = TRUE),
                        plotOutput("proba_occ")),
                    box(title = "Model Evaluation",
                        selectInput('model_ev', 'Please select the metric to evaluate the model', c("ROC","density","boxplot","kappa","FPR","prevalence"), multiple = FALSE, selectize = TRUE),
                        plotOutput("eval")),
                    box(title = 'Variable response',
                        selectInput('response_var', 'Please select the variable to get its ecological response', names(data$enfa), multiple = FALSE, selectize = TRUE),
                        plotOutput("response_eco")),
                    box(title = 'Variable importance',
                        plotOutput("var_importance"))))
      out
      
    })
    ### end Combining model contents
    output$ui_Models_main <- renderUI({
      out <- NULL
      val <- obj$cur_selection_results
      ## Categorical (defined in controller/ui_results_imputation.R)
      if (val=="btn_Models_results_1") {
        return(uiOutput("ui_bioclim"))
      }
      if (val=="btn_Models_results_2") {
        return(uiOutput("ui_domain"))
      }
      if (val=="btn_Models_results_3") {
        return( uiOutput("ui_mahal"))
      }
      if (val=="btn_Models_results_4") {
        return(uiOutput("ui_GLM"))
      }
      if (val=="btn_Models_results_5") {
        return(uiOutput("ui_GAM"))
      }
      if (val=="btn_Models_results_6") {
        return( uiOutput("ui_MaxEnt"))
      }
      if (val=="btn_Models_results_7") {
        return(uiOutput("ui_BRT"))
      }
      if (val=="btn_Models_results_8") {
        return(uiOutput("ui_RF"))
      }
      if (val=="btn_Models_results_9") {
        return( uiOutput("ui_SVM"))
      }
      if (val=="btn_Models_results_10") {
        return( uiOutput("ui_Combining"))
      }
    })	
    output$ui_Models_sidebar_left <- renderUI({
      output$ui_sel_Models_btns <- renderUI({
        cc1 <- c("Bioclim","Domain","Mahalanobis distance")
        cc2 <- c("Generalized Linear Models", "Generalized Additive Models")
        cc3 <- c("MaxEnt", "Boosted Regression Trees", "Random Forest","Support Vector Machines")
        cc4 <- c("Combining model")
        df <- data.frame(lab=c(cc1,cc2,cc3,cc4), header=NA)
        df$header[1] <- "Profile methods"
        df$header[4] <- "Classical regression models"
        df$header[6] <- "Machine learning methods"
        df$header[10] <- "Combining model predictions"
        out <- NULL
        for (i in 1:nrow(df)) {
          id <- paste0("btn_Models_results_",i)
          if (obj$cur_selection_results==id) {
            style <- "primary"
          } else {
            style <- "default"
          }
          if (!is.na(df$header[i])) {
            out <- list(out, fluidRow(column(12, h4(df$header[i]), align="center")))
          }
          out <- list(out, fluidRow(
            column(12, bsButton(id, label=df$lab[i], block=TRUE, size="extra-small", style=style))
          ))
        }
        out
      })
      # required observers that update the color of the active button!
      eval(parse(text=genObserver_menus(pat="btn_Models_results_", n=1:10, updateVal="cur_selection_results")))
      return(uiOutput("ui_sel_Models_btns"))
    })
    output$ui_Models_noproblem <- renderUI({
      return(list(
        noInputData(uri="ui_Models"),
        fluidRow(column(12, tags$br(), p(""), align="center"))
        #fluidRow(column(12, myActionButton("nodata_anonymize_uploadproblem", label="Upload a previously saved problem", btn.style="primary"), align="center"))
      ))
    })
    output$ui_Models <- renderUI({
      if(length(input$Occ)==0){
        return(uiOutput("ui_Models_noproblem"))}
      {
      fluidRow(
        #column(12,offset=0,radioButtons('choice_block', 'Please Choose your model technic', choices=c(No_sp_blocking="Modelling without spatial blocking",sp_blocking="Modelling with spatial blocking"), selected="Modelling without spatial blocking", inline=TRUE),align="center", class="wb-header" ),
        column(2, uiOutput("ui_Models_sidebar_left"), class="wb_sidebar"),
        column(10, uiOutput("ui_Models_main"), class="wb-maincolumn")
        )
      }
    }
    )
    ######################################## end biomod2 contents ##########################################
    ########################################################################
  }
  shinyApp(ui, server)
  }
cirad_hema()