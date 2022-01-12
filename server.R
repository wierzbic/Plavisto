# Plavisto - Plastid Genome Visualization Tool
# Version 1.0
# Copyright (C) 2022  Andrzej Wierzbicki
#   
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

rm(list=ls())
source("pgvt_functions_v1_0.R", local = TRUE)

shinyServer(function(input, output, session) {

#define variables
position <- reactiveValues(data = NULL)
position$start <- 54000
position$end <- 59000
position$size <- 5000
outputs <- reactiveValues(data = NULL)
outputs$position <- 0
outputs$yvalue <- 0
click_equalize_button <- 0

#general constants
plot_cutoff <- 0
coverage_plot_lim <- data.frame()
tmp <- data.frame()
#coverage_plot_globallist <- list()

coverage_plot_globallist <- vector(mode = "list", length = number_of_plots)
wc_globallist <- list()

#Load annotations datasets
annotation <- read.csv('plastid_genes.csv', header=FALSE) 
annotation$V8[annotation$V3 == '+'] <- 1
annotation$V8[annotation$V3 == '-'] <- -1
operons <- read.csv('plastid_operons.csv', header=FALSE) 
operons$V7[operons$V3 == '+'] <- 1
operons$V7[operons$V3 == '-'] <- -1

#showModal(modalDialog("Please wait while datasets are loading", footer=NULL))

#Load initial dataset
data_to_load <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE, na.strings=c("NA", ""), colClasses = c("character", "character", "character", "character", "character", "integer"))
for(i in 1:number_of_plots){    
  coverage_plot <- data.frame()
  #browser()
  default_datasets <- data_to_load[data_to_load$default == i ,]
  if(nrow(default_datasets) > 0){                                         #Test if any default datasets are defined for this lane
    for(dd in 1:nrow(default_datasets)) {                                 #Go through all rows representing default datasets; use for loop to return dataframe ready for ggplot2
      if(!is.na(default_datasets$filename_plus[dd])){                                                                     #if filename for plus strand exists, load the file
        if(nrow(coverage_plot) == 0) coverage_plot <- load_data(default_datasets$filename_plus[dd], "+", default_datasets$samplename[dd])     #If first dataset just load
          else coverage_plot <- rbind(coverage_plot, load_data(default_datasets$filename_plus[dd], "+", default_datasets$samplename[dd]))     #If some data exist, add
        cat(paste("Loaded dataset: ", default_datasets$samplename[dd], " file name: ",default_datasets$filename_plus[dd], "\n"))              #Print message
      }
      if(!is.na(default_datasets$filename_minus[dd])){                                                                     #if filename for minus strand exists, load the file
        if(nrow(coverage_plot) == 0) coverage_plot <- load_data(default_datasets$filename_minus[dd], "-", default_datasets$samplename[dd])
        else coverage_plot <- rbind(coverage_plot, load_data(default_datasets$filename_minus[dd], "-", default_datasets$samplename[dd]))
        cat(paste("Loaded dataset: ", default_datasets$samplename[dd], " file name: ",default_datasets$filename_minus[dd], "\n"))
      }
    }
    coverage_plot_tmp <- group_by(coverage_plot, name)
    wc <- summarize(coverage_plot_tmp, count = sum(data))
    coverage_plot_globallist[[i]] <- coverage_plot
    wc_globallist[[i]] <- wc
  }
}
  
#Scrolling and zooming of the plots
  observeEvent(input$ScrollL, {                     #scroll left
    position$start <- position$start - input$size
    position$end <- position$end - input$size
    correctPosition()
  })
  observeEvent(input$ScrollR, {                     #scroll right
    position$start <- position$start + input$size
    position$end <- position$end + input$size
    correctPosition()
  })
  observeEvent(input$size, {                        #zoom in or out
    req(input$size)
    position$start <- position$end - (position$end - position$start)/2 - (input$size / 2)
    position$end <- position$start + input$size
    #browser()
    if(input$size == 0) {position$end <- position$start + 100}
    correctPosition()
  })
  observeEvent(input$Zoom100,   { updateSliderInput(session, "size", value = 100) })
  observeEvent(input$Zoom1000,  { updateSliderInput(session, "size", value = 1000)  })  
  observeEvent(input$Zoom10000, { updateSliderInput(session, "size", value = 10000)  })   #define actions of zoom buttons
  observeEvent(input$Zoom50000, { updateSliderInput(session, "size", value = 50000)  })
  observeEvent(input$ZoomAll,   { updateSliderInput(session, "size", value = genome_length)  })
  correctPosition <- function()   #function to assure that display is within the covered genome sequence
  {
    if(position$start <= 0){
      position$start <- 1
      position$end <- 1 + input$size
    }
    if(position$end > genome_length){
      position$start <- genome_length - input$size
      position$end <- genome_length 
    }
  }
  observeEvent(input$plot_click, {        #scroll by clicking on plot
      position$start <- input$plot_click$x - (input$size / 2)
      position$end <- input$plot_click$x + (input$size / 2)
      if(position$start <= 0){
        position$start <- 1
        position$end <- 1 + input$size
      }
      if(position$end > genome_length){
        position$start <- genome_length - input$size
        position$end <- genome_length 
      }
  })

#Scroll based on click location
  observeEvent(input$big_plot_click, {
      position$start <- input$big_plot_click$x - (input$size / 2)
      position$end <- input$big_plot_click$x + (input$size / 2)
      if(position$start <= 0){
        position$start <- 1
        position$end <- 1 + input$size
      }
      if(position$end > genome_length){
        position$start <- genome_length - input$size
        position$end <- genome_length 
      }
  }) 

#Dataset selection using trees  
lapply(1:number_of_plots, function(i){  
  output[[paste0('tree', i)]] <- jsTree::renderJsTree({         #load tree created with jsTree with data from infiles.csv
    tree_data <- read.table('infiles.csv', sep=",", header=TRUE, stringsAsFactors = FALSE)
    nested_data <- apply(tree_data[,1:3],1,paste,collapse='/')  #generate data for the tree, separated with /
    nodestate1 <- tree_data$default == i                        #list check boxes to be checked on startup
    jsTree(nested_data, nodestate = nodestate1)
  })  
  })  
  
#React to changes in dataset selection on trees
lapply(1:number_of_plots, function(i){  
  observeEvent(input[[paste0('tree', i, '_update')]], {            #load the list os selected datasets from the tree created with jsTree
    #browser()
    current_selection <- input[[paste0('tree', i, '_update')]]$.current_tree
    selections <- as.data.frame(jsonlite::fromJSON(current_selection))
    if(nrow(selections) > 0){             #test if anything is checked; if all boxes are unchecked do nothing
      colnames(selections) = "col1"
      selections <- separate(selections, col1, c("groupname", "foldername", "samplename"), sep = "/", extra = "drop", fill = "right")
      selections <- drop_na(selections)   #remove rows with NAs - contain partial tree info
      dataset_list <- reselect_datasets_tree(i, selections$samplename, coverage_plot_globallist, wc_globallist)
      coverage_plot_globallist <<- dataset_list[[1]]                                                              #Watch for scoping issues!
      wc_globallist <<- dataset_list[[2]]                                                                         #Should update variables in server but not globally
      generate_one_plot(i, coverage_plot_globallist, wc_globallist)
    }
  })
})
  
#respond to changing plot type
lapply(1:number_of_plots, function(i){                                #Generate code for all plots
    observeEvent(input[[paste0('plot', i, '_type')]], {  
      #browser()
      #cat(ls(envir=.GlobalEnv))
      generate_one_plot(i, coverage_plot_globallist, wc_globallist)
    })
})

#Generate a main plot  
generate_one_plot <- function(plot_number, coverage_plot_globallist, wc_globallist) {
     if(input[[paste0('plot', plot_number, '_type')]] == 1) output[[paste0('Plot', plot_number)]] <- renderPlot({   #RNA plot
       req(input$bin_size)
       #plot_chl(plot_number, position$start, position$end, input$bin_size, input$zoom, input[[paste0('plot', plot_number, '_scale')]], input[[paste0('plot', plot_number, '_norm')]], click_equalize_button, outputs$position, outputs$yvalue, input[[paste0('plot', plot_number, '_diff')]], coverage_plot_globallist, wc_globallist, annotation, operons)
     })
     if(input[[paste0('plot', plot_number, '_type')]] == 2) output[[paste0('Plot', plot_number)]] <- renderPlot({   #DNA plot
        req(input$bin_size)
        plot_MNase(plot_number, position$start, position$end, input$bin_size, input$zoom, input[[paste0('plot', plot_number, '_scale')]], input[[paste0('plot', plot_number, '_norm')]], click_equalize_button, outputs$position, outputs$yvalue, input[[paste0('plot', plot_number, '_diff')]], coverage_plot_globallist, wc_globallist, annotation, operons)
     })
}

#Generate scatterplot
#output$scatterPlot <- renderPlotly({
#  plot_scatter(1, position$start, position$end, input$bin_size, input$zoom, input$plot1_scale, input$plot1_norm, click_equalize_button, outputs$position, outputs$yvalue, input$plot1_diff, coverage_plot_globallist, wc_globallist)
#})

#Generate navigation plot
output$naviPlot <- renderPlot({             #navigation plot
  plot_navi(position$start, position$end, annotation, operons)
})
#removeModal()
}) 