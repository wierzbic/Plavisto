# Plavisto - Plastid Genome Visualization Tool
# Version 1.0 (01.19.2022)
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

library(shiny)

shinyUI(fluidPage( 
  title = "Plavisto",
  h2("Plavisto"),
  h4(a("Wierzbicki Lab", href="https://sites.lsa.umich.edu/wierzbicki-lab/", target="_blank"), HTML('<u>Pl</u>astid Genome <u>Vis</u>ualization <u>To</u>ol')),
  tags$head(includeHTML(("analytics.html"))),
  tabsetPanel(     
      tabPanel("Browser",
          fluidRow(  
              column(12,
                  lapply(1:number_of_plots, function(i){                                                                 #Generate the desired number of plot panels
                    conditionalPanel(
                      condition = paste0('input.plot', i, '_type != 3'),                                                 #Only show selected plots
                      plotOutput(paste0('Plot', i), height=plot_height, click = "big_plot_click", brush = "big_plot_brush")      #Define plot height and events
                    )  
                  }),   
                  plotOutput('naviPlot', width="auto", height=40, click = "plot_click")                                  #Define navigation plot
                  )
              ), 
          fluidRow( 
              column(4,
                  h4("Shown Region Size"),
                  actionButton("Zoom100", "100 bp"),
                  actionButton("Zoom1000", "1 kb"),
                  actionButton("Zoom10000", "10 kb"),
                  actionButton("Zoom50000", "50 kb"),
                  actionButton("ZoomAll", "All"),
                  numericInput("size", NULL, 5000, 100, genome_length, 100),
                  h4("Number of bins"),
                  numericInput("bin_size", NULL, 500, 10, 2000, 100),
                  h4("Zoom"),
                  numericInput("zoom", NULL, 1, 1, 100, 1)
                  ),
              column(3,
                  fluidRow(
                    h4("Navigation"),
                    actionButton("ScrollL", "Scroll Left"),
                    actionButton("ScrollR", "Scroll Right"),
                    h4(paste0('Select type of plot ', 1)),
                    radioButtons(paste0('plot', 1, '_type'), NULL, choices = c("Show plot" = 2, "No plot" = 3), inline = 1, selected = 2),
                  ),
                  lapply(2:number_of_plots, function(i){ 
                    fluidRow(
                      h4(paste0('Select type of plot ', i)),
                      radioButtons(paste0('plot', i, '_type'), NULL, choices = c("Show plot" = 2, "No plot" = 3), inline = 1, selected = 3)
                    )
                  })
                  ),
              column(5,
                  lapply(1:number_of_plots, function(i){                                                                 #Generate the desired number of radio buttons
                    fluidRow(
                      #radioButtons(paste0('plot', i, '_type'), paste0('Select type of plot ', i), choices = c("RNA plot" = 1, "DNA plot" = 2, "None" = 3), inline = 1, selected = 2),
                      conditionalPanel(                                                                                    #display the rest only if a plot is selected
                        condition = paste0('input.plot', i, '_type != 3'),
                        h4(paste0("Options for plot ", i)),
                        radioButtons(paste0('plot', i, '_scale'), label = NULL, choices = c("Linear" = 0, "Log" = 1), inline = 1, selected = 0),
                        radioButtons(paste0('plot', i, '_diff'), label = NULL, choices = c("Absolute values" = 0, "Differential" = 1), inline = 1, selected = 0),
                        radioButtons(paste0('plot', i, '_norm'), label = NULL, choices = c("Global scaling" = 0, "Local scaling" = 1, "No scaling" = 3), inline = 1, selected = 3)
                      )
                    )
                  }) 
                  
              )    
          )   
      ),  
      tabPanel("Datasets", 
          lapply(1:number_of_plots, function(i){     
            column(5, 
                   h4(paste0("Plot ", i)),
                   jstreeOutput(paste0('tree', i), height = "100%", width = "100%")  #Generate the desired number of trees in the Datasets panel
            )
          })
      ), 
      tabPanel("Instructions", 
          includeHTML("www/info3.html"),
          uiOutput("users_number")
      ) 
    ) 
))    
