library(shinydashboard)
#library(shinyIncubator)
library(shiny)
library(plotly)
library(d3heatmap)
library(shinyjs)

dashboardPage(
  dashboardHeader(title = "STAR summary",titleWidth = 300),
  dashboardSidebar(width = 300,
                   div(style="overflow-y: scroll"),
                   tags$head(tags$style(HTML(".sidebar { height: 90vh; overflow-y: auto; }" ))),
                   sidebarMenu(
                     menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))),
                   uiOutput("projects"),
                   sidebarMenu(
                     menuItemOutput("menuitem_loaddata"),
                     menuItem("PhenoData", tabName = "pheno", icon = icon("hand-o-right")),
                     menuItem("FastQC Report", tabName = "fastqc", icon = icon("hand-o-right")),
                     menuItem("Bargraph-Samples", tabName = "bargraph", icon = icon("bar-chart")),
                     menuItem("Library Complexity Summary", tabName = "libcomp", icon = icon("hand-o-right")),
                     menuItem("Metrics", tabName = "Metrics", icon = icon("hand-o-right")),
                     menuItem("Mark Duplicates", tabName = "markdup", icon = icon("hand-o-right"))
                   )#end of sidebar menu

),

  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    useShinyjs(),
     tabItems(
       tabItem(tabName = "dashboard",
               box(
                 width = 12, status = "primary",solidHeader = TRUE,
                 title = "Mapped Summary",
                 uiOutput("plotUI")
               )
       ),
       tabItem(tabName = "pheno",DT::dataTableOutput('anno')),
       tabItem(tabName = "fastqc",htmlOutput("fastqc",height=800,width=1100)),
       tabItem(tabName = "bargraph",uiOutput("indoptions"),uiOutput("allsamp"),plotlyOutput("barplotsind_out",width = 1100, height = 800)),
  tabItem(tabName = "libcomp",
          DT::dataTableOutput('libcomplex'),hr(),
          fluidRow(
            column(6,uiOutput("xoptions")),
            column(6,uiOutput("yoptions"))),
          plotlyOutput("libc_bplot",width = 1200, height = 500)
  ),
  tabItem(tabName = "Metrics",
          DT::dataTableOutput('metrics'),
          fluidRow(
            column(6,uiOutput("mxoptions")),
            column(6,uiOutput("myoptions"))),
          plotlyOutput("metr_bplot",width = 1200, height = 500)
          
  ),
  tabItem(tabName = "markdup",
          DT::dataTableOutput('mrkdup'),
          fluidRow(
            column(6,uiOutput("dxoptions")),
            column(6,uiOutput("dyoptions"))),
          plotlyOutput("mrkdup_bplot",width = 1200, height = 500)
  ))))
