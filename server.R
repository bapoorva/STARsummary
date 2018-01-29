#library(shinyIncubator)
library(shiny)
library(shinyjs)
library(tidyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(NMF)
library(Biobase)
library(reshape2)
library(d3heatmap)
library(plotly)
library(shinyjs)
library(htmlwidgets)
library(DT)
library(FactoMineR)
library(factoextra)
library(shinyRGL)
library(rgl)
library(ReactomePA)
library(limma)
library(ggrepel)
library(dplyr)
shinyServer(function(input, output,session) {
  
  #Read the parameter file
  readexcel = reactive({
    file = read.csv("data/param.csv")
  })
  
  #Get Project list and populate drop-down
  output$projects = renderUI({
    excel=readexcel()
    prj=excel$projects
    selectInput("projects","Select a project",as.list(sort(as.character(prj))))
  })
  
  output$datasetTable<- renderTable({
    read.csv('data/param.csv',stringsAsFactors = F)
  }, digits = 1)
  
  
  #load Rdata corresponding to selected project
  fileload <- reactive({
    inFile = paste('data/',as.character(input$projects),'.RData',sep = '')
    load(inFile)
    loaddata=STAR
    return(loaddata)
  })
  
  #Read the star_summary file
  inputfile = reactive({
    STAR=fileload()
    inFile = STAR$starsummary
#     file <- inFile %>% mutate(library= gsub('.*STAR/','',library))  %>%
#       separate(library,c('id','run','lane','barcode'),sep='_') %>% mutate(pool=paste(run,lane,sep='_')) %>%
#       dplyr::select(-Startedjobon:-MappingspeedMillionofreadsperhour)
    file <- inFile %>% mutate(id= gsub('.*STAR/','',library)) %>%
            dplyr::select(-Startedjobon:-MappingspeedMillionofreadsperhour,-library)
          file=file %>% dplyr::select(id,Numberofinputreads:ofchimericreads)
  })
  
  output$menuitem_loaddata <- renderMenu({
    menuItem("Summary",  tagList(checkboxGroupInput("stats",label="Mapped Statistics",choices=list("Unique Reads"='unique',"Multi-Mapping Reads"='multi',"Unmapped Reads"='unmapped'),selected = 'unique')
    ),
    icon = icon("database"),tabName='mapsum')
    
  })
  
  #cleanupthe results file
  readfile = reactive({
    file=inputfile()
    #colnames(file)=c("id","run","lane","barcode","Numberofinputreads","Averageinputreadlength","Uniquelymappedreads_number","Uniquelymappedreads_percentage","Averagemappedlength","Numberofsplices_Total","Numberofsplices_Annotated_sjdb","Numberofsplices_GT_AG","Numberofsplices_GC_AG","Numberofsplices_AT_AC","Numberofsplices:Non-canonical","Mismatchrateperbase","Deletionrateperbase","Deletionaveragelength","Insertionrateperbase","Insertionaveragelength","Numberofreadsmappedtomultipleloci","percentageofreadsmappedtomultipleloci","Numberofreadsmappedtotoomanyloci","percentageofreadsmappedtotoomanyloci","percentageofreadsunmapped_toomanymismatches","percentageofreadsunmapped_tooshort","percentageofreadsunmapped_other","Numberofchimericreads","ofchimericreads","pool")
    colnames(file)=c("id","Numberofinputreads","Averageinputreadlength","Uniquelymappedreads_number","Uniquelymappedreads_percentage","Averagemappedlength","Numberofsplices_Total","Numberofsplices_Annotated_sjdb","Numberofsplices_GT_AG","Numberofsplices_GC_AG","Numberofsplices_AT_AC","Numberofsplices:Non-canonical","Mismatchrateperbase","Deletionrateperbase","Deletionaveragelength","Insertionrateperbase","Insertionaveragelength","Numberofreadsmappedtomultipleloci","percentageofreadsmappedtomultipleloci","Numberofreadsmappedtotoomanyloci","percentageofreadsmappedtotoomanyloci","percentageofreadsunmapped_toomanymismatches","percentageofreadsunmapped_tooshort","percentageofreadsunmapped_other","Numberofchimericreads","ofchimericreads")
    return(file)
  })
  
  #read phenofile and merge it with the summary files
  readpheno = reactive({
    inFile = paste('data/',as.character(input$projects),'/Phenodata.csv',sep = '')
    pheno=read.csv(inFile)
    return(pheno)
  })
  
  output$pheno = DT::renderDataTable({
    DT::datatable(readpheno(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "Library complex")
  })
  


  addResourcePath("library", "/srv/shiny-server")
  #addResourcePath("library", "/Users/bapoorva/Desktop/Shiny/rna-web_git")

  output$fastqc = renderUI({
    tags$iframe(seamless="seamless",src=paste0("library/STARsummary/data/www/",input$projects,"_multiqc_report.html"), height=1400, width=1300)
  })
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  #populate drop down (samples)
  output$indoptions <- renderUI({
    d<- inputfile()
    ind=sort(unique(d$id))
    selectInput("ind", "Select ID",as.list(as.character(ind)))
  })
  
  output$allsamp <- renderUI({
  checkboxInput("allsamp", label = "Click to view all samples", value = FALSE)
  })
  
  #create bar graph for selected samples
  barplotind = reactive({
    d<- inputfile()
    if(input$allsamp == TRUE){
      d=d
    }else{
    d=d[d$id==input$ind,]}
    p=d %>% dplyr::select(id,`Uniquelymappedreads`,`ofreadsmappedtomultipleloci`,`ofreadsmappedtotoomanyloci`,`ofreadsunmapped.toomanymismatches`,`ofreadsunmapped.tooshort`,`ofreadsunmapped.other`) %>%
      gather("maptype","perc",-id) 
    p=plot_ly(p,x=~id,y=~perc,color=~maptype,type='bar') %>%
      layout(title = "Bargraph - Samples",barmode='stack',xaxis = list(title = "ID"),yaxis = list(title = "percentage"),margin=list(b=120,pad=4))
    
   (gg <- ggplotly(p))
  })
  
 
  output$barplotsind_out = renderPlotly({
    barplotind()
  })
  
  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  #Create table for uniquely mapped reads with link to FASTQC html files
  table_unique = reactive({
    dt = readfile()
#     run=dt$run
#     lane=dt$lane
#     barcode=dt$barcode
#     link1_name=paste0(run,"_s_",lane,"_1_",barcode,"_fastqc.html")
#     link2_name=paste0(run,"_s_",lane,"_2_",barcode,"_fastqc.html")
#     dt$run=paste0("<a href='",link1_name,"'target='_blank'>",dt$run,"_1","</a>","\n","<a href='",link2_name,"'target='_blank'>",dt$run,"_2","</a>")
#     
    dt=as.data.frame(dt)
    dt=dt %>% dplyr::select(id:Insertionaveragelength)
    return(dt)
  })
  #Create table for multi-mapped reads with link to FASTQC html files
  table_multi = reactive({
    dt = readfile()
#     run=dt$run
#     lane=dt$lane
#     barcode=dt$barcode
#     link1_name=paste0("/fujfs/d3/MAGnet_RNAseq_v2/fastQC/",run,"_s_",lane,"_1_",barcode,"_fastqc.html")
#     link2_name=paste0("/fujfs/d3/MAGnet_RNAseq_v2/fastQC/",run,"_s_",lane,"_2_",barcode,"_fastqc.html")
#     dt$run=paste0("<a href='",link1_name,"'target='_blank'>",dt$run,"_1","</a>","\n","<a href='",link2_name,"'target='_blank'>",dt$run,"_2","</a>")
#     dt=as.data.frame(dt)
    dt=dt %>% dplyr::select(id,Numberofinputreads,Numberofreadsmappedtomultipleloci:percentageofreadsmappedtotoomanyloci)
    return(dt)
  })
  #Create table for unmapped reads with link to FASTQC html files
  table_unmapped = reactive({
     dt = readfile()
#     run=dt$run
#     lane=dt$lane
#     barcode=dt$barcode
#     link1_name=paste0("/fujfs/d3/MAGnet_RNAseq_v2/fastQC/",run,"_s_",lane,"_1_",barcode,"_fastqc.html")
#     link2_name=paste0("/fujfs/d3/MAGnet_RNAseq_v2/fastQC/",run,"_s_",lane,"_2_",barcode,"_fastqc.html")
#     dt$run=paste0("<a href='",link1_name,"'target='_blank'>",dt$run,"_1","</a>","\n","<a href='",link2_name,"'target='_blank'>",dt$run,"_2","</a>")
    dt=as.data.frame(dt)
    dt=dt %>% dplyr::select(id,Numberofinputreads,percentageofreadsunmapped_toomanymismatches:ofchimericreads)
    return(dt)
  })
  #~~~~~~~~~~~~~~~~~~~~
  output$table_unique = DT::renderDataTable({
    DT::datatable(table_unique(),
                  extensions =  c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),selection = list(mode = 'single', selected =1),caption= "Unique Reads",escape=FALSE)
  })
          
          output$table_unmapped = DT::renderDataTable({
            DT::datatable(table_unmapped(),
                          extensions =  c('Buttons','Scroller'),
                          options = list(dom = 'Bfrtip',
                            searchHighlight = TRUE,
                            pageLength = 10,
                            lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                            scrollX = TRUE,
                            buttons = c('copy','print')
                          ),selection = list(mode = 'single', selected =1),caption= "Unmapped Reads",escape=FALSE)
          })
          
          output$table_multi = DT::renderDataTable({
            DT::datatable(table_multi(),
                          extensions =  c('Buttons','Scroller'),
                          options = list(dom = 'Bfrtip',
                            searchHighlight = TRUE,
                            pageLength = 10,
                            lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                            scrollX = TRUE,
                            buttons = c('copy', 'print')
                          ),selection = list(mode = 'single', selected =1),caption= "Multi-Mapped Reads",escape=FALSE)
          })
  #~~~~~~~~~~~~~~~~~~~~
  #Get column list and populate drop-down for each tab
  output$ui_unique = renderUI({
    cols1=c("Uniquelymappedreads_percentage","Uniquelymappedreads_number","Numberofsplices_Total","Numberofsplices_GT_AG","Numberofsplices_GC_AG","Numberofsplices_AT_AC")
    selectInput("attr1","Select an attribute for the boxplot",as.list(as.character(cols1)))
  })
  
  output$ui_multi = renderUI({
    cols2=c("Numberofreadsmappedtomultipleloci", "percentageofreadsmappedtomultipleloci","Numberofreadsmappedtotoomanyloci","percentageofreadsmappedtotoomanyloci")
    selectInput("attr2","Select an attribute for the boxplot",as.list(as.character(cols2)))
  })
  
  output$ui_unmapped = renderUI({
    cols3=c("percentageofreadsunmapped_toomanymismatches","percentageofreadsunmapped_tooshort","percentageofreadsunmapped_other")
    selectInput("attr3","Select an attribute for the boxplot",as.list(as.character(cols3)))
  })
  #~~~~~~~~~~~~~~~~~~~~
  #Generate box-plot (uniquely mapped)
  boxplot1_out = reactive({
    d=readfile()
    attr1=input$attr1
    v=paste("d$",attr1,sep="")
    v=eval(parse(text=v))
    p <- plot_ly(d,x=~id,y=~v,color=~id) %>%
      layout(title = "BOX PLOT",xaxis = list(title ="ID"),yaxis = list(title = as.character(attr1)),margin=list(b=120,pad=4))
  })
  #Generate box-plot (multi-mapped)
  boxplot2_out = reactive({
    d=readfile()
    attr2=input$attr2
    v=paste("d$",attr2,sep="")
    v=eval(parse(text=v))
    p <- plot_ly(d,x=~id,y=~v,color=~id) %>%
      layout(title = "BOX PLOT",xaxis = list(title ="ID"),yaxis = list(title = as.character(attr2)),margin=list(b=120,pad=4)) 

  })
  #Generate box-plot (unmapped)
  boxplot3_out = reactive({
    d=readfile()
    attr3=input$attr3
    v=paste("d$",attr3,sep="")
    v=eval(parse(text=v))
    p <- plot_ly(d,x=~id,y=~v,color=~id) %>%
      layout(title = "BOX PLOT",xaxis = list(title ="ID"),yaxis = list(title = as.character(attr3)),margin=list(b=120,pad=4))
  })


  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  #########################DISPLAY Merged data libcomplex and metrics data #########################################
  ##################################################################################################################
  ##################################################################################################################
# 
  #Read file based on id selected
  libcomplex = reactive({
    STAR=fileload()
    df =STAR$libcomplex
    validate(
      need(nrow(df) != 0, "Information not available")
    )
return(df)
  })
  
  metrics = reactive({
    STAR=fileload()
    df2=STAR$metrics
    validate(
      need(nrow(df2) != 0, "Information not available")
    )
    df2$Sample=rownames(df2)
    df=as.data.frame(df2)
    rownames(df)=df$Sample
    df=as.data.frame(df)
  })
  
  mrkdup = reactive({
    STAR=fileload()
    df=STAR$markdups
    validate(
      need(nrow(df) != 0, "Information not available")
    )
    df$Sample=rownames(df)
    df=as.data.frame(df)
    rownames(df)=df$Sample
    return(df)
  })
  
  #print raw star summary report in tab3
  output$libcomplex = DT::renderDataTable({
    DT::datatable(libcomplex(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=TRUE,caption= "Library complex")
  })
  output$mrkdup = DT::renderDataTable({
    DT::datatable(mrkdup(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=TRUE,caption= "Mark duplicates")
  })
  
  output$metrics = DT::renderDataTable({
    DT::datatable(metrics(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=TRUE,caption= "Metrics")
  })

  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  #########################DISPLAY annotation data #########################################
  ##################################################################################################################
  ##################################################################################################################

  anno = reactive({
    v = fileload()
    validate(
      need(nrow(v$pData) !=0, "No phenotype information. Please check RData for pheno file")
    )
    pData<-v$pData
    return(pData)
    
  })
  
  output$anno = DT::renderDataTable({
    DT::datatable(anno(),
                  extensions = c('Buttons','Scroller'),
                  options = list(dom = 'Bfrtip',
                                 searchHighlight = TRUE,
                                 pageLength = 10,
                                 lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                 scrollX = TRUE,
                                 buttons = c('copy', 'print')
                  ),rownames=FALSE,caption= "PhenoData")
  })
  

#   ##################################################################################################################
#   ##################################################################################################################
#   ##################################################################################################################
#   #Populate dropdowns for distribution
  output$yoptions <- renderUI({
    d<- libcomplex()
    d=d %>% dplyr::select(-Sample,-LIBRARY)
    fac=colnames(d)
    selectInput("yoptions", "Select y attribute",as.list(as.character(fac)))
  })
  
  output$xoptions <- renderUI({
    fac=c("Sample","Library_Pool","Tissue_Source","CHF_Etiology","Gender","Race","Afib","VTVF","Diabetes","Hypertension","Random_Pool","RIN")
    selectInput("xoptions", "Select x attribute",as.list(as.character(fac)))
  })
  
  output$mxoptions <- renderUI({
    fac=c("Sample","Library_Pool","Tissue_Source","CHF_Etiology","Gender","Race","Afib","VTVF","Diabetes","Hypertension","Random_Pool","RIN")
    selectInput("mxoptions", "Select x attribute",as.list(as.character(fac)))
  })
  
  output$myoptions <- renderUI({
    d<- metrics()
    d=d %>% dplyr::select(-Sample)
    fac=colnames(d)
    selectInput("myoptions", "Select y attribute",as.list(as.character(fac)))
  })
  
  output$dxoptions <- renderUI({
    #pheno=readpheno()
    #run=sort(unique(d$run))
    #fac=c("Library_Pool","Tissue_Source","CHF_Etiology","Gender","Race","Afib","VTVF","Diabetes","Hypertension","Random_Pool","RIN")
    #fac=c("Library_Pool","Tissue_Source","CHF_Etiology","Gender","Race","Afib","VTVF","Diabetes","Hypertension","Random_Pool","RIN")
    fac=c("Sample","Library_Pool","Tissue_Source","CHF_Etiology","Gender","Race","Afib","VTVF","Diabetes","Hypertension","Random_Pool","RIN")
    selectInput("dxoptions", "Select x attribute",as.list(as.character(fac)))
  })
  
  output$dyoptions <- renderUI({
    d<- mrkdup()
    #run=sort(unique(d$run))
    d=d %>% dplyr::select(-Sample,-LIBRARY)
    fac=colnames(d)
    selectInput("dyoptions", "Select y attribute",as.list(as.character(fac)))
  })
  
  #create bar graph for selected samples
  libcplot = reactive({
    d= libcomplex()
    #pheno=readpheno()
    #all <- inner_join(d,pheno,by='Sample') %>% dplyr::mutate(Pool=paste(Flowcell,Lane,sep='_'))
    all <-d
    xoptions=input$xoptions
    yoptions=input$yoptions
    yop=paste("all$",yoptions,sep="")
    yop=eval(parse(text=yop))
    xop=paste("all$",xoptions,sep="")
    xop=eval(parse(text=xop))
    p <- plot_ly(all,x=xop,y=yop,color=xop) %>%
    layout(title = "BOX PLOT",xaxis = list(title =as.character(input$xoptions)),yaxis = list(title = as.character(input$yoptions)))
    (gg <- ggplotly(p))
  })

  
  metricsplot = reactive({
    d= metrics()
    all=d
    xoptions=input$mxoptions
    yoptions=input$myoptions
    yop=paste("all$",yoptions,sep="")
    yop=eval(parse(text=yop))
    xop=paste("all$",xoptions,sep="")
    xop=eval(parse(text=xop))
    p <- plot_ly(all,x=xop,y=yop,color=xop) %>%
      layout(title = "BOX PLOT",xaxis = list(title =as.character(input$mxoptions)),yaxis = list(title = as.character(input$myoptions)))
    (gg <- ggplotly(p))
  })
  
  mrkdupsplot = reactive({
    d= mrkdup()
    all=d
#     pheno=readpheno()
#     all <- inner_join(d,pheno,by='Sample') %>% dplyr::mutate(Pool=paste(Flowcell,Lane,sep='_'))
    xoptions=input$dxoptions
    yoptions=input$dyoptions
    yop=paste("all$",yoptions,sep="")
    yop=eval(parse(text=yop))
    xop=paste("all$",xoptions,sep="")
    xop=eval(parse(text=xop))
    p <- plot_ly(all,x=xop,y=yop,color=xop) %>%
      layout(title = "BOX PLOT",xaxis = list(title =as.character(input$dxoptions)),yaxis = list(title = as.character(input$dyoptions)))
    (gg <- ggplotly(p))
  })
  output$libc_bplot = renderPlotly({
    libcplot()
  })
  
  output$metr_bplot = renderPlotly({
    metricsplot()
  })
  
  output$mrkdup_bplot = renderPlotly({
    mrkdupsplot()
  })

  ##################################################################################################################
  ##################################################################################################################
  ##################################################################################################################
  ###############################################   USER-InTERFACE ################################################
  ##################################################################################################################
  ##################################################################################################################
  #User-Interface

  output$plotUI = renderUI({
    do.call(tabsetPanel,
            lapply(input$stats,function(s){ #for either upregulated/downregulated selected
              call("tabPanel",s,
                   call('uiOutput',paste0("ui_",s)),
                   call("plotlyOutput",paste0("boxplot_",s),width="1000",height="500"),hr(),
                   call('dataTableOutput',paste0("table_",s)),
                   call("downloadButton",paste0("save_",s),'Save as csv')
              )
            })
    )
  })
#   output$plotUI = renderUI({
#     fac=c("apoo","prachu                                                                                                                                                                            ")
#     selectInput("kt", "Select y attribute",as.list(as.character(fac)))
#   })
  
  observe({
    
    lapply(input$stats, function(s){
      ## Add a DataTable 
      output[['table_unique']] <- DT::renderDataTable(table_unique(),options=list(iDisplayLength=10))
      output[['table_multi']] <- DT::renderDataTable(table_multi(),options=list(iDisplayLength=10))
      output[['table_unmapped']] <- DT::renderDataTable(table_unmapped(),options=list(iDisplayLength=10))
      
      
      output$table_unique = DT::renderDataTable({
        DT::datatable(table_unique(),
                      extensions = c('Buttons','Scroller'),
                      options = list(dom = 'Bfrtip',
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                     scrollX = TRUE,
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                      ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE,caption= "Unique Reads")
      })
      
      output$table_unmapped = DT::renderDataTable({
        DT::datatable(table_unmapped(),
                      extensions = c('Buttons','Scroller'),
                      options = list(dom = 'Bfrtip',
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                     scrollX = TRUE,
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                      ),rownames=TRUE,selection = list(mode = 'single', selected =1),escape=FALSE,caption= "Unmapped Reads")
      })
      
      output$table_multi = DT::renderDataTable({
        DT::datatable(table_multi(),
                      extensions = c('Buttons','Scroller'),
                      options = list(dom = 'Bfrtip',
                                     searchHighlight = TRUE,
                                     pageLength = 10,
                                     lengthMenu = list(c(30, 50, 100, 150, 200, -1), c('30', '50', '100', '150', '200', 'All')),
                                     scrollX = TRUE,
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                      ),selection = list(mode = 'single', selected =1),caption= "Multi-Mapped Reads",escape=FALSE)
      })
      
      ## Add the boxplots
      output[['boxplot_unique']] = renderPlotly({boxplot1_out()})
      output[['boxplot_multi']] = renderPlotly({boxplot2_out()})
      output[['boxplot_unmapped']] = renderPlotly({boxplot3_out()})
      
      
      output[['save_unique']] <- downloadHandler(
        filename = function() { paste0("unique",".csv") },
        content = function(file) { write.csv(table_unique(),file=file,row.names=FALSE) }
      )
      output[['save_unmapped']] <- downloadHandler(
        filename = function() { paste0("unmapped",".csv") },
        content = function(file) { write.csv(table_unmapped(),file=file,row.names=FALSE) }
      )
      output[['save_multi']] <- downloadHandler(
        filename = function() { paste0("multi",".csv") },
        content = function(file) { write.csv(table_multi(),file=file,row.names=FALSE) }
      )
      
      return(s)
    })
  })
 
  
  
})
