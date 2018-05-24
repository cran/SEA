server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=-1)
  
  
  upda1<-observeEvent(input$Resolution,{
    
    if(input$Resolution=="General resolution"){
      widthG<-960
      heightG<-600
      pointG<-20
      ppiG<-72 
    }else if(input$Resolution=="High resolution"){
      widthG<-10000
      heightG<-6000
      pointG<-30
      ppiG<-300 
    }
    updateTextInput(session, "wid", value=widthG)
    updateTextInput(session, "hei", value=heightG)
    updateTextInput(session, "ppi1", value=pointG)
    updateTextInput(session, "ppi2", value=ppiG)
  }) 
  
  
  
  observeEvent(input$RunSo, {
    updateTabsetPanel(session, "inTabset",
                      selected="RE"
                      
    )
  })
  
  
  observeEvent(input$PoPr, {
    updateTabsetPanel(session, "inTabset",
                      selected="POP"
    )
  })
  
  observeEvent(input$drawPl, {
    updateTabsetPanel(session, "inTabset",
                      selected="DC"
    )
  })
  
  manual<-eventReactive(input$manl,{
    RShowDoc("Instruction",package="SEA") 
  })
  output$manll<-renderUI(manual())
  
  
  ReData<-reactive({
    req(input$fileDataset$datapath)
    dataRaw<-fread(input$fileDataset$datapath,header = FALSE,stringsAsFactors=T)
  })
  
  data<-reactive({
    dataRaw<-ReData()
    data_show<-dataRaw[-1,]
    colnames(data_show)<-as.matrix(dataRaw[1,])
    as.data.frame(data_show)
  })
  
  output$datashow<-renderDataTable({
    data()    
  })
  
  SEA<-eventReactive(input$RunSo,{
    
    id <- showNotification("Calculation in progress, please be patient...", duration = NULL,type="message")
    data<-ReData()
    
    result=switch(input$PopulationType,"F2" = F2Fun(data,input$F2select1),
                  "F2:3"=F23Fun(data,input$F3select1,as.numeric(input$F23text2)),
                  "DH"=DHFun(data,input$DHselect1),
                  "BIL"=BILFun(data,input$BILselect1,input$BILfr),
                  "BC (B1 B2)"=BCFun(data,input$BCselect1),
                  "BCF (B1:2 B2:2)"=BCFFun(data,input$BCFselect1,as.numeric(input$BCFtext2)),
                  "G4F2 (P1 P2 F1 F2)"=G4F2Fun(data,input$G4F2select1),
                  "G4F3 (P1 P2 F1 F2:3)"=G4F3Fun(data,input$G4F3select1,as.numeric(input$G4F3text2)),
                  "G3DH (P1 P2 DH)"= G3DHFun(data,input$G3DHselect1,as.numeric(input$G3DHtext2)),
                  "G5BC (P1 P2 F1 B1 B2)"=G5BCFun(data,input$G5BCselect1),
                  "G5BCF (P1 P2 F1 B1:2 B2:2)"=G5BCFFun(data,input$G5BCFselect1,as.numeric(input$G5BCFtext2)),
                  "G5 (P1 P2 F1 F2 F2:3)"=G5Fun(data,input$G5select1,as.numeric(input$G5text2)),
                  "G6 (P1 P2 F1 F2 B1 B2)"=G6Fun(data,input$G6select1),
                  "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=G6FFun(data,input$G6Fselect1,as.numeric(input$G6Ftext2)))
    removeNotification(id)   
    result
  })
  
  output$RunResult<-renderDataTable({
    as.data.frame(SEA()[[1]])
  })
  
  output$downloadresult<-downloadHandler(
    filename = function(){paste("result", ".csv", sep = "")},
    content = function(file) {
      results<-as.data.frame(SEA()[[1]])
      fwrite(results, file, row.names = FALSE)
    }
  )
  
  modelup<-reactive({
    
    modelselect1=switch(input$PopulationType,"F2" = input$F2select1,
                        "F2:3"=input$F3select1,
                        "DH"=input$DHselect1,
                        "BIL"=input$BILselect1,
                        "BC (B1 B2)"=input$BCselect1,
                        "BCF (B1:2 B2:2)"=input$BCFselect1,
                        "G4F2 (P1 P2 F1 F2)"=input$G4F2select1,
                        "G4F3 (P1 P2 F1 F2:3)"=input$G4F3select1,
                        "G3DH (P1 P2 DH)"= input$G3DHselect1,
                        "G5BC (P1 P2 F1 B1 B2)"=input$G5BCselect1,
                        "G5BCF (P1 P2 F1 B1:2 B2:2)"=input$G5BCFselect1,
                        "G5 (P1 P2 F1 F2 F2:3)"=input$G5select1,
                        "G6 (P1 P2 F1 F2 B1 B2)"=input$G6select1,
                        "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=input$G6Fselect1)
    
    plmodelselect2=switch(input$PopulationType,"F2" = "plF2select2",
                          "F2:3"="plF3select2",
                          "DH"="plDHselect2",
                          "BIL"="plBILselect2",
                          "BC (B1 B2)"="plBCselect2",
                          "BCF (B1:2 B2:2)"="plBCFselect2",
                          "G4F2 (P1 P2 F1 F2)"="plG4F2select2",
                          "G4F3 (P1 P2 F1 F2:3)"="plG4F3select2",
                          "G3DH (P1 P2 DH)"= "plG3DHselect2",
                          "G5BC (P1 P2 F1 B1 B2)"="plG5BCselect2",
                          "G5BCF (P1 P2 F1 B1:2 B2:2)"="plG5BCFselect2",
                          "G5 (P1 P2 F1 F2 F2:3)"="plG5select2",
                          "G6 (P1 P2 F1 F2 B1 B2)"="plG6select2",
                          "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"="plG6Fselect2") 
    
    
    modelse<-list(modelselect1,plmodelselect2) 
    
  })
  
  
  upda2<-observe({
    
    resultu<-as.matrix(SEA()[[1]])
    
    OrPp<-order(as.numeric(resultu[,3]))
    
    upvalue<-resultu[OrPp,1][1:5]
    
    modelsel<-modelup()
    
    modelselect1<-modelsel[[1]];plmodelselect2<-modelsel[[2]]
    
    
    modelselect2=switch(input$PopulationType,"F2" = "F2select2",
                        "F2:3"="F3select2",
                        "DH"="DHselect2",
                        "BIL"="BILselect2",
                        "BC (B1 B2)"="BCselect2",
                        "BCF (B1:2 B2:2)"="BCFselect2",
                        "G4F2 (P1 P2 F1 F2)"="G4F2select2",
                        "G4F3 (P1 P2 F1 F2:3)"="G4F3select2",
                        "G3DH (P1 P2 DH)"= "G3DHselect2",
                        "G5BC (P1 P2 F1 B1 B2)"="G5BCselect2",
                        "G5BCF (P1 P2 F1 B1:2 B2:2)"="G5BCFselect2",
                        "G5 (P1 P2 F1 F2 F2:3)"="G5select2",
                        "G6 (P1 P2 F1 F2 B1 B2)"="G6select2",
                        "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"="G6Fselect2")
    
    if(modelselect1=="All models"){
      updateSelectInput(session, modelselect2, choices=upvalue)
      updateSelectInput(session, plmodelselect2, choices=upvalue)
    }else{
      updateSelectInput(session, modelselect2, choices=modelselect1)
      updateSelectInput(session, plmodelselect2, choices=modelselect1)
    }
  })
  
  
  SEAPP<-eventReactive(input$CaPop,{
    
    id <- showNotification("Calculation in progress,please be patient...", duration = NULL,type = "message")
    data<-ReData()
    modelsel<-modelup();model<-modelsel[[1]];modelselect22<-modelsel[[2]]
    
    if(model!="All models"){
      resultP<-SEA()
    }else{
      resultP=switch(input$PopulationType,"F2" = F2Fun(data,input$F2select2),
                     "F2:3"=F23Fun(data,input$F3select2,as.numeric(input$F23text2)),
                     "DH"=DHFun(data,input$DHselect2),
                     "BIL"=BILFun(data,input$BILselect2,input$BILfr),
                     "BC (B1 B2)"=BCFun(data,input$BCselect2),
                     "BCF (B1:2 B2:2)"=BCFFun(data,input$BCFselect2,as.numeric(input$BCFtext2)),
                     "G4F2 (P1 P2 F1 F2)"=G4F2Fun(data,input$G4F2select2),
                     "G4F3 (P1 P2 F1 F2:3)"=G4F3Fun(data,input$G4F3select2,as.numeric(input$G4F3text2)),
                     "G3DH (P1 P2 DH)"= G3DHFun(data,input$G3DHselect2,as.numeric(input$G3DHtext2)),
                     "G5BC (P1 P2 F1 B1 B2)"=G5BCFun(data,input$G5BCselect2),
                     "G5BCF (P1 P2 F1 B1:2 B2:2)"=G5BCFFun(data,input$G5BCFselect2,as.numeric(input$G5BCFtext2)),
                     "G5 (P1 P2 F1 F2 F2:3)"=G5Fun(data,input$G5select2,as.numeric(input$G5text2)),
                     "G6 (P1 P2 F1 F2 B1 B2)"=G6Fun(data,input$G6select2),
                     "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=G6FFun(data,input$G6Fselect2,as.numeric(input$G6Ftext2)))
    }
    
    modelselect22=switch(input$PopulationType,"F2" = input$F2select2,
                         "F2:3"=input$F3select2,
                         "DH"=input$DHselect2,
                         "BIL"=input$BILselect2,
                         "BC (B1 B2)"=input$BCselect2,
                         "BCF (B1:2 B2:2)"=input$BCFselect2,
                         "G4F2 (P1 P2 F1 F2)"=input$G4F2select2,
                         "G4F3 (P1 P2 F1 F2:3)"=input$G4F3select2,
                         "G3DH (P1 P2 DH)"= input$G3DHselect2,
                         "G5BC (P1 P2 F1 B1 B2)"=input$G5BCselect2,
                         "G5BCF (P1 P2 F1 B1:2 B2:2)"=input$G5BCFselect2,
                         "G5 (P1 P2 F1 F2 F2:3)"=input$G5select2,
                         "G6 (P1 P2 F1 F2 B1 B2)"=input$G6select2,
                         "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=input$G6Fselect2)
    
    if(modelselect22!="0MG"&&modelselect22!="PG-AD"&&modelselect22!="PG-ADI"&&modelselect22!="PG-AI"&&modelselect22!="PG-A"){
      resultpp<-PosPro(input$PopulationType,resultP,data)
      removeNotification(id)
      as.data.frame(resultpp)
    }else{
      showModal(modalDialog(title = "Warning", strong("No posterior probability!"), easyClose = TRUE)) 
      removeNotification(id)
    }
  })
  
  output$RePop<-renderDataTable({
    SEAPP()
  })
  
  output$PopSave<-downloadHandler(
    filename = function(){paste("Posterior Probability", ".csv", sep = "")},
    content = function(file) {
      popresults<-SEAPP()
      fwrite(popresults, file, row.names = FALSE)
    }
  )
  
  
  SEAPl<-eventReactive(input$DrPl,{
    
    data<-ReData()
    
    resultPl=switch(input$PopulationType,"F2" = F2Fun(data,input$plF2select2),
                    "F2:3"=F23Fun(data,input$plF3select2,as.numeric(input$F23text2)),
                    "DH"=DHFun(data,input$plDHselect2),
                    "BIL"=BILFun(data,input$plBILselect2,input$BILfr),
                    "BC (B1 B2)"=BCFun(data,input$plBCselect2),
                    "BCF (B1:2 B2:2)"=BCFFun(data,input$plBCFselect2,as.numeric(input$BCFtext2)),
                    "G4F2 (P1 P2 F1 F2)"=G4F2Fun(data,input$plG4F2select2),
                    "G4F3 (P1 P2 F1 F2:3)"=G4F3Fun(data,input$plG4F3select2,as.numeric(input$G4F3text2)),
                    "G3DH (P1 P2 DH)"= G3DHFun(data,input$plG3DHselect2,as.numeric(input$G3DHtext2)),
                    "G5BC (P1 P2 F1 B1 B2)"=G5BCFun(data,input$plG5BCselect2),
                    "G5BCF (P1 P2 F1 B1:2 B2:2)"=G5BCFFun(data,input$plG5BCFselect2,as.numeric(input$G5BCFtext2)),
                    "G5 (P1 P2 F1 F2 F2:3)"=G5Fun(data,input$plG5select2,as.numeric(input$G5text2)),
                    "G6 (P1 P2 F1 F2 B1 B2)"=G6Fun(data,input$plG6select2),
                    "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=G6FFun(data,input$plG6Fselect2,as.numeric(input$G6Ftext2)))
    
    Population2=switch(input$PopulationType,
                       "F2" = NULL,
                       "F2:3"=NULL,
                       "DH"=NULL,
                       "BIL"=NULL,
                       "G4F2 (P1 P2 F1 F2)"=NULL,
                       "G4F3 (P1 P2 F1 F2:3)"=NULL,
                       "G3DH (P1 P2 DH)"=NULL,
                       "BC (B1 B2)"= input$BCPlotSelect,
                       "BCF (B1:2 B2:2)"=input$BCFPlotSelect,
                       "G5BC (P1 P2 F1 B1 B2)"=input$G5BCPlotSelect,
                       "G5BCF (P1 P2 F1 B1:2 B2:2)"=input$G5BCFPlotSelect,
                       "G5 (P1 P2 F1 F2 F2:3)"=input$G5PlotSelect,
                       "G6 (P1 P2 F1 F2 B1 B2)"=input$G6PlotSelect,
                       "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=input$G6FPlotSelect)
    
    
    plotfun(input$PopulationType,Population2,data,resultPl,as.numeric(input$bins),input$pr,input$colour)
  })
  
  
  output$SEAPlot<-renderPlot({
    
    SEAPl()
    
  })
  
  
  output$downloadplot <- downloadHandler(
    filename = function() {
      paste("plot", sep = ".", switch(
        input$plformat, "*.png"=".png", "*.tiff"=".tiff", "*.jpeg"=".jpeg","*.pdf=.pdf"
      ))
    },
    content = function(file) {
      if(input$plformat=="*.png"){
        png(file,width=as.numeric(input$wid), height=as.numeric(input$hei), units= "px", pointsize =as.numeric(input$ppi1),res=as.numeric(input$ppi2))
      }else if(input$plformat=="*.tiff"){
        tiff(file,width=as.numeric(input$wid), height=as.numeric(input$hei), units= "px", pointsize =as.numeric(input$ppi1),res=as.numeric(input$ppi2))
      }else if(input$plformat=="*.jpeg"){
        jpeg(file,width=as.numeric(input$wid), height=as.numeric(input$hei), units= "px", pointsize =as.numeric(input$ppi1),res=as.numeric(input$ppi2))
      }else if(input$plformat=="*.pdf"){
        pdf(file,width=12)
        
      }
      data<-ReData()
      
      resultPl=switch(input$PopulationType,"F2" = F2Fun(data,input$plF2select2),
                      "F2:3"=F23Fun(data,input$plF3select2,as.numeric(input$F23text2)),
                      "DH"=DHFun(data,input$plDHselect2),
                      "BIL"=BILFun(data,input$plBILselect2,input$BILfr),
                      "BC (B1 B2)"=BCFun(data,input$plBCselect2),
                      "BCF (B1:2 B2:2)"=BCFFun(data,input$plBCFselect2,as.numeric(input$BCFtext2)),
                      "G4F2 (P1 P2 F1 F2)"=G4F2Fun(data,input$plG4F2select2),
                      "G4F3 (P1 P2 F1 F2:3)"=G4F3Fun(data,input$plG4F3select2,as.numeric(input$G4F3text2)),
                      "G3DH (P1 P2 DH)"= G3DHFun(data,input$plG3DHselect2,as.numeric(input$G3DHtext2)),
                      "G5BC (P1 P2 F1 B1 B2)"=G5BCFun(data,input$plG5BCselect2),
                      "G5BCF (P1 P2 F1 B1:2 B2:2)"=G5BCFFun(data,input$plG5BCFselect2,as.numeric(input$G5BCFtext2)),
                      "G5 (P1 P2 F1 F2 F2:3)"=G5Fun(data,input$plG5select2,as.numeric(input$G5text2)),
                      "G6 (P1 P2 F1 F2 B1 B2)"=G6Fun(data,input$plG6select2),
                      "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=G6FFun(data,input$plG6Fselect2,as.numeric(input$G6Ftext2)))
      
      Population2=switch(input$PopulationType,
                         "F2" = NULL,
                         "F2:3"=NULL,
                         "DH"=NULL,
                         "BIL"=NULL,
                         "G4F2 (P1 P2 F1 F2)"=NULL,
                         "G4F3 (P1 P2 F1 F2:3)"=NULL,
                         "G3DH (P1 P2 DH)"=NULL,
                         "BC (B1 B2)"= input$BCPlotSelect,
                         "BCF (B1:2 B2:2)"=input$BCFPlotSelect,
                         "G5BC (P1 P2 F1 B1 B2)"=input$G5BCPlotSelect,
                         "G5BCF (P1 P2 F1 B1:2 B2:2)"=input$G5BCFPlotSelect,
                         "G5 (P1 P2 F1 F2 F2:3)"=input$G5PlotSelect,
                         "G6 (P1 P2 F1 F2 B1 B2)"=input$G6PlotSelect,
                         "G6F (P1 F1 P2 B1:2 B2:2 F2:3)"=input$G6FPlotSelect)
      
      
      plotfun(input$PopulationType,Population2,data,resultPl,as.numeric(input$bins),input$pr,input$colour)
      dev.off()
      
    })
  
}