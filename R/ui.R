ui <-tagList(
  navbarPage(
    "",id = "tabs",
    tabPanel(strong("SEA"),
             h2("SEA (Segregation Analysis)",align="center"),
             h2("or",align="center"),
             h2("Mixed major-genes plus polygenes inheritance analysis",align="center"),
             br(),
             h4("Version 2.0, Realeased February 2021",align="center"),
             br(),
             column(10,
                    h4("1. This second version of the software (SEA.R) was modified and validated via Monte Carlo 
					   simulation studies by Mr Wang Jing-Tian in my lab, at College of Plant Science and Technology
					   of Huazhong Agricultural University, from the 1st version of SEA.R."),
                    br(),
                    h4("2. The first version of SEA.R was developed by Mr Zhang Ya-Wen, Ms Du Ying-Wen, and Dr Ren
                       Wen-Long in my lab, at College of Plant Science and Technology of Huazhong Agricultural 
					   University, from the VC++ codes of mixed major genes plus polygenes inheritance analysis of
					   quantitative traits, which also named SEgregation Analysis (SEA). Registration number for this 
					   software copyright in China: 2017SR579357."),
                    br(),
                    h4(
                      "3. The VC++ codes for the SEA software was developed by Mr Liu Bing and Mr Cao Xi-wen in my 
					  lab at College of Agriculture of Nanjing Agricultural University, based on the turbo C++ 
					  version of the SEA software. The copyright registration numbers for several software programs 
					  had been obtained in China."),
                    br(),
                    h4("4. The turbo C++ version of the SEA software was developed by Dr Yuan-Ming Zhang. This software
					   is based on the EIM algorithm proposed by Zhang et al. Genetical Research 2003, 81(2): 157-163. 
					   This version is widely used in China. Individual codes are derived from the EM algorithm based 
					   SEA codes, written by Dr Wang Jiankang."),
                    br(),
                    offset=1)
                    ),
    tabPanel(strong("Start"),
             titlePanel("SEA (Segregation Analysis)"),
             sidebarLayout(
               sidebarPanel(
                 selectInput("PopulationType","Select population:",choices=c("F2","F2:3","DH","BIL","BC (B1 B2)","BCF (B1:2 B2:2)","G4F2 (P1 P2 F1 F2)","G4F3 (P1 P2 F1 F2:3)","G3DH (P1 P2 DH)","G5BC (P1 P2 F1 B1 B2)","G5BCF (P1 P2 F1 B1:2 B2:2)","G5 (P1 P2 F1 F2 F2:3)","G6 (P1 P2 F1 F2 B1 B2)","G6F (P1 F1 P2 B1:2 B2:2 F2:3)")),
                 fileInput("fileDataset", "Input dataset",multiple = TRUE),
                 conditionalPanel("input.PopulationType == 'F2'",
                                  selectInput(inputId = "F2select1",label = "Model Selection",
                                              choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'F2:3'",
                                  selectInput(inputId = "F3select1",label = "Model Selection",
                                              choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","All models"),selected = "All models"),
                                  textInput("F23text2", label = "No. of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'DH'",
                                  selectInput(inputId = "DHselect1",label = "Model Selection",
                                              choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'BIL'",
                                  selectInput(inputId = "BILselect1",label = "Model Selection",
                                              choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","2MG-IE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","All models"),selected = "All models"),
                                  selectInput("BILfr", label = "The BIL type",choices = c("BIL1(F1xP1)","BIL2(F1xP2)"))
                 ),
                 conditionalPanel("input.PopulationType == 'BC (B1 B2)'",
                                  selectInput(inputId = "BCselect1",label = "Model Selection",
                                              choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'BCF (B1:2 B2:2)'",
                                  selectInput(inputId = "BCFselect1",label = "Model Selection",
                                              choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","All models"),selected = "All models"),
                                  textInput("BCFtext2", label = "No.of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'G4F2 (P1 P2 F1 F2)'",
                                  selectInput(inputId = "G4F2select1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'G4F3 (P1 P2 F1 F2:3)'",
                                  selectInput(inputId = "G4F3select1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models"),
                                  textInput("G4F3text2", label = "No. of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'G3DH (P1 P2 DH)'",
                                  selectInput(inputId = "G3DHselect1",label = "Model Selection",
                                              choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","2MG-IE","PG-AI","PG-A","MX1-A-AI","MX1-A-A","MX2-AI-AI","MX2-AI-A","MX2-A-A","MX2-EA-A","MX2-ED-A","MX2-ER-A","MX2-AE-A","MX2-CE-A","MX2-DE-A","MX2-IE-A","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","MX3-AI-AI","MX3-AI-A","MX3-A-A","MX3-CEA-A","MX3-PEA-A","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA","All models"),selected = "All models"),
                                  textInput("G3DHtext2", label = "No. of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'G5BC (P1 P2 F1 B1 B2)'",
                                  selectInput(inputId = "G5BCselect1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'G5BCF (P1 P2 F1 B1:2 B2:2)'",
                                  selectInput(inputId = "G5BCFselect1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models"),
                                  textInput("G5BCFtext2", label = "No.of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'G5 (P1 P2 F1 F2 F2:3)'",
                                  selectInput(inputId = "G5select1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models"),
                                  textInput("G5text2", label = "No.of plants measured in each family",value = "1")
                 ),
                 conditionalPanel("input.PopulationType == 'G6 (P1 P2 F1 F2 B1 B2)'",
                                  selectInput(inputId = "G6select1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models")
                 ),
                 conditionalPanel("input.PopulationType == 'G6F (P1 F1 P2 B1:2 B2:2 F2:3)'",
                                  selectInput(inputId = "G6Fselect1",label = "Model Selection",
                                              choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD","All models"),selected = "All models"),
                                  textInput("G6Ftext2", label = "No.of plants measured in each family",value = "1")
                 ),
                 br(),
                 actionButton("RunSo", label = "Run",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 br(),
                 br(),
                 actionButton("PoPr", label = "Posterior Probability",width = 250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 br(),
                 actionButton("drawPl", label = "Distribution curves",width = 250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 br(),
                 actionButton("manl", label = "User manual",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 uiOutput("manll")
               ),
               mainPanel(
                 tabsetPanel(id="inTabset",
                             tabPanel("Dataset",value = "DA",

                                      dataTableOutput("datashow")
                             ),
                             tabPanel("Result",value = "RE",
                                      fluidRow(
                                        br(),
                                        column(12,
                                               h3(strong("Result")),
                                               br(),
                                               downloadButton("downloadresult", "Download result",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        ),
                                        column(12,
                                               br(),
                                               br(),
                                               dataTableOutput("RunResult")
                                        )
                                      )
                             ),
                             tabPanel("Posterior Probability",value = "POP",
                                      fluidRow(
                                        column(12,
                                               br(),
                                               h3(strong("Posterior Probability"))
                                        ),
                                        column(4,
                                               br(),
                                               conditionalPanel("input.PopulationType == 'F2'",
                                                                selectInput(inputId = "F2select2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'F2:3'",
                                                                selectInput(inputId = "F3select2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'DH'",
                                                                selectInput(inputId = "DHselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BIL'",
                                                                selectInput(inputId = "BILselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","2MG-IE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BC (B1 B2)'",
                                                                selectInput(inputId = "BCselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BCF (B1:2 B2:2)'",
                                                                selectInput(inputId = "BCFselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G4F2 (P1 P2 F1 F2)'",
                                                                selectInput(inputId = "G4F2select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G4F3 (P1 P2 F1 F2:3)'",
                                                                selectInput(inputId = "G4F3select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G3DH (P1 P2 DH)'",
                                                                selectInput(inputId = "G3DHselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","2MG-IE","PG-AI","PG-A","MX1-A-AI","MX1-A-A","MX2-AI-AI","MX2-AI-A","MX2-A-A","MX2-EA-A","MX2-ED-A","MX2-ER-A","MX2-AE-A","MX2-CE-A","MX2-DE-A","MX2-IE-A","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","MX3-AI-AI","MX3-AI-A","MX3-A-A","MX3-CEA-A","MX3-PEA-A","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BC (P1 P2 F1 B1 B2)'",
                                                                selectInput(inputId = "G5BCselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BCF (P1 P2 F1 B1:2 B2:2)'",
                                                                selectInput(inputId = "G5BCFselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5 (P1 P2 F1 F2 F2:3)'",
                                                                selectInput(inputId = "G5select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6 (P1 P2 F1 F2 B1 B2)'",
                                                                selectInput(inputId = "G6select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6F (P1 F1 P2 B1:2 B2:2 F2:3)'",
                                                                selectInput(inputId = "G6Fselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               )
                                        ),
                                        column(4,
                                               br(),
                                               br(),
                                               actionButton("CaPop", "Calculate Posterior Probability ",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")

                                        ),
                                        column(4,
                                               br(),
                                               br(),
                                               downloadButton("PopSave", "Save Posterior Probability ",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        ),
                                        column(12,
                                               br(),
                                               br(),
                                               br(),
                                               dataTableOutput("RePop")
                                        )

                                      )
                             ),
                             tabPanel("Distribution curves",value = "DC",
                                      fluidRow(
                                        column(3,
                                               br(),
                                               sliderInput(inputId = "bins",label = "Number of groups:",min = 1,max = 30,value = 10)
                                        ),
                                        column(3,
                                               br(),
                                               br(),
                                               conditionalPanel("input.PopulationType == 'F2'",
                                                                selectInput(inputId = "plF2select2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'F2:3'",
                                                                selectInput(inputId = "plF3select2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'DH'",
                                                                selectInput(inputId = "plDHselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE",
                                                                                        "2MG-DE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BIL'",
                                                                selectInput(inputId = "plBILselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE","2MG-DE","2MG-IE","3MG-AI","3MG-A","3MG-CEA","3MG-PEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BC (B1 B2)'",
                                                                selectInput(inputId = "plBCselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BCF (B1:2 B2:2)'",
                                                                selectInput(inputId = "plBCFselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G4F2 (P1 P2 F1 F2)'",
                                                                selectInput(inputId = "plG4F2select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G4F3 (P1 P2 F1 F2:3)'",
                                                                selectInput(inputId = "plG4F3select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G3DH (P1 P2 DH)'",
                                                                selectInput(inputId = "plG3DHselect2",label = "Optimal Model Selection",
                                                                            choices = c("0MG","1MG-A","2MG-AI","2MG-A","2MG-EA","2MG-ED","2MG-ER","2MG-AE","2MG-CE",
                                                                                        "2MG-DE","2MG-IE","PG-AI","PG-A","MX1-A-AI","MX1-A-A","MX2-AI-AI","MX2-AI-A","MX2-A-A","MX2-EA-A","MX2-ED-A","MX2-ER-A","MX2-AE-A","MX2-CE-A","MX2-DE-A","MX2-IE-A",
                                                                                        "3MG-AI","3MG-A","3MG-CEA","3MG-PEA","MX3-AI-AI","MX3-AI-A","MX3-A-A","MX3-CEA-A","MX3-PEA-A","4MG-AI","4MG-CEA","4MG-EEA","4MG-EEEA"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BC (P1 P2 F1 B1 B2)'",
                                                                selectInput(inputId = "plG5BCselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BCF (P1 P2 F1 B1:2 B2:2)'",
                                                                selectInput(inputId = "plG5BCFselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5 (P1 P2 F1 F2 F2:3)'",
                                                                selectInput(inputId = "plG5select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD",
                                                                                        "MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6 (P1 P2 F1 F2 B1 B2)'",
                                                                selectInput(inputId = "plG6select2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6F (P1 F1 P2 B1:2 B2:2 F2:3)'",
                                                                selectInput(inputId = "plG6Fselect2",label = "Optimal Model Selection",
                                                                            choices = c("1MG-AD","1MG-A","1MG-EAD","1MG-NCD","2MG-ADI","2MG-AD","2MG-A","2MG-EA","2MG-CD","2MG-EAD","PG-ADI","PG-AD","MX1-AD-ADI","MX1-AD-AD","MX1-A-AD","MX1-EAD-AD","MX1-NCD-AD","MX2-ADI-ADI","MX2-ADI-AD","MX2-AD-AD","MX2-A-AD","MX2-EA-AD","MX2-CD-AD","MX2-EAD-AD"))
                                               )
                                        ),
                                        column(3,
                                               br(),
                                               br(),
                                               conditionalPanel("input.PopulationType == 'BC (B1 B2)'",
                                                                selectInput(inputId = "BCPlotSelect",label = "Generation Selection",
                                                                            choices = c("B1","B2"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'BCF (B1:2 B2:2)'",
                                                                selectInput(inputId = "BCFPlotSelect",label = "Generation Selection",
                                                                            choices = c("B1:2","B2:2"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BC (P1 P2 F1 B1 B2)'",
                                                                selectInput(inputId = "G5BCPlotSelect",label = "Generation Selection",
                                                                            choices = c("B1","B2"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5BCF (P1 P2 F1 B1:2 B2:2)'",
                                                                selectInput(inputId = "G5BCFPlotSelect",label = "Generation Selection",
                                                                            choices = c("B1:2","B2:2"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G5 (P1 P2 F1 F2 F2:3)'",
                                                                selectInput(inputId = "G5PlotSelect",label = "Generation Selection",
                                                                            choices = c("F2","F2:3"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6 (P1 P2 F1 F2 B1 B2)'",
                                                                selectInput(inputId = "G6PlotSelect",label = "Generation Selection",
                                                                            choices = c("B1","B2","F2"))
                                               ),
                                               conditionalPanel("input.PopulationType == 'G6F (P1 F1 P2 B1:2 B2:2 F2:3)'",

                                                                selectInput(inputId = "G6FPlotSelect",label = "Generation Selection",
                                                                            choices = c("B1:2","B2:2","F2:3"))
                                               )
                                        ),
                                        column(12,
                                               radioButtons("moselect", " ", c("Parameter Settings", "Download plot"),inline = TRUE)
                                        ),
                                        conditionalPanel("input.moselect == 'Parameter Settings'",
                                                         column(3,
                                                                textInput("pr", label = "The number of intervals for right vertical axis",value = "5")
                                                         ),
                                                         column(3,
                                                                br(),
                                                                selectInput(inputId = "colour",label = "Curve colour",choices = c("red","black","blue","yellow","green","pink","purple","gray","brown"))
                                                         ),
                                                         column(3,
                                                                br(),
                                                                br(),
                                                                actionButton("DrPl", "To draw distribution curves",width=250,style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                         ),
                                                         column(8,
                                                                plotOutput("SEAPlot")
                                                         )
                                        ),
                                        conditionalPanel("input.moselect == 'Download plot'",
                                                         column(12,
                                                                radioButtons("Resolution", "Select resolution of plot", c("General resolution", "High resolution"),inline = TRUE)
                                                         ),
                                                         column(3,
                                                                textInput("wid", label = "Figure width",value = "960")
                                                         ),
                                                         column(3,
                                                                textInput("hei", label = "Figure height",value = "600")
                                                         ),
                                                         column(3,
                                                                textInput("ppi1", label = "Word resolution (1/72 inch, ppi):",value = "20")
                                                         ),
                                                         column(3,
                                                                textInput("ppi2", label = "Figure resolution (ppi):",value = "72")
                                                         ),
                                                         column(3,
                                                                br(),
                                                                radioButtons("plformat", "Plot format", c("*.png", "*.tiff", "*.jpeg","*.pdf"),inline = TRUE)
                                                         ),
                                                         column(3,
                                                                br(),
                                                                br(),
                                                                downloadButton("downloadplot", "Save Distribution curves",style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                                         )
                                        )
                                      )
                             )
                 )
               )
             )
    ),
    tabPanel(strong("Reference"),

             column(10,

                    h4("Please cite the below references for various types of populations.",align="center"),
                    br(),
                    br(),

                    h4("1. Zhang Yuan-Ming. Segregation analysis of quantitative traits and R software. Golden Light
                       Academic Publishing, 2017(All)"),
                    br(),
                    h4("2. Gai Junyi, Zhang Yuan-Ming, Wang Jiankang. Genetic system of quantitative traits in plants.
                       Beijing: Science Press, 2003.(All)"),
                    br(),
                    h4("3. Zhang Yuan-Ming, Gai Junyi, Yang Yonghua. The EIM algorithm in the joint segregation
                       analysis of quantitative traits. Genetical Research 2003, 81(2): 157-163.(All)"),
                    br(),
                    h4("4.Wang Jiankang, Gai Junyi. Identification of major gene and polygene mixed inheritance
                       model and estimation of genetic parameters of a quantitative trait from F2 Progeny. Acta Genetica
                       Sinica 1997, 24(5): 432-440.(F2)"),
                    br(),
                    h4("5.Zhang Yuan-Ming, Gai Junyi, Wang Jiankang. Identification of two major genes plus polygenes 
					   mixed inheritance model of quantitative traits in B1 and B2, and F2. Journal of Biomathematics 2000,
					   15(3):358-366.(F2 and BC (B1 B2))"),
                    br(),
                    h4("6.Zhang Yuan-Ming, Gai Junyi, Qi Cunkou. Detection of genetic system of quantitative traits using
					   backcross and selfing families. Hereditas (Beijing) 2001, 23(4):329-776.(F2:3)"),
                    br(),
                    h4("7.Zhang Yuan-Ming, Gai Junyi, Wang Yongjun. An expansion of joint segregation analysis of 
					   quantitative trait for using P1, P2 and DH or RIL populations. Hereditas (Beijing) 2001, 
					   23(5):467-470.(DH and G3DH (P1 P2 DH))"),
                    br(),
                    h4("8.Wang Jinshe, Zhao Tuanjie, Gai Junyi. Establishment of segregation analysis of mixed inheritance
            		   model with four major genes plus polygenes in backcross inbred lines (BIL) populations. Acta Agron 
					   Sin 2013, 39(2):198-206.(BIL and BCF (B1:2 B2:2))"),
                    br(),
                    h4("9.Zhang Yuan-Ming, Gai Junyi, Zhang Mengchen. Jointly segregating analysis of P1 P2 F1 and F2 or F2:3
					   families. Journal of Southwest Agricultural University 2000, 42(1):6-9.
					   (G4F2 (P1 F1 P2 F2) and G4F3 (P1 F1 P2 F2:3))"),
                    br(),
                    h4("10.Zhang Yuan-Ming, Gai Junyi. The IECM algorithm for estimation of component distribution parameters
 					   in segregating analysis of quantitative traits. Acta Agron Sin 2000, 26(6):699-706.
					   (G5BC (P1 P2 F1 B1 B2) and G5BCF (P1 P2 F1 B1:2 B2:2))"),
                    br(),
                    h4("11.Wang Jiankang, Gai Junyi. Identification of major gene and polygene mixed inheritance model of
					   quantitative traits by using joint analysis of P1, F1, P2, F2 and F2:3 generations. Acta Agron Sin 
					   1998, 24(6): 651-659.(G5 (P1 P2 F1 F2 F2:3))"),
                    br(),
                    h4("12.Gai Junyi, Wang Jiankang. Identification and estimation of a QTL model and its effects. Theor Appl
					   Genet 1998, 97(7): 1162-1168.(G6 (P1 P2 F1 F2 B1 B2))"),
                    br(),
                    h4("13.Gai Junyi, Zhang Yuan-Ming, Wang Jiankang. A joint analysis of multiple generations for QTL models
 					   extended to mixed two major genes plus polygene. Acta Agron Sin 2000, 26(4):385-391.
					   (G6 (P1 P2 F1 F2 B1 B2))"),
                    br(),
                    h4("14.Zhang Yuan-Ming, Gai Junyi, Qi Cunkou. The precision of segregating analysis of quantitative trait 
					   and its improving methods. Acta Agron Sin 2001, 27(6):787-793.
					   (G6F (P1 F1 P2 B1:2 B2:2 F2:3))"),
                    offset=1)
                    )
                    )
                    )