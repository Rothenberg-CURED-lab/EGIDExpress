#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#*#*#*#*#*#*#*#published site#*#*#*#*#*#*#*#*#
library(shiny)
library(ggplot2)
library(shinyjs)
library(shinydashboard)
library(dashboardthemes)
library(Seurat)
library(gridExtra)
options(scipen = 999)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 5 = Sherrill EoE transcriptome
# 6 = Rochman TE-7 IL-13 timecourse
# 7 = SPINK7 transcriptome
# 8 = EPC2 IL-13 
# 10 = Blanchard Esophageal Biopsy Microarray *only published*
# 11 = Julie EoG Microarray
# 12 = Blanchard Primary Esophageal Epithelial cells + IL13
# 13 = Tetsuo EDP
# 14 = Ting Single Cell 
# 15 = Ben Baruch-Morgenstern PIR-B KO Microarray 
# 19 = Rochman 3' PPIs 
# 20 = Rochman PPI Multiplex
# 23 = Rochman EPC2 +/- IL-13 & VitD #
# 24 = Tetsuo EoG RNAseq 
# 25 = Bouffi Mouse Eos IL-4 & IL-33
# 26 = Tetsuo EC RNA Seq 2021
# 28 = Julie Mouse IL13Tg esophagus
# 34 = Tetsuo EoD RNAseq
# 35 = Netali Single Cell
# 36 = Netali Mast Cell Only Single Cell
# 37 = Rochman Eso Epi in Allergic Inflammation
load("data.RData")

#load data 5
colname5 <- c("group", "gene")
group5 <- c("Control", "Control", "Control", "Control", "Control", "Control", "EoE", "EoE", "EoE", "EoE", "EoE", "EoE", "EoE", "EoE", "EoE", "EoE")
gr5 <- as.numeric(length(unique(group5)))
n5 <- as.numeric(NCOL(Genes5))
# "Control", "EoE"

#load data 6
Genes6 <- cbind(GENENAME6, Genes6)
colname6 <- c("group", "gene")
group6 <- c("Control_2h", "IL-13_2h", "Control_6h", "IL-13_6h", "Control_24h", "IL-13_24h")
gr6 <- as.numeric(length(unique(group6)))
n6 <- as.numeric(NCOL(Genes6))
# "Control_2h", "IL-13_2h", "Control_6h", "IL-13_6h", "Control_24h", "IL-13_24h"

#load data 7
colname7 <- c("group", "gene")
group7 <- c("Undif_NSC", "Undif_NSC", "Undif_NSC", "NSC", "NSC", "NSC", "SPINK7_silenced", "SPINK7_silenced", "SPINK7_silenced")
gr7 <- as.numeric(length(unique(group7)))
n7 <- as.numeric(NCOL(Genes7))
# "Undif_NSC", "NSC", "SPINK7_silenced"

#load data 8
colname8 <- c("group", "gene")
group8 <- c("Day_0", "Day_0", "Day_0", "Day_6", "Day_6", "Day_6", "Day_6_IL-13", "Day_6_IL-13", "Day_6_IL-13")
gr8 <- as.numeric(length(unique(group8)))
n8 <- as.numeric(NCOL(Genes8))
# "Day_0", "Day_6", "Day_6_IL-13"

#load data 10
Genes10 <- cbind(GENENAME10, Genes10)
colname10 <- c("group", "gene")
group10 <- c(rep("NL", 14), rep("CE", 8), rep("EoE", 18), rep("EoE+FP+R", 14), rep("EoE+FP+NR", 10))
gr10 <- as.numeric(length(unique(group10)))
n10 <- as.numeric(NCOL(Genes10))
# "NL", "CE", "EoE", "EoE+FP+R", "EoE+FP+NR"

#load data 11
Genes11 <- cbind(GENENAME11, Genes11)
colname11 <- c("group", "gene")
group11 <- c("Control", "Control", "Control", "Control", "Control", "EoG", "EoG", "EoG", "EoG", "EoG")
gr11 <- as.numeric(length(unique(group11)))
n11 <- as.numeric(NCOL(Genes11))
# "Control", "EoG"

#load data 12
Genes12 <- cbind(GENENAME12, Genes12)
colname12 <- c("group", "gene")
group12 <- c("Untreated", "Untreated", "Untreated", "IL-13", "IL-13", "IL-13")
gr12 <- as.numeric(length(unique(group12)))
n12 <- as.numeric(NCOL(Genes12))
# "Untreated", "IL-13"

#load data 13
Genes13 <- cbind(GENENAME13, Genes13)
colname13 <- c("group", "gene")
group13 <- c(rep("EoEe1", 30), rep("EoEe2", 25), rep("EoEe3", 31))
gr13 <- as.numeric(length(unique(group13)))
n13 <- as.numeric(NCOL(Genes13))
# "EoEe1", "EoEe2", "EoEe3"

#load data 14

#load data 15
colname15 <- c("group", "gene")
group15 <- c("bm_WT_noDox", "bm_WT_noDox", "bm_KO_noDox", "bm_KO_noDox", "bm_KO_noDox", "bm_WT_Dox", "bm_WT_Dox", "bm_WT_Dox", "bm_KO_Dox", "bm_KO_Dox", "bm_KO_Dox", "eso_WT_Dox", "eso_WT_Dox", "eso_WT_Dox", "eso_KO_Dox", "eso_KO_Dox", "eso_KO_Dox")
gr15 <- as.numeric(length(unique(group15)))
n15 <- as.numeric(NCOL(Genes15))
# "bm_WT_noDox", "bm_KO_noDox", "bm_WT_Dox", "bm_KO_Dox", "eso_WT_Dox", "eso_KO_Dox"

#load data 19
Genes19 <- cbind(GENENAME19, Genes19)
colname19 <- c("group", "gene")
group19 <- c( rep("Untreated", 4), rep("IL-13", 4), rep("Omeprazole", 4), rep("Esomeprazole", 4), rep("IL-13+Omeprazole", 4), rep("IL-13+Esomeprazole", 4))
gr19 <- as.numeric(length(unique(group19)))
n19 <- as.numeric(NCOL(Genes19))
#"Untreated", "IL-13", "Omeprazole", "Esomeprazole", "IL-13+Omeprazole", "IL-13+Esomeprazole"

#load data 20
Genes20 <- cbind(GENENAME20, Genes20)
colname20 <- c("group", "gene")
group20 <- c( rep("Untreated", 2), rep("Omeprazole", 2), rep("Esomeprazole", 2), rep("IL-13", 2), rep("IL-13+Omeprazole", 2), rep("IL-13+Esomeprazole", 2))
gr20 <- as.numeric(length(unique(group20)))
n20 <- as.numeric(NCOL(Genes20))
#"Untreated", "Omeprazole", "Esomeprazole", "IL-13", "IL-13+Omeprazole", "IL-13+Esomeprazole"

#load data 23
Genes23 <- cbind(GENENAME23, Genes23)
colname23 <- c("group", "gene")
group23 <- c(rep("Untreated", 3), rep("IL-13", 3), rep("VitD", 3), rep("IL-13_&_VitD", 3))
gr23 <- as.numeric(length(unique(group23)))
n23 <- as.numeric(NCOL(Genes23))

#load data 24
colname24 <- c("group", "gene")
group24 <- c(rep("Non-EoG", 12), rep("EoG", 9))
gr24 <- as.numeric(length(unique(group24)))
n24 <- as.numeric(NCOL(Genes24))
# "Non-EoG", "EoG"

#load data 25
Genes25 <- cbind(GENENAME25, Genes25)
colname25 <- c("group", "gene")
group25 <- c("untreated_1hr", "untreated_4hr", "IL4_1hr", "IL4_4hr", "IL33_1hr", "IL33_4hr")
gr25 <- as.numeric(length(unique(group25)))
n25 <- as.numeric(NCOL(Genes25))
# "untreated_1hr", "untreated_4hr", "IL4_1hr", "IL4_4hr", "IL33_1hr", "IL33_4hr"

#load data 26
Genes26 <- cbind(GENENAME26, Genes26)
colname26 <- c("group", "gene")
group26 <- c(rep("Non_EoC", 8), rep("Active_EoC", 6))
gr26 <- as.numeric(length(unique(group26)))
n26 <- as.numeric(NCOL(Genes26))
# "Active_EoC", "Non_EoC"

#load data 28
Genes28 <- cbind(GENENAME28, Genes28)
colname28 <- c("group", "gene")
group28 <- c("NoDox", "NoDox", "DOX", "DOX", "DOX")
gr28 <- as.numeric(length(unique(group28)))
n28 <- as.numeric(NCOL(Genes28))
# "NoDox", "DOX"

#load data 34
Genes34 <- cbind(GENENAME34, Genes34)
colname34 <- c("group", "gene")
group34 <- c(rep("NL", 8), rep("DE", 11), rep("EoD", 8))
gr34 <- as.numeric(length(unique(group34)))
n34 <- as.numeric(NCOL(Genes34))
# "NL", "DE", "EoD"

# Define UI 
#DEFINE SIDEBAR and BODY
sidebar <- dashboardSidebar(width = 237, 
    sidebarMenu(
      menuItem((div("RNAseq", style = "font-weight: bold")),
        menuSubItem("EoE Transcriptome", tabName = "tab5"),
        menuSubItem("EPC2 +/- SPINK7 Monolayer/ALI", tabName = "tab7"),
        menuSubItem("EPC2 +/- IL-13 Monolayer/ALI", tabName = "tab8"),
        menuSubItem("TE7 +/- IL-13", tabName = "tab6"),
        menuSubItem("EoC RNASeq", tabName = "tab26"),
        menuSubItem("EoG Transcriptome", tabName = "tab24"),
        menuSubItem("Mouse Eos", tabName = "tab25"),
        menuSubItem("Effect of PPIs on EPC2 cells", tabName = "tab19"),
        menuSubItem("EPC2 IL-13 & VitD", tabName = "tab23"),
        menuSubItem("EoD RNAseq", tabName = "tab34")),
      menuItem((div("Single-cell RNAseq", style = "font-weight: bold")),
        menuSubItem("T-cell Single Cell RNAseq", tabName = "tab14"),
        menuSubItem("Whole Biopsy Single Cell", tabName = "tab35"),
        menuSubItem("Mast Cell Single Cell", tabName = "tab36"),
        menuSubItem("Eso Epi Allergic Inflammation", tabName = "tab37")),
      menuItem((div("Microarray", style = "font-weight: bold")),
        menuSubItem("EoE Transcriptome", tabName = "tab10"),
        menuSubItem("EoG Transcriptome", tabName = "tab11"),
        menuSubItem("Eso Epi +/- IL-13", tabName = "tab12"),
        menuSubItem("M. musculus PIR-B KO", tabName = "tab15"),
        menuSubItem("Mouse IL13Tg esophagus", tabName = "tab28")),
      menuItem((div("qPCR", style = "font-weight: bold")),
        menuSubItem("EDP: EoE Endotypes", tabName = "tab13")),
      menuItem((div("Multiplex", style = "font-weight: bold")),
        menuSubItem("EPC2 +/- PPI Multiplex", tabName = "tab20")),
      menuItem(div("Abbreviations", style = "font-weight: bold"), tabName = "tababbreviation")
    )
)
body <- dashboardBody(

  shinyDashboardThemes(
    theme = "blue_gradient"
  ),
    tags$head(tags$link(rel = "shortcut icon", href = "favicon.ico")),
    tabItems(
        tabItem(tabName = "tab5",
                fluidPage(
                    
                    # Application title
                    titlePanel("EoE Transcriptome by RNA Sequencing (2014)"),
                    
                    # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                    sidebarLayout(
                        sidebarPanel(
                            selectizeInput("yvariable5",
                                           label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                        choices = NULL),
                            radioButtons("graph5", 
                                         label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                         choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                            p("Primary publication:", tags$br(),
                              "Sherrill et al, 2014", tags$br(), 
                              a(href = "https://doi.org/10.1038/gene.2014.27", "https://doi.org/10.1038/gene.2014.27")),
                            p("Experimental Design:"),
                            p("RNA sequencing was performed to quantify transcripts present in esophageal biopsy specimens obtained from 6 healthy controls (Control) and 10 patients with active EoE (EoE). ")
                        ),
                        
                        # Show the generated plot
                        mainPanel(
                            tabsetPanel(type = "tabs",
                                        tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph5"),
                                                 h5("p-value"),
                                                 verbatimTextOutput("selectedanova5")),
                                        tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata5")),
                                        tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo5")),
                                       tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData5", "Download"))
                                        
                            )
                        )
                    )
                )),
        tabItem(tabName = "tab6",
                fluidPage(
                    
                    # Application title
                    titlePanel("Kinetics of esophageal epithelial cell line (TE-7) transcriptional response to IL-13 (2014)"),
                    
                    # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                    sidebarLayout(
                        sidebarPanel(
                            selectizeInput("yvariable6",
                                           label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                           choices = NULL),
                            radioButtons("graph6", 
                                         label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                         choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                            p("Primary publication:", tags$br(),
                            "Rochman et al, 2014", tags$br(),
                            a(href = "https://doi.org/10.1038/mi.2014.109", "https://doi.org/10.1038/mi.2014.109")),
                            p("Experimental Design:"),
                            p("TE-7 cells (a human esophageal epithelial cell line) were stimulated with IL-13 (100 ng/ml) for 2h, 6h, or 24h.  Transcript levels were quantified by RNA sequencing.  Untreated cells (Control_2h, Control_4h, Control_24h); IL-13-treated cells (IL-13_2h, IL-13_6h, IL-13_24h).")
                        ),
                        # Show the generated plot
                        mainPanel(
                            tabsetPanel(type = "tabs",
                                        tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph6")),
                                        #h5("p-value"),
                                        #verbatimTextOutput("selectedanova6")),
                                        tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata6")),
                                        tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo6")),
                                       tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData6", "Download"))
                                        
                            )
                        )
                    )
                )
        ),
        tabItem(tabName = "tab7",
                fluidPage(
                    
                    # Application title
                    titlePanel("Gene expression of undifferentiated EPC2 (monolayer), differentiated EPC2 (ALI), and differentiated SPINK7-silenced EPC2 (ALI) (2018)"),
                    
                    # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                    sidebarLayout(
                        sidebarPanel(
                            selectizeInput("yvariable7",
                                           label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                           choices = NULL),
                            radioButtons("graph7", 
                                         label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                         choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                            p("Primary publication:", tags$br(),
                            "Azouz et al, 2018", tags$br(),
                            a(href = "https://doi.org/10.1126/scitranslmed.aap9736", "https://doi.org/10.1126/scitranslmed.aap9736")),
                            p("Experimental Design:"),
                            p("A telomerase-immortalized human esophageal epithelial cell line (EPC2 cells) was transfected with either non-silencing control shRNA or SPINK7 shRNA.  Cells were grown either in monolayer culture (undifferentiated) or in air liquid interface (ALI) cultures (differentiated).  Transcript levels of the following types of cells were determined:  cells transfected with non-silencing control shRNA and grown in monolayer culture (Undif_NSC); cells transfected with non-silencing control shRNA grown at the ALI (NSC); cells transfected with SPINK7 shRNA grown at the ALI (SPINK7_silenced).")
                            ),
                        # Show the generated plot
                        mainPanel(
                            tabsetPanel(type = "tabs",
                                        tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph7"),
                                                 h5("p-value"),
                                                 verbatimTextOutput("selectedanova7")),
                                        tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata7")),
                                        tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo7")),
                                       tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData7", "Download"))
                            )
                        )
                        
                    )
                )
        ),
        tabItem(tabName = "tab8",
                fluidPage(
                    
                    # Application title
                    titlePanel("Gene expression of undifferentiated EPC2 (monolayer), differentiated EPC2 (ALI), and differentiated EPC2 (ALI) +/- IL-13 (2014)"),
                    
                    # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                    sidebarLayout(
                        sidebarPanel(
                            selectizeInput("yvariable8",
                                           label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                        choices = NULL),
                            radioButtons("graph8", 
                                         label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                         choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                            p("Primary publication:", tags$br(),
                            "Sherrill et al, 2014", tags$br(),
                            a(href = "https://doi.org/10.1038/mi.2013.90", "https://doi.org/10.1038/mi.2013.90")),
                            p("Experimental Design:"),
                            p("A telomerase-immortalized human esophageal epithelial cell line (EPC2 cells) was grown submerged in monolayer culture (undifferentiated) for 8 days, or were grown submerged in monolayer culture for 8 days followed by 6 days of culture at the air liquid interface (differentiated).  In some cases, cells were treated with IL-13 (100 ng/ml) during the 6 days of ALI culture.  Transcript levels of the following conditions were assessed by RNA sequencing:  Day_0:  cells grown submerged for 8 days; Day_6:  cells grown submerged for 8 days and then grown at the ALI for 6 additional days; Day_6_IL-13:  cells grown submerged for 8 days and then grown at the ALI for 6 additional days in the presence of IL-13.")
                        ),
                        # Show the generated plot
                        mainPanel(
                            tabsetPanel(type = "tabs",
                                        tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph8"),
                                                 h5("p-value"),
                                                 verbatimTextOutput("selectedanova8")),
                                        tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata8")),
                                        tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo8")),
                                       tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData8", "Download"))
                            )
                        )
                        
                    )
                )
        ),
        tabItem(tabName = "tab10",
                fluidPage(
                  
                  # Application title
                  titlePanel("Esophageal Gene Expression (2006)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable10",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                  choices = NULL),
                      radioButtons("graph10", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                      "Blanchard et al, 2006", tags$br(),
                      a(href = "https://doi.org/10.1172/JCI26679", "https://doi.org/10.1172/JCI26679"), tags$br(),
                      "Blanchard et al, 2007", tags$br(),
                      a(href = "https://doi.org/10.1016/j.jaci.2007.10.024", "https://doi.org/10.1016/j.jaci.2007.10.024"), tags$br(),
                      "Caldwell et al, 2010", tags$br(),
                      a(href = "https://doi.org/10.1016/j.jaci.2010.01.038", "https://doi.org/10.1016/j.jaci.2010.01.038"),tags$br()
                      ),
                      p("Experimental Design:"),
                      p("RNA isolated from distal esophageal tissue was subjected to genome-wide transcript profile microarray analysis.  Experimental groups:  NL:  control; CE:  chronic esophagitis; EoE:  active EoE; EoE+FP+R:  fluticasone responder; EoE+FP+NR:  fluticasone non-responder.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph10"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova10")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata10")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo10")),
                                 tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData10", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        tabItem(tabName = "tab11",
                fluidPage(
                  
                  # Application title
                  titlePanel("Gastric Gene Expression (2007)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable11",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                  choices = NULL),
                      radioButtons("graph11", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                      "Blanchard et al, 2007", tags$br(),
                      a(href = "https://doi.org/10.1016/j.jaci.2014.07.026", "https://doi.org/10.1016/j.jaci.2014.07.026")),
                      p("Experimental Design:"),
                      p("RNA isolated from gastric antrum tissue was subjected to genome-wide transcript analysis by microarray analysis.  Experimental groups:  control and active eosinophilic gastritis (EoG).")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph11"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova11")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata11")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo11")),
                                 tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData11", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        tabItem(tabName = "tab12",
                fluidPage(
                  
                  # Application title
                  titlePanel("Transcriptional response of primary esophageal epithelial cells to IL-13 (2014)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable12",
                                  label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                  choices = NULL),
                      radioButtons("graph12", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                      "Caldwell et al, 2014", tags$br(),
                      a(href = "https://doi.org/10.1016/j.jaci.2007.10.024", "https://doi.org/10.1016/j.jaci.2007.10.024")),
                      p("Experimental Design:"),
                      p("Human primary esophageal epithelial cells from patients with EoE were cultured in the presence or absence of IL-13 (100 ng/ml) for 48 h, and the mRNA was subjected to genome-wide transcript profile microarray analysis. The data from the IL-13-stimulated cells were normalized pairwise to unstimulated controls.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph12"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova12")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata12")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo12")),
                                 tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData12", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        tabItem(tabName = "tab13",
                fluidPage(
                  
                  # Application title
                  titlePanel("EoE Diagnostic Panel: EoE Endotypes (2018)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable13",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph13", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Shoda et al, 2018", tags$br(),
                        a(href = "https://doi.org/10.1016/S2468-1253(18)30096-7", "https://doi.org/10.1016/S2468-1253(18)30096-7")),
                      p("Experimental Design:"),
                      p("RNA was isolated from distal esophageal tissue of 86 individuals with active EoE.  The RNA was subjected to EDP analysis, which quantifies 96 individual gene transcripts.  Graphs depict the -deltaCt value -(ct of GAPDH – ct of gene of interest) for each individual. Experimental groups: EoEe1: EoE endotype 1; EoEe2: EoE endotype 2; EoEe3: EoE endotype 3.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph13"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova13")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata13")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo13")),
                                 tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData13", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        tabItem(tabName = "tab14",
                fluidPage(
                  titlePanel("Single cell RNA sequencing of T-cells (2019)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable14",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph14", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Violin Plot" = 1, "Ridge Plot" = 2)),
                      checkboxInput("sep14", "Separate by EoE/Remission/Normal", value = TRUE),
                      p("Primary publication:", tags$br(),
                        "Wen et al, 2019", tags$br(),
                        a(href = "https://doi.org/10.1172/JCI125917", "https://doi.org/10.1172/JCI125917")),
                      p("Experimental Design:"),
                      p("Single-cell transcriptome profililng of esophagus biopsies from patients with active eosinophilic esophagitis (EoE, n=11), in remission (n=6), and without disease (n=5).")
                      
                      
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph14"))
                      )
                      
                    )
                  )
                )),
        tabItem(tabName = "tab15",
                fluidPage(
                  
                  # Application title
                  titlePanel("Gene expression profiling of bone marrow and esophagus from PIR-B KO mice by microarray (2016)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable15",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph15", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Ben Baruch-Morgenstern et al, 2016", tags$br(),
                        a(href = "https://doi.org/10.4049/jimmunol.1501873", "https://doi.org/10.4049/jimmunol.1501873")),
                      p("Experimental Design:"),
                      p("17 samples were analyzed in this experiment. The experiment was designed with 3 replicates for each treatment/genotype/tissue (Dox and no Dox/wildtype and knockout/bone marrow and esophagus), with the exception of the sample wildtype bone marrow no Dox where 1 replicate was dropped due to low hybridization signal. The no Dox and wildtype samples are controls for the treatment and background of the mice. Experimental groups: bm_WT_noDox: bone marrow, wildtype, no Dox; bm_KO_noDox: bone marrow, PirB knockout, no Dox; bm_WT_Dox: bone marrow, wildtype, Dox; bm_KO_Dox: bone marrow, PirB knockout, Dox; eso_WT_Dox: esophagus, wildtype, Dox; eso_KO_Dox: esophagus, PirB knockout, Dox.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph15"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova15")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata15")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo15")),
                                 tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData15", "Download"))
                      )
                    )
                    
                  )
                )),
        #UI 26
        tabItem(tabName = "tab26",
                fluidPage(
                  titlePanel("EoC RNASeq (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable26",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph26", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Shoda et al., 2022", tags$br(),
                        a(href = "https://doi.org/10.1053/j.gastro.2022.01.022", "https://doi.org/10.1053/j.gastro.2022.01.022")),
                      p("Experimental Design:"),
                      p("RNA sequencing was performed to quantify transcripts present in colon biopsy specimens obtained from 8 control individuals (Non_EoC) and 6 patients with active EoC (Active_EoC).")
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph26"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova26")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata26")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData26", "Download"))
                      )
                    )
                  )
                )),
        tabItem(tabName = "tab24",
                fluidPage(
                  
                  # Application title
                  titlePanel("EoG Transcriptome by RNA Sequencing (2020)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable24",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph24", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Shoda et al, 2020", tags$br(),
                        a(href = "https://doi.org/10.1016/j.jaci.2019.11.007", "https://doi.org/10.1016/j.jaci.2019.11.007")),
                      p("Experimental Design:"),
                      p("RNA sequencing was performed to quantify transcripts present in gastric biopsy specimens obtained from 12 healthy controls (Non-EoG) and 9 patients with EG (EoG). ")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph24"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova24")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata24")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo24")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData24", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        tabItem(tabName = "tab25",
                fluidPage(
                  
                  # Application title
                  titlePanel("Mouse Eosinophils +/- IL4 or IL33 (2013)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable25",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph25", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Bouffi et al, 2013", tags$br(),
                        a(href = "https://doi.org/10.4049/jimmunol.1301465", "https://doi.org/10.4049/jimmunol.1301465")),
                      p("Experimental Design:"),
                      p("RNA sequencing was performed to quantify transcripts present in LDBM Mouse Eosinophils untreated, treated with IL-4, or treated with IL-33 at 1 hour and 4 hours.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph25")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata25")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo25")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData25", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        #UI 19
        tabItem(tabName = "tab19",
                fluidPage(
                  titlePanel("Effect of PPIs omeprazole and esomeprazole on gene expression in EPC2 cells (2020)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable19",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph19", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Rochman et al, 2020", tags$br(),
                        a(href = "https://doi.org/10.1016/j.jaci.2020.09.039", "https://doi.org/10.1016/j.jaci.2020.09.039")),
                      p("Experimental Design:"),
                      p("For the monolayer cultures, EPC2 cells were seeded at a density of 2.5x10^5 cells/well in KSFM in a 24-well plate. The next day, the medium was replenished. After 24 hours, stimulants were added in a total of 350 to 500mL of fresh medium for 24 hours. IL-13 was added to a final concentration of 100 ng/mL unless otherwise indicated. Omeprazole and esomeprazole were dissolved in DMSO to the stock concentration of 100 mM, aliquoted into 10-mL amounts, and stored at –80°C. PPIs were thawed once and used at a working concentration of 100mM; cells were pretreated with the PPIs for 1 hour before adding IL-13 to the medium. GNF-351 was dissolved in DMSO to a 20 mM stock concentration. The stock was aliquoted in 10-mL doses and stored at –80°C.  Cells were pretreated with GNF-351 at a final concentration of 2mM for 30 minutes before adding PPIs."),
                      p("RNA sequencing was performed with high-quality RNA (RNA integrity number > 8) by using the QuantSeq 3' mRNA-Seq Library Prep Kit FWD for Illumina (Lexogen, Vienna, Austria, catalog no 015.96)."),
                      p("Libraries were diluted to final concentrations of 5 nM and sequenced on a HiSeq 4000 Illumina sequencing machine at the Genomics and Cell Characterization Core Facility at the University of Oregon with 100- to 150-bp-length single reads.")
                      
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph19"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova19")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata19")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo19")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData19", "Download"))
                      )
                      
                    )
                  )
                )),
        #UI 20
        tabItem(tabName = "tab20",
                fluidPage(
                  titlePanel("EPC2 +/- PPI Multiplex (2020)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable20",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph20", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Rochman et al, 2020", tags$br(),
                        a(href = "https://doi.org/10.1016/j.jaci.2020.09.039", "https://doi.org/10.1016/j.jaci.2020.09.039")),
                      p("Experimental Design:"),
                      p("For the monolayer cultures, EPC2 cells were seeded at a density of 2.5x10^5 cells/well in KSFM in a 24-well plate. The next day, the medium was replenished. After 24 hours, stimulants were added in a total of 350 to 500mL of fresh medium for 24 hours. IL-13 was added to a final concentration of 100 ng/mL unless otherwise indicated. Omeprazole and esomeprazole were dissolved in DMSO to the stock concentration of 100 mM, aliquoted into 10-mL amounts, and stored at –80°C. PPIs were thawed once and used at a working concentration of 100mM; cells were pretreated with the PPIs for 1 hour before adding IL-13 to the medium. GNF-351 was dissolved in DMSO to a 20 mM stock concentration. The stock was aliquoted in 10-mL doses and stored at –80°C.  Cells were pretreated with GNF-351 at a final concentration of 2mM for 30 minutes before adding PPIs."),
                      p("Multiplex analysis was performed by Eve Technologies (Calgary, Canada) using the Human Cytokine Array/Chemokine Array 65-Plex Panel (HD65).")
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph20"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova20")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata20")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo20")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData20", "Download"))
                      )
                      
                    )
                  )
                )),
        #UI 23
        tabItem(tabName = "tab23",
                fluidPage(
                  titlePanel("Effect of IL-13 and Vitamin D on gene expression in EPC2 cells (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable23",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph23", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Bar graph:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Brusilovsky et al, 2022", tags$br(),
                        a(href = "https://doi.org/10.1136/gutjnl-2022-327276", "https://doi.org/10.1136/gutjnl-2022-327276")),
                      p("Experimental Design:"),
                      p("The contribution of vitamin D (VD) deficiency to the pathogenesis of allergic diseases remains elusive. We aimed to define the impact of VD on esophageal allergic inflammation. To this end, we examined the potential clinical relevance of the esophageal epithelial cell (EPC2) transcriptional responses to IL- 13 or VD. Cells were grown in the monolayer culture at a super-confluent state (1×106 cells/cm2). Cells were treated with recombinant human IL- 13 (100 ng/ mL) and/or VD (1,25- dihydroxycholecalciferol; 100 nM) and subjected to bulk RNA sequencing."),
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph23"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova23")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata23")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo23")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData23", "Download"))
                      )
                      
                    )
                  )
                )),
        #UI 28
        tabItem(tabName = "tab28",
                fluidPage(
                  titlePanel("Mouse IL13Tg Esophagus (2010)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable28",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph28", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Primary publication:", tags$br(),
                        "Zuo et al. 2010, JI", tags$br(),
                        a(href = "https://www.jimmunol.org/content/185/1/660", "https://www.jimmunol.org/content/185/1/660")),
                      p("Experimental Design:"),
                      p("Bi-transgenic mice (CC10-iIL-13) were generated in which IL-13 was expressed in a lung-specific manner that allowed for external regulation of the transgene expression, as previously described(13).", tags$br(), tags$br(), 
                        "RNA from the esophagus or lung of IL-13 inducible transgenic mice obtained after 30 days of DOX or NO-DOX treatment was subjected to gene chip analysis using MOE 430 2.0 chips as previously reported (15). Gene transcript levels were determined using algorithms in the Microarray Analysis Suite and GeneSpring software (Silicon Genetics). For comparison with human EoE, genes were translated to human homologues represented on the human U133 chip using GeneSpring software and compared to microarray results previously described(8, 15). Welch T test and fold change cut-off were performed. The correlation between IL-13 expression and numbers of dysregulated genes were also analyzed using GeneSpring software.", tags$br(), tags$br(), 
                        "2 NO-DOX samples and 3 DOX samples", tags$br(), tags$br(),
                        "MOE 430 2.0 chips - Affymetrix", tags$br(),
                        "Composition of probe sets can be identified/analyzed on NetAffx", tags$br(),
                        a(href = "https://www.affymetrix.com/site/login/login.affx?toURL=/analysis/netaffx/index.affx", "https://www.affymetrix.com/site/login/login.affx?toURL=/analysis/netaffx/index.affx")
                      ),
                      p("References:", tags$br(),
                        "8 -", a(href = "https://pubmed.ncbi.nlm.nih.gov/18073124/", "https://pubmed.ncbi.nlm.nih.gov/18073124/"), tags$br(),
                        "13 -", a(href = "https://pubmed.ncbi.nlm.nih.gov/14757645/", "https://pubmed.ncbi.nlm.nih.gov/14757645/"), tags$br(),
                        "15 -", a(href = "https://pubmed.ncbi.nlm.nih.gov/16453027/", "https://pubmed.ncbi.nlm.nih.gov/16453027/")
                      )),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph28"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova28")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata28")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo28")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData28", "Download"))
                      )
                      
                    )
                  )
                )),
        #UI 34
        tabItem(tabName = "tab34",
                fluidPage(
                  
                  # Application title
                  titlePanel("EoD RNAseq (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable34",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph34", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Bar graph:  group average +/- SEM" = 1, "Box and whisker plot" = 2, "Dot plot:  individual samples" = 3)),
                      p("Experimental Design:"),
                      p("RNA sequencing was performed to quantify transcripts present in duodenal biopsy specimens obtained from non-EoD controls (NL), patients with duodenal eosinophilia (DE), or patients with eosinophilic duodenitis (EoD).  Graphs were depicted by the TPM for each individual.")
                    ),
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph34"),
                                           h5("p-value"),
                                           verbatimTextOutput("selectedanova34")),
                                  tabPanel(div(icon("fas fa-table", style = "color:black"), "Expression data"), tableOutput("selecteddata34")),
                                  tabPanel(div(icon("fas fa-clipboard", style = "color:black"), "Sample Information"), tableOutput("sampleinfo34")),
                                  tabPanel(div(icon("fas fa-download", style = "color:black"), "Download data"), downloadButton("downloadData34", "Download"))
                      )
                    )
                    
                  )
                )
        ),
        #UI 35
        tabItem(tabName = "tab35",
                fluidPage(
                  titlePanel("Single Cell Whole Biopsy RNAseq (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable35",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph35", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Feature Plot" = 1, "Violin Plot" = 2, "Ridge Plot" = 3)),
                      p("Primary publication:", tags$br(),
                        "Ben-Baruch Morgenstern et al. 2022, JACI", tags$br(),
                        a(href = "https://doi.org/10.1016/j.jaci.2022.02.025", "https://doi.org/10.1016/j.jaci.2022.02.025")),
                      p("Rochman et al. 2022, JCI", tags$br(),
                        a(href = "https://doi.org/10.1172/jci.insight.159093", "https://doi.org/10.1172/jci.insight.159093")),
                      p("Experimental Design:"),
                      p("Esophageal biopsies were obtained from patients with active EoE (n = 5), patients with EoE in histologic remission (n = 3), and control individuals with histologically normal esophageal biopsies and no history of esophageal disease (n = 2) and were subjected to 10X scRNA-sequencing.")
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph35"))
                      )
                      
                    )
                  )
                )),
        #UI 36
        tabItem(tabName = "tab36",
                fluidPage(
                  titlePanel("Mast Cell Only Single Cell RNAseq (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable36",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph36", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Feature Plot" = 1, "Violin Plot" = 2, "Ridge Plot" = 3)),
                      p("Primary publication:", tags$br(),
                        "Ben-Baruch Morgenstern et al. 2022, JACI", tags$br(),
                        a(href = "https://doi.org/10.1016/j.jaci.2022.02.025", "https://doi.org/10.1016/j.jaci.2022.02.025")),
                      p("Experimental Design:"),
                      p("Esophageal biopsies were obtained from patients with active EoE (n = 5), patients with EoE in histologic remission (n = 3), and control individuals with histologically normal esophageal biopsies and no history of esophageal disease (n = 2) and were subjected to 10X scRNA-sequencing.  This dataset includes mast cells only.")
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph36"))
                      )
                      
                    )
                  )
                )),
        #UI 37
        tabItem(tabName = "tab37",
                fluidPage(
                  titlePanel("Esophageal Epithelium in Allergic Inflammation (2022)"),
                  
                  # Sidebar with a set of radio buttons to choose comparison and drop down box to choose the y axis variable 
                  sidebarLayout(
                    sidebarPanel(
                      selectizeInput("yvariable37",
                                     label=div(icon("fas fa-filter", style = "color: black"), "Gene:"),
                                     choices = NULL),
                      radioButtons("graph37", 
                                   label = div(icon("chart-line", style = "color:black"), "Type of graph:"),
                                   choices = list("Feature Plot" = 1, "Violin Plot" = 2, "Ridge Plot" = 3)),
                      p("Primary publication:", tags$br(),
                        "Rochman et al. 2022, JCI", tags$br(),
                        a(href = "https://doi.org/10.1172/jci.insight.159093", "https://doi.org/10.1172/jci.insight.159093")),
                      p("Experimental Design:"),
                      p("Esophageal biopsies were obtained from patients with active EoE (n = 5), patients with EoE in histologic remission (n = 3), and control individuals with histologically normal esophageal biopsies and no history of esophageal disease (n = 2) and were subjected to 10X scRNA-sequencing.  This dataset includes epithelial cells only.")
                    ),
                    
                    # Show the generated plot
                    mainPanel(
                      tabsetPanel(type = "tabs",
                                  tabPanel(div(icon("chart-bar", style = "color:black"), "Graph"), plotOutput("selectedgraph37"))
                      )
                      
                    )
                  )
                )),
        #UI Abbreviations
        tabItem(tabName = "tababbreviation",
                fluidPage(
                  titlePanel("Abbreviation Glossary"),
                    mainPanel(
                      tableOutput("abbreviations"))
                      )
                )
#this goes after last tab UI        
))
ui <- dashboardPage(
    dashboardHeader(title = "EGID Express", titleWidth = 237,                
                    tags$li(a(shiny::icon("fas fa-database", style = "color: black"),"Rothenberg CURED Lab Members", style = "color: black", href = 'https://egidexpress.research.cchmc.org/lab',
                            title = "Rothenberg Lab Members"),
                            class = "dropdown")),
#                    tags$li(a(href = 'https://www.cincinnatichildrens.org/research/divisions/a/allergy-immunology/labs/rothenberg',
#                    tags$img(src = 'logo.png',
#                            title = "Cincinnati Children's", height = "40px"),
#                            style = "padding-top:5px; padding-bottom:5px;"),
#                            class = "dropdown")),
    sidebar,
    body
)                
server <- function(session, input, output) {
  updateSelectizeInput(session, "yvariable5", choices = c(GENENAME5), server = TRUE)
  updateSelectizeInput(session, "yvariable6", choices = c(GENENAME6), server = TRUE)
  updateSelectizeInput(session, "yvariable7", choices = c(GENENAME7), server = TRUE)
  updateSelectizeInput(session, "yvariable8", choices = c(GENENAME8), server = TRUE)
  updateSelectizeInput(session, "yvariable10", choices = c(GENENAME10), server = TRUE)
  updateSelectizeInput(session, "yvariable11", choices = c(GENENAME11), server = TRUE)
  updateSelectizeInput(session, "yvariable12", choices = c(GENENAME12), server = TRUE)
  updateSelectizeInput(session, "yvariable13", choices = c(GENENAME13), server = TRUE)
  updateSelectizeInput(session, "yvariable14", choices = c(GENENAME14), server = TRUE)
  updateSelectizeInput(session, "yvariable15", choices = c(GENENAME15), server = TRUE)
  updateSelectizeInput(session, "yvariable19", choices = c(GENENAME19), server = TRUE)
  updateSelectizeInput(session, "yvariable20", choices = c(GENENAME20), server = TRUE)
  updateSelectizeInput(session, "yvariable23", choices = c(GENENAME23), server = TRUE)
  updateSelectizeInput(session, "yvariable24", choices = c(GENENAME24), server = TRUE)
  updateSelectizeInput(session, "yvariable25", choices = c(GENENAME25), server = TRUE)
  updateSelectizeInput(session, "yvariable26", choices = c(GENENAME26), server = TRUE) 
  updateSelectizeInput(session, "yvariable28", choices = c(GENENAME28), server = TRUE)
  updateSelectizeInput(session, "yvariable34", choices = c(GENENAME34), server = TRUE)
  updateSelectizeInput(session, "yvariable35", choices = c(GENENAME35), server = TRUE)
  updateSelectizeInput(session, "yvariable36", choices = c(GENENAME36), server = TRUE)
  updateSelectizeInput(session, "yvariable37", choices = c(GENENAME37), server = TRUE)
  
  #server for 5
    y5 <- reactive({
        cbind(group5, as.data.frame(t(Genes5[input$yvariable5, 2:17])))
    })
    toplot5 <- reactive({
        dat5 <- y5()
        colnames(dat5) <- c("group", "gene")
        dat5$group <- factor(dat5$group, levels = c("Control", "EoE"))
        dat5
    })
    g5 <- reactive({
        graph5 <- ggplot(data = toplot5()) +
            ggtitle(input$yvariable5) +
            theme(aspect.ratio = as.numeric(2*gr5/n5), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
        graph5
    })
    output$selectedgraph5 <- renderPlot({
        graph5 <- g5()
        if (input$graph5 == 1){graph5 <- g5() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
        if (input$graph5 == 2){graph5 <- g5() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
        if (input$graph5 == 3){graph5 <- g5() + aes(x=row.names(toplot5()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot5())) + scale_fill_manual(values = cbPalette)}
        graph5 + xlab(NULL) + ylab("Expression, RPKM")
    })
    output$selectedanova5 <- renderPrint({
        x5 <- data.frame(toplot5())
        df5 <- round(as.data.frame(pairwise.t.test(x5$gene, x5$group, p.adjust.method = "fdr")$p.value), digits = 3)
        if (df5[1,1] == "0") {print("<0.001")}
        else {print(df5[1,1])}
    })
    output$selecteddata5 <- renderTable({
        d5 <- as.data.frame(Genes5[input$yvariable5, 2:17])
        d5 <- rbind(colnames(d5), d5)
        row.names(d5) <- c("Sample", "RPKM")
        d5 <- as.data.frame(t(d5))
    })
    output$downloadData5 <- downloadHandler(
        filename = paste("EoEtranscriptome.csv", sep = ""),
        content = function(file){
            write.csv(data_raw5, file, row.names = TRUE)}
    )
    output$sampleinfo5 <- renderTable({
      sample_info5
    })
    
    #server for 6
    y6 <- reactive({
        cbind(group6, as.data.frame(t(Genes6[input$yvariable6, 2:7])))
    })
    toplot6 <- reactive({
        dat6 <- y6()
        colnames(dat6) <- c("group", "gene")
        dat6$group <- factor(dat6$group, levels = c("Control_2h", "Control_6h", "Control_24h", "IL-13_2h", "IL-13_6h", "IL-13_24h"))
        dat6
    })
    g6 <- reactive({
        graph6 <- ggplot(data = toplot6()) +
            ggtitle(input$yvariable6) +
            theme(aspect.ratio = as.numeric(2*gr6/n6), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
        graph6
    })
    output$selectedgraph6 <- renderPlot({
        graph6 <- g6()
        if (input$graph6 == 1){graph6 <- g6() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
        if (input$graph6 == 2){graph6 <- g6() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
        if (input$graph6 == 3){graph6 <- g6() + aes(x=row.names(toplot6()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot6())) + scale_fill_manual(values = cbPalette)}
        graph6 + xlab(NULL) + ylab("Expression, RPKM")
    })
    #output$selectedanova6 <- renderPrint({
    #  req(credentials()$user_auth)
    #  x6 <- data.frame(toplot6())
    #  df6 <- round(as.data.frame(pairwise.t.test(x6$gene, x6$group, p.adjust.method = "fdr")$p.value), digits = 3)
    #  df6 <- sapply(df6, as.character)
    #  df6[is.na(df6)] <- " "
    #  df6[df6 == "0"] <- "<0.001"
    #  as.data.frame(df6)
    #      })
    output$selecteddata6 <- renderTable({
        d6 <- as.data.frame(Genes6[input$yvariable6, 2:7])
        d6 <- rbind(colnames(d6), d6)
        row.names(d6) <- c("Sample", "RPKM")
        d6 <- as.data.frame(t(d6))
    }) 
    output$downloadData6 <- downloadHandler(
        filename = paste("TE-7_IL-13_response_timecourse.csv", sep = ""),
        content = function(file){
            write.csv(data_raw6, file, row.names = TRUE)}
    )
    output$sampleinfo6 <- renderTable({
      sample_info6
    })
    
    #server for 7
    y7 <- reactive({
        cbind(group7, as.data.frame(t(Genes7[input$yvariable7, 2:10])))
    })
    toplot7 <- reactive({
        dat7 <- y7()
        colnames(dat7) <- c("group", "gene")
        dat7$group <- factor(dat7$group, levels = c("Undif_NSC", "NSC", "SPINK7_silenced"))
        dat7
    })
    g7 <- reactive({
        graph7 <- ggplot(data = toplot7()) +
            ggtitle(input$yvariable7) +
            theme(aspect.ratio = as.numeric(2*gr7/n7), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
        graph7
    })
    output$selectedgraph7 <- renderPlot({
        graph7 <- g7()
        if (input$graph7 == 1){graph7 <- g7() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
        if (input$graph7 == 2){graph7 <- g7() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
        if (input$graph7 == 3){graph7 <- g7() + aes(x=row.names(toplot7()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot7())) + scale_fill_manual(values = cbPalette)}
        graph7 + xlab(NULL) + ylab("Expression, RPKM")
    })
    output$selectedanova7 <- renderPrint({
        x7 <- data.frame(toplot7())
        df7 <- round(as.data.frame(pairwise.t.test(x7$gene, x7$group, p.adjust.method = "fdr")$p.value), digits = 3)
        df7 <- sapply(df7, as.character)
        df7[is.na(df7)] <- " "
        df7[df7 == "0"] <- "<0.001"
        row.names(df7) <- c("NSC", "SPINK7_silenced")
        colnames(df7) <- c("Undif_NSC", "NSC")
        as.data.frame(df7)
    })
    output$selecteddata7 <- renderTable({
        d7 <- as.data.frame(Genes7[input$yvariable7, 2:10])
        d7 <- rbind(colnames(d7), d7)
        row.names(d7) <- c("Sample", "RPKM")
        d7 <- as.data.frame(t(d7))
    })
    output$downloadData7 <- downloadHandler(
        filename = paste("SPINK7_transcriptome.csv", sep = ""),
        content = function(file){
            write.csv(data_raw7, file, row.names = TRUE)}
    )
    output$sampleinfo7 <- renderTable({
      sample_info7
    })
    
    #server for 8
    y8 <- reactive({
        cbind(group8, as.data.frame(t(Genes8[input$yvariable8, 2:10])))
    })
    toplot8 <- reactive({
        dat8 <- y8()
        colnames(dat8) <- c("group", "gene")
        dat8$group <- factor(dat8$group, levels = c("Day_0", "Day_6", "Day_6_IL-13"))
        dat8
    })
    g8 <- reactive({
        graph8 <- ggplot(data = toplot8()) +
            ggtitle(input$yvariable8) +
            theme(aspect.ratio = as.numeric(2*gr8/n8), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
        graph8
    })
    output$selectedgraph8 <- renderPlot({
        graph8 <- g8()
        if (input$graph8 == 1){graph8 <- g8() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
        if (input$graph8 == 2){graph8 <- g8() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
        if (input$graph8 == 3){graph8 <- g8() + aes(x=row.names(toplot8()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot8())) + scale_fill_manual(values = cbPalette)}
        graph8 + xlab(NULL) + ylab("Expression, RPKM")
    })
    output$selectedanova8 <- renderPrint({
        x8 <- data.frame(toplot8())
        df8 <- round(as.data.frame(pairwise.t.test(x8$gene, x8$group, p.adjust.method = "fdr")$p.value), digits = 3)
        df8 <- sapply(df8, as.character)
        df8[is.na(df8)] <- " "
        df8[df8 == "0"] <- "<0.001"
        row.names(df8) <- c("Day_6", "Day_6_IL-13")
        colnames(df8) <- c("Day_0", "Day_6")
        as.data.frame(df8)
    })
    output$selecteddata8 <- renderTable({
        d8 <- as.data.frame(Genes8[input$yvariable8, 2:10])
        d8 <- rbind(colnames(d8), d8)
        row.names(d8) <- c("Sample", "RPKM")
        d8 <- as.data.frame(t(d8))
    })
    output$downloadData8 <- downloadHandler(
        filename = paste("EPC2_IL13.csv", sep = ""),
        content = function(file){
            write.csv(data_raw8, file, row.names = TRUE)}
    )
    output$sampleinfo8 <- renderTable({
      sample_info8
    })
    
    #server for 10
    y10 <- reactive({
      cbind(group10, as.data.frame(t(Genes10[input$yvariable10, 2:65])))
    })
    toplot10 <- reactive({
      dat10 <- y10()
      colnames(dat10) <- c("group", "gene")
      dat10$group <- factor(dat10$group, levels = c("NL", "CE", "EoE", "EoE+FP+R", "EoE+FP+NR"))
      dat10
    })
    g10 <- reactive({
      graph10 <- ggplot(data = toplot10()) +
        ggtitle(input$yvariable10) +
        theme(aspect.ratio = as.numeric(2*gr10/n10), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph10
    })
    output$selectedgraph10 <- renderPlot({
      graph10 <- g10()
      if (input$graph10 == 1){graph10 <- g10() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph10 == 2){graph10 <- g10() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph10 == 3){graph10 <- g10() + aes(x=row.names(toplot10()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot10())) + scale_fill_manual(values = cbPalette) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6))
}
      graph10 + xlab(NULL) + ylab("Expression")
    })
    output$selectedanova10 <- renderPrint({
      x10 <- data.frame(toplot10())
      df10 <- round(as.data.frame(pairwise.t.test(x10$gene, x10$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df10 <- sapply(df10, as.character)
      df10[is.na(df10)] <- " "
      df10[df10 == "0"] <- "<0.001"
      row.names(df10) <- c("CE", "EoE", "EoE+FP+R", "EoE+FP+NR")
      colnames(df10) <- c("NL", "CE", "EoE", "EoE+FP+R")
      as.data.frame(df10)
    })
    output$selecteddata10 <- renderTable({
      d10 <- as.data.frame(Genes10[input$yvariable10, 2:65])
      d10 <- rbind(colnames(d10), d10)
      row.names(d10) <- c("Sample", "Expression")
      d10 <- as.data.frame(t(d10))
    })
    output$downloadData10 <- downloadHandler(
      filename = paste("Esophageal_microarray.csv", sep = ""),
      content = function(file){
        write.csv(data_raw10, file, row.names = TRUE)}
    )
    output$sampleinfo10 <- renderTable({
      sample_info10
    })

    #server for 11
    y11 <- reactive({
      cbind(group11, as.data.frame(t(Genes11[input$yvariable11, 2:11])))
    })
    toplot11 <- reactive({
      dat11 <- y11()
      colnames(dat11) <- c("group", "gene")
      dat11$group <- factor(dat11$group, levels = c("Control", "EoG"))
      dat11
    })
    g11 <- reactive({
      graph11 <- ggplot(data = toplot11()) +
        ggtitle(input$yvariable11) +
        theme(aspect.ratio = as.numeric(2*gr11/n11), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph11
    })
    output$selectedgraph11 <- renderPlot({
      graph11 <- g11()
      if (input$graph11 == 1){graph11 <- g11() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph11 == 2){graph11 <- g11() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph11 == 3){graph11 <- g11() + aes(x=row.names(toplot11()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot11())) + scale_fill_manual(values = cbPalette)}
      graph11 + xlab(NULL) + ylab("Expression")
    })
   output$selectedanova11 <- renderPrint({
      x11 <- data.frame(toplot11())
      df11 <- round(as.data.frame(pairwise.t.test(x11$gene, x11$group, p.adjust.method = "fdr")$p.value), digits = 3)
      if (df11[1,1] == "0") {print("<0.001")}
      else {print(df11[1,1])}
    })
    output$selecteddata11 <- renderTable({
      d11 <- as.data.frame(Genes11[input$yvariable11, 2:11])
      d11 <- rbind(colnames(d11), d11)
      row.names(d11) <- c("Sample", "Expression")
      d11 <- as.data.frame(t(d11))
    })
    output$downloadData11 <- downloadHandler(
      filename = paste("EG_microarray.csv", sep = ""),
      content = function(file){
        write.csv(data_raw11, file, row.names = TRUE)}
    )
    output$sampleinfo11 <- renderTable({
      sample_info11
    })
    
    #server for 12
    y12 <- reactive({
      cbind(group12, as.data.frame(t(Genes12[input$yvariable12, 2:7])))
    })
    toplot12 <- reactive({
      dat12 <- y12()
      colnames(dat12) <- c("group", "gene")
      dat12$group <- factor(dat12$group, levels = c("Untreated", "IL-13"))
      dat12
    })
    g12 <- reactive({
      graph12 <- ggplot(data = toplot12()) +
        ggtitle(input$yvariable12) +
        theme(aspect.ratio = as.numeric(2*gr12/n12), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph12
    })
    output$selectedgraph12 <- renderPlot({
      graph12 <- g12()
      if (input$graph12 == 1){graph12 <- g12() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph12 == 2){graph12 <- g12() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph12 == 3){graph12 <- g12() + aes(x=row.names(toplot12()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot12())) + scale_fill_manual(values = cbPalette)}
      graph12 + xlab(NULL) + ylab("Expression")
    })
    output$selectedanova12 <- renderPrint({
      x12 <- data.frame(toplot12())
      df12 <- round(as.data.frame(pairwise.t.test(x12$gene, x12$group, p.adjust.method = "fdr")$p.value), digits = 3)
      if (df12[1,1] == "0") {print("<0.001")}
      else {print(df12[1,1])}
    })
    output$selecteddata12 <- renderTable({
      d12 <- as.data.frame(Genes12[input$yvariable12, 2:7])
      d12 <- rbind(colnames(d12), d12)
      row.names(d12) <- c("Sample", "Expression")
      d12 <- as.data.frame(t(d12))
    })
    output$downloadData12 <- downloadHandler(
      filename = paste("primary_epithelial_microarray.csv", sep = ""),
      content = function(file){
        write.csv(data_raw12, file, row.names = TRUE)}
    )
    output$sampleinfo12 <- renderTable({
      sample_info12
    })
    
    #server for 13
    y13 <- reactive({
      cbind(group13, as.data.frame(t(Genes13[input$yvariable13, 2:87])))
    })
    toplot13 <- reactive({
      dat13 <- y13()
      colnames(dat13) <- c("group", "gene")
      dat13$group <- factor(dat13$group, levels = c("EoEe1", "EoEe2", "EoEe3"))
      dat13
    })
    g13 <- reactive({
      graph13 <- ggplot(data = toplot13()) +
        ggtitle(input$yvariable13) +
        theme(aspect.ratio = as.numeric(4*gr13/n13), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph13
    })
    output$selectedgraph13 <- renderPlot({
      graph13 <- g13()
      if (input$graph13 == 1){graph13 <- g13() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph13 == 2){graph13 <- g13() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph13 == 3){graph13 <- g13() + aes(x=row.names(toplot13()), y=gene, color=group) + geom_point(stat="identity") + scale_x_discrete(limits=row.names(toplot13())) + scale_color_manual(values = cbPalette) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6))}
      graph13 + xlab(NULL) + ylab("Expression")
    })
    output$selectedanova13 <- renderPrint({
      x13 <- data.frame(toplot13())
      df13 <- round(as.data.frame(pairwise.t.test(x13$gene, x13$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df13 <- sapply(df13, as.character)
      df13[is.na(df13)] <- " "
      df13[df13 == "0"] <- "<0.001"
      row.names(df13) <- c("EoEe2", "EoEe3")
      colnames(df13) <- c("EoEe1", "EoEe2")
      as.data.frame(df13)
    })
    output$selecteddata13 <- renderTable({
      d13 <- as.data.frame(Genes13[input$yvariable13, 2:87])
      d13 <- rbind(colnames(d13), d13)
      row.names(d13) <- c("Sample", "Expression")
      d13 <- as.data.frame(t(d13))
    })
    output$downloadData13 <- downloadHandler(
      filename = paste("EDP.csv", sep = ""),
      content = function(file){
        write.csv(data_raw13, file, row.names = TRUE)}
    )
    output$sampleinfo13 <- renderTable({
      sample_info13
    })
    
    #server for 14
    output$selectedgraph14 <- renderPlot({
      req(input$yvariable14)
      if (input$sep14 == TRUE){
        if (input$graph14 == 1){graph14 <- grid.arrange(VlnPlot(EoE14, features = input$yvariable14, group.by = "clusters") + ggtitle("EoE") + scale_x_discrete(name = ""), VlnPlot(Remission14, features = input$yvariable14, group.by = "clusters") + ggtitle("Remission") + scale_x_discrete(name = "") + theme(legend.position = 'none'), VlnPlot(Normal14, features = input$yvariable14, group.by = "clusters", cols = my_color_palette) + ggtitle("Normal") + scale_x_discrete(name = "") + theme(legend.position = 'none'), nrow = 3, left = input$yvariable14)}
        else if (input$graph14 == 2){graph14 <- grid.arrange(RidgePlot(EoE14, features = input$yvariable14, group.by = "clusters") + ggtitle("EoE"), RidgePlot(Remission14, features = input$yvariable14, group.by = "clusters") + ggtitle("Remission"), RidgePlot(Normal14, features = input$yvariable14, group.by = "clusters") + ggtitle("Normal"), nrow = 3, top = input$yvariable14)}
      }
      else if (input$sep14 == FALSE) {
        if (input$graph14 == 1){graph14 <- VlnPlot(data14, features = input$yvariable14, group.by = "clusters")}
        else if (input$graph14 == 2){graph14 <- RidgePlot(data14, features = input$yvariable14, group.by = "clusters")}
      }
      graph14
    })

    #server for 15
    y15 <- reactive({
      cbind(group15, as.data.frame(t(Genes15[input$yvariable15, 2:18])))
    })
    toplot15 <- reactive({
      dat15 <- y15()
      colnames(dat15) <- c("group", "gene")
      dat15$group <- factor(dat15$group, levels = c("bm_WT_noDox", "bm_KO_noDox", "bm_WT_Dox", "bm_KO_Dox", "eso_WT_Dox", "eso_KO_Dox"))
      dat15
    })
    g15 <- reactive({
      graph15 <- ggplot(data = toplot15()) +
        ggtitle(input$yvariable15) +
        theme(aspect.ratio = as.numeric(4*gr15/n15), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph15
    })
    output$selectedgraph15 <- renderPlot({
      graph15 <- g15()
      if (input$graph15 == 1){graph15 <- g15() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph15 == 2){graph15 <- g15() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph15 == 3){graph15 <- g15() + aes(x=row.names(toplot15()), y=gene, color=group) + geom_point(stat="identity") + scale_x_discrete(limits=row.names(toplot15())) + scale_color_manual(values = cbPalette) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6))}
      graph15 + xlab(NULL) + ylab("Expression")
    })
    output$selectedanova15 <- renderPrint({
      x15 <- data.frame(toplot15())
      df15 <- round(as.data.frame(pairwise.t.test(x15$gene, x15$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df15 <- sapply(df15, as.character)
      df15[is.na(df15)] <- " "
      df15[df15 == "0"] <- "<0.001"
      row.names(df15) <- c("bm_KO_noDox", "bm_WT_Dox", "bm_KO_Dox", "eso_WT_Dox", "eso_KO_Dox")
      colnames(df15) <- c("bm_WT_noDox", "bm_KO_noDox", "bm_WT_Dox", "bm_KO_Dox", "eso_WT_Dox")
      as.data.frame(df15)
    })
    output$selecteddata15 <- renderTable({
      d15 <- as.data.frame(Genes15[input$yvariable15, 2:18])
      d15 <- rbind(colnames(d15), d15)
      row.names(d15) <- c("Sample", "Expression")
      d15 <- as.data.frame(t(d15))
    })
    output$downloadData15 <- downloadHandler(
      filename = paste("PIR_B.csv", sep = ""),
      content = function(file){
        write.csv(data_raw15, file, row.names = TRUE)}
    )  
    output$sampleinfo15 <- renderTable({
      sample_info15
    })
    
    #server for 19
    y19 <- reactive({
      cbind(group19, as.data.frame(t(Genes19[input$yvariable19, 2:25])))
    })
    toplot19 <- reactive({
      dat19 <- y19()
      colnames(dat19) <- c("group", "gene")
      dat19$group <- factor(dat19$group, levels = c("Untreated", "IL-13", "Omeprazole", "Esomeprazole", "IL-13+Omeprazole", "IL-13+Esomeprazole"))
      dat19
    })
    g19 <- reactive({
      graph19 <- ggplot(data = toplot19()) +
        ggtitle(input$yvariable19) +
        theme(aspect.ratio = as.numeric(2*gr19/n19), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph19
    })
    output$selectedgraph19 <- renderPlot({
      graph19 <- g19()
      if (input$graph19 == 1){graph19 <- g19() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph19 == 2){graph19 <- g19() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph19 == 3){graph19 <- g19() + aes(x=row.names(toplot19()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot19())) + scale_fill_manual(values = cbPalette)}
      graph19 + xlab(NULL) + ylab("Expression, TPM")
    })
    output$selectedanova19 <- renderPrint({
      x19 <- data.frame(toplot19())
      df19 <- round(as.data.frame(pairwise.t.test(x19$gene, x19$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df19 <- sapply(df19, as.character)
      df19[is.na(df19)] <- " "
      df19[df19 == "0"] <- "<0.001"
      row.names(df19) <- c("IL-13", "Omeprazole", "Esomeprazole", "IL-13+Omeprazole", "IL-13+Esomeprazole")
      colnames(df19) <- c("Untreated", "IL-13", "Omeprazole", "Esomeprazole", "IL-13+Omeprazole")
      as.data.frame(df19)
    })
    output$selecteddata19 <- renderTable({
      d19 <- as.data.frame(Genes19[input$yvariable19, 2:25])
      d19 <- rbind(colnames(d19), d19)
      row.names(d19) <- c("Sample", "Expression")
      d19 <- as.data.frame(t(d19))
    })
    output$downloadData19 <- downloadHandler(
      filename = paste("PPIsinEos.csv", sep = ""),
      content = function(file){
        write.csv(data_raw19, file, row.names = TRUE)}
    )
    output$sampleinfo19 <- renderTable({
      sample_info19
    })
    
    
    #server for 20
    y20 <- reactive({
      cbind(group20, as.data.frame(t(Genes20[input$yvariable20, 2:13])))
    })
    toplot20 <- reactive({
      dat20 <- y20()
      colnames(dat20) <- c("group", "gene")
      dat20$group <- factor(dat20$group, levels = c("Untreated", "Omeprazole", "Esomeprazole", "IL-13", "IL-13+Omeprazole", "IL-13+Esomeprazole"))
      dat20
    })
    g20 <- reactive({
      graph20 <- ggplot(data = toplot20()) +
        ggtitle(input$yvariable20) +
        theme(aspect.ratio = as.numeric(2*gr20/n20), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph20
    })
    output$selectedgraph20 <- renderPlot({
      graph20 <- g20()
      if (input$graph20 == 1){graph20 <- g20() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph20 == 2){graph20 <- g20() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph20 == 3){graph20 <- g20() + aes(x=row.names(toplot20()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot20())) + scale_fill_manual(values = cbPalette)}
      graph20 + xlab(NULL) + ylab("Expression")
    })
    output$selectedanova20 <- renderPrint({
      x20 <- data.frame(toplot20())
      df20 <- round(as.data.frame(pairwise.t.test(x20$gene, x20$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df20 <- sapply(df20, as.character)
      df20[is.na(df20)] <- " "
      df20[df20 == "0"] <- "<0.001"
      row.names(df20) <- c("Omeprazole", "Esomeprazole", "IL-13", "IL-13+Omeprazole", "IL-13+Esomeprazole")
      colnames(df20) <- c("Untreated", "Omeprazole", "Esomeprazole", "IL-13", "IL-13+Omeprazole")
      as.data.frame(df20)
    })
    output$selecteddata20 <- renderTable({
      d20 <- as.data.frame(Genes20[input$yvariable20, 2:13])
      d20 <- rbind(colnames(d20), d20)
      row.names(d20) <- c("Sample", "Expression")
      d20 <- as.data.frame(t(d20))
    })
    output$downloadData20 <- downloadHandler(
      filename = paste("PPIsinEosMultiplex.csv", sep = ""),
      content = function(file){
        write.csv(data_raw20, file, row.names = TRUE)}
    )
    output$sampleinfo20 <- renderTable({
      sample_info20
    })
#server for 23
    y23 <- reactive({
      cbind(group23, as.data.frame(t(Genes23[input$yvariable23, 2:13])))
    })
    toplot23 <- reactive({
      dat23 <- y23()
      colnames(dat23) <- c("group", "gene")
      dat23$group <- factor(dat23$group, levels = c("Untreated", "IL-13", "VitD", "IL-13_&_VitD"))
      dat23
    })
    g23 <- reactive({
      graph23 <- ggplot(data = toplot23()) +
        ggtitle(input$yvariable23) +
        theme(aspect.ratio = as.numeric(2*gr23/n23), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
      graph23
    })
    output$selectedgraph23 <- renderPlot({
      graph23 <- g23()
      if (input$graph23 == 1){graph23 <- g23() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
      if (input$graph23 == 2){graph23 <- g23() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
      if (input$graph23 == 3){graph23 <- g23() + aes(x=row.names(toplot23()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot23())) + scale_fill_manual(values = cbPalette)}
      graph23 + xlab(NULL) + ylab("Expression, RPKM")
    })
    output$selectedanova23 <- renderPrint({
      x23 <- data.frame(toplot23())
      df23 <- round(as.data.frame(pairwise.t.test(x23$gene, x23$group, p.adjust.method = "fdr")$p.value), digits = 3)
      df23 <- sapply(df23, as.character)
      df23[is.na(df23)] <- " "
      df23[df23 == "0"] <- "<0.001"
      row.names(df23) <- c("IL-13", "VitD", "IL-13_&_VitD")
      colnames(df23) <- c("Untreated", "IL-13", "VitD")
      as.data.frame(df23)
    })
    output$downloadData23 <- downloadHandler(
      filename = paste("EPC2_IL-13_&_VitD.csv", sep = ""),
      content = function(file){
        write.csv(data_raw23, file, row.names = TRUE)}
    )
    output$selecteddata23 <- renderTable({
      d23 <- as.data.frame(Genes23[input$yvariable23, 2:13])
      d23 <- rbind(colnames(d23), d23)
      row.names(d23) <- c("Sample", "Expression")
      d23 <- as.data.frame(t(d23))
    })
    output$sampleinfo23 <- renderTable({
      sample_info23
    })
#server for 24
  y24 <- reactive({
    cbind(group24, as.data.frame(t(Genes24[input$yvariable24, 2:22])))
  })
  toplot24 <- reactive({
    dat24 <- y24()
    colnames(dat24) <- c("group", "gene")
    dat24$group <- factor(dat24$group, levels = c("Non-EoG", "EoG"))
    dat24
  })
  g24 <- reactive({
    graph24 <- ggplot(data = toplot24()) +
    ggtitle(input$yvariable24) +
    theme(aspect.ratio = as.numeric(2*gr24/n24), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
  graph24
  })
  output$selectedgraph24 <- renderPlot({
    graph24 <- g24()
    if (input$graph24 == 1){graph24 <- g24() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
    if (input$graph24 == 2){graph24 <- g24() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
    if (input$graph24 == 3){graph24 <- g24() + aes(x=row.names(toplot24()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot24())) + scale_fill_manual(values = cbPalette)}
    graph24 + xlab(NULL) + ylab("Expression, TPM")
  })
  output$selectedanova24 <- renderPrint({
    x24 <- data.frame(toplot24())
    df24 <- round(as.data.frame(pairwise.t.test(x24$gene, x24$group, p.adjust.method = "fdr")$p.value), digits = 3)
    if (df24[1,1] == "0") {print("<0.001")}
    else {print(df24[1,1])}
  })
  output$selecteddata24 <- renderTable({
    d24 <- as.data.frame(Genes24[input$yvariable24, 2:22])
    d24 <- rbind(colnames(d24), d24)
    row.names(d24) <- c("Sample", "TPM")
    d24 <- as.data.frame(t(d24))
  })
  output$downloadData24 <- downloadHandler(
    filename = paste("EGtranscriptome.csv", sep = ""),
    content = function(file){
      write.csv(data_raw24, file, row.names = TRUE)}
  )
  output$sampleinfo24 <- renderTable({
    sample_info24
  })

#server for 25
  y25 <- reactive({
    cbind(group25, as.data.frame(t(Genes25[input$yvariable25, 2:7])))
  })
  toplot25 <- reactive({
    dat25 <- y25()
    colnames(dat25) <- c("group", "gene")
    dat25$group <- factor(dat25$group, levels = c("untreated_1hr", "untreated_4hr", "IL4_1hr", "IL4_4hr", "IL33_1hr", "IL33_4hr"))
    dat25
  })
  g25 <- reactive({
    graph25 <- ggplot(data = toplot25()) +
      ggtitle(input$yvariable25) +
      theme(aspect.ratio = as.numeric(2*gr25/n25), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
    graph25
  })
  output$selectedgraph25 <- renderPlot({
    graph25 <- g25()
    if (input$graph25 == 1){graph25 <- g25() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
    if (input$graph25 == 2){graph25 <- g25() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
    if (input$graph25 == 3){graph25 <- g25() + aes(x=row.names(toplot25()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot25())) + scale_fill_manual(values = cbPalette)}
    graph25 + xlab(NULL) + ylab("Expression, RPKM")
  })
  output$selecteddata25 <- renderTable({
    d25 <- as.data.frame(Genes25[input$yvariable25, 2:7])
    d25 <- rbind(colnames(d25), d25)
    row.names(d25) <- c("Sample", "RPKM")
    d25 <- as.data.frame(t(d25))
  })
  output$downloadData25 <- downloadHandler(
    filename = paste("mouseeos.csv", sep = ""),
    content = function(file){
      write.csv(data_raw25, file, row.names = TRUE)}
  )
  output$sampleinfo25 <- renderTable({
    sample_info25
  })
  
  #server for 26
  y26 <- reactive({
    cbind(group26, as.data.frame(t(Genes26[input$yvariable26, 2:15])))
  })
  toplot26 <- reactive({
    dat26 <- y26()
    colnames(dat26) <- c("group", "gene")
    dat26$group <- factor(dat26$group, levels = c("Non_EoC", "Active_EoC"))
    dat26
  })
  g26 <- reactive({
    graph26 <- ggplot(data = toplot26()) +
      ggtitle(input$yvariable26) +
      theme(aspect.ratio = as.numeric(2*gr26/n26), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
    graph26
  })
  output$selectedgraph26 <- renderPlot({
    graph26 <- g26()
    if (input$graph26 == 1){graph26 <- g26() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
    if (input$graph26 == 2){graph26 <- g26() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
    if (input$graph26 == 3){graph26 <- g26() + aes(x=row.names(toplot26()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot26())) + scale_fill_manual(values = cbPalette)}
    graph26 + xlab(NULL) + ylab("Expression, TPM")
  })
  output$selectedanova26 <- renderPrint({
    x26 <- data.frame(toplot26())
    df26 <- round(as.data.frame(pairwise.t.test(x26$gene, x26$group, p.adjust.method = "fdr")$p.value), digits = 3)
    if (df26[1,1] == "0") {print("<0.001")}
    else {print(df26[1,1])}
  })
  output$downloadData26 <- downloadHandler(
    filename = paste("EoC_RNAseq.csv", sep = ""),
    content = function(file){
      write.csv(data_raw26, file, row.names = TRUE)}
  )
  output$selecteddata26 <- renderTable({
    d26 <- as.data.frame(Genes26[input$yvariable26, 2:15])
    d26 <- rbind(colnames(d26), d26)
    row.names(d26) <- c("Sample", "Expression")
    d26 <- as.data.frame(t(d26))
  })
  
  #server for 28
  y28 <- reactive({
    cbind(group28, as.data.frame(t(Genes28[input$yvariable28, 2:6])))
  })
  toplot28 <- reactive({
    dat28 <- y28()
    colnames(dat28) <- c("group", "gene")
    dat28$group <- factor(dat28$group, levels = c("NoDox", "DOX"))
    dat28
  })
  g28 <- reactive({
    graph28 <- ggplot(data = toplot28()) +
      ggtitle(input$yvariable28) +
      theme(aspect.ratio = as.numeric(2*gr28/n28), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
    graph28
  })
  output$selectedgraph28 <- renderPlot({
    graph28 <- g28()
    if (input$graph28 == 1){graph28 <- g28() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
    if (input$graph28 == 2){graph28 <- g28() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
    if (input$graph28 == 3){graph28 <- g28() + aes(x=row.names(toplot28()), y=gene, fill=group) + geom_bar(stat="identity") + scale_x_discrete(limits=row.names(toplot28())) + scale_fill_manual(values = cbPalette)}
    graph28 + xlab(NULL) + ylab("Expression")
  })
  output$selectedanova28 <- renderPrint({
    x28 <- data.frame(toplot28())
    df28 <- round(as.data.frame(pairwise.t.test(x28$gene, x28$group, p.adjust.method = "fdr")$p.value), digits = 3)
    if (df28[1,1] == "0") {print("<0.001")}
    else {print(df28[1,1])}
  })
  output$selecteddata28 <- renderTable({
    d28 <- as.data.frame(Genes28[input$yvariable28, 2:6])
    d28 <- rbind(colnames(d28), d28)
    row.names(d28) <- c("Sample", "Expression")
    d28 <- as.data.frame(t(d28))
  })
  output$downloadData28 <- downloadHandler(
    filename = paste("Mouse_IL13Tg_esophagus.csv", sep = ""),
    content = function(file){
      write.csv(data_raw28, file, row.names = TRUE)}
  )
  output$sampleinfo28 <- renderTable({
    sample_info28
  })
  
  #server for 34
  y34 <- reactive({
    cbind(group34, as.data.frame(t(Genes34[input$yvariable34, 2:28])))
  })
  toplot34 <- reactive({
    dat34 <- y34()
    colnames(dat34) <- c("group", "gene")
    dat34$group <- factor(dat34$group, levels = c("NL", "DE", "EoD"))
    dat34
  })
  g34 <- reactive({
    graph34 <- ggplot(data = toplot34()) +
      ggtitle(input$yvariable34) +
      theme(aspect.ratio = as.numeric(10*gr34/n34), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14), plot.title = element_text(size = 20, face = "bold"))
    graph34
  })
  output$selectedgraph34 <- renderPlot({
    graph34 <- g34()
    if (input$graph34 == 1){graph34 <- g34() + aes(x=group,y=gene) + stat_summary(fun.y="mean", geom="bar") + stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.3)}
    if (input$graph34 == 2){graph34 <- g34() + aes(x=group,y=gene) + geom_boxplot() + geom_point()}
    if (input$graph34 == 3){graph34 <- g34() + aes(x=row.names(toplot34()), y=gene, color=group) + geom_point(stat="identity") + scale_x_discrete(limits=row.names(toplot34())) + scale_color_manual(values = cbPalette) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6))}
    graph34 + xlab(NULL) + ylab("Expression - TPM")
  })
  output$selectedanova34 <- renderPrint({
    x34 <- data.frame(toplot34())
    df34 <- round(as.data.frame(pairwise.t.test(x34$gene, x34$group, p.adjust.method = "fdr")$p.value), digits = 3)
    df34 <- sapply(df34, as.character)
    df34[is.na(df34)] <- " "
    df34[df34 == "0"] <- "<0.001"
    row.names(df34) <- c("DE", "EoD")
    colnames(df34) <- c("NL", "DE")
    as.data.frame(df34)
  })
  output$selecteddata34 <- renderTable({
    d34 <- as.data.frame(Genes34[input$yvariable34, 2:28])
    d34 <- rbind(colnames(d34), d34)
    row.names(d34) <- c("Sample", "Expression")
    d34 <- as.data.frame(t(d34))
  })
  output$downloadData34 <- downloadHandler(
    filename = paste("EoD.csv", sep = ""),
    content = function(file){
      write.csv(Genes34, file, row.names = TRUE)}
  )
  output$sampleinfo34 <- renderTable({
    sample_info34
  })
  
  #server for 35
  output$selectedgraph35 <- renderPlot({
    req(input$yvariable35)
    if (input$graph35 == 1){graph35 <- FeaturePlot(data35, features = input$yvariable35, label = TRUE, split.by = "Diagnosis")}
    if (input$graph35 == 2){graph35 <- VlnPlot(data35, features = input$yvariable35)}
    if (input$graph35 == 3){graph35 <- RidgePlot(data35, features = input$yvariable35)}
    graph35
  })
  
  #server for 36
  output$selectedgraph36 <- renderPlot({
    req(input$yvariable36)
    if (input$graph36 == 1){graph36 <- FeaturePlot(data36, features = input$yvariable36, label = TRUE, split.by = "Diagnosis")}
    if (input$graph36 == 2){graph36 <- VlnPlot(data36, features = input$yvariable36)}
    if (input$graph36 == 3){graph36 <- RidgePlot(data36, features = input$yvariable36)}
    graph36
  })

  #server for 37
  output$selectedgraph37 <- renderPlot({
    req(input$yvariable37)
    if (input$graph37 == 1){graph37 <- FeaturePlot(data37, features = input$yvariable37, label = TRUE, split.by = "group", min.cutoff = "q1")}
    if (input$graph37 == 2){graph37 <- VlnPlot(data37, features = input$yvariable37)}
    if (input$graph37 == 3){graph37 <- RidgePlot(data37, features = input$yvariable37)}
    graph37
  })
  
  #server for abbreviations
  output$abbreviations <- renderTable({
    abbreviations
  })
#this goes after last tab server 
}

# Run the application 
shinyApp(ui = ui, server = server)
