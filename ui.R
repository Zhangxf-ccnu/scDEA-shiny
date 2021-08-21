# ui.R
library(shiny)
library(shinyBS)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(scDEA)
library(graphics)
library(SummarizedExperiment)
library(Seurat)
# library(BiocInstall)
# library(clusterProfiler)
# library(networkD3)
# source("Visulization.R")
# source("tsne_cluster.R")
shinyUI(navbarPage("scDEA", id="navbar",
                   theme = shinytheme("flatly"),
                   fluidPage(
                     sidebarLayout(
                       sidebarPanel(

                         ###Parameters of twelve individual methods
                         fileInput("count",
                                   label = "Count matrix",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv")),
                         fileInput("cell_label",
                                   label = "Cell label",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv")),
                         numericInput("trimmed",
                                      label = "trimmed",
                                      value = 0.2, min =  0, max = 0.5, step = 0.1),
                         checkboxInput("weight", "weight", value = TRUE),


                         checkboxInput("BPSC", "BPSC", value = TRUE),
                         checkboxInput("DEsingle", "DEsingle", value = TRUE),
                         checkboxInput("DESeq2", "DESeq2", value = TRUE),
                         checkboxInput("edgeR", "edgeR", value = TRUE),
                         checkboxInput("MAST", "MAST", value = TRUE),
                         checkboxInput("monocle", "monocle", value = TRUE),
                         checkboxInput("scDD", "scDD", value = TRUE),
                         checkboxInput("Ttest", "Ttest", value = TRUE),
                         checkboxInput("Wilcoxon", "Wilcoxon", value = TRUE),
                         checkboxInput("limma", "limma", value = TRUE),
                         checkboxInput("Seurat", "Seurat", value = TRUE),
                         checkboxInput("zingeR.edgeR", "zingeR.edgeR", value = TRUE),


                         checkboxInput("is.normalized", "Is normalized data?", value = FALSE),
                         checkboxInput("verbose", "Save the DE results?", value = TRUE),
                         selectInput('BPSC.normalize', 'BPSC.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         numericInput("BPSC.coef", label = "BPSC.coef", value = 2, step = 1),
                         checkboxInput("BPSC.parallel", "BPSC.parallel", value = TRUE),

                         selectInput('DEsingle.normalize', 'DEsingle.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         checkboxInput("DEsingle.parallel", "DEsingle.parallel", value = TRUE),


                         selectInput('DESeq2.test', 'DESeq2.test', c("Wald", "LRT"), selected = "Wald"),
                         checkboxInput("DESeq2.parallel", "DESeq2.parallel", value = TRUE),
                         checkboxInput("DESeq2.beta.prior", "DESeq2.beta.prior", value = TRUE),
                         selectInput('DESeq2.fitType', 'DESeq2.fitType', c("parametric","local", "mean"), selected = "parametric"),
                         selectInput('DESeq2.normalize', 'DESeq2.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "TMM"),

                         selectInput('edgeR.Test', 'edgeR.Test', c("LRT", "QLFT"), selected = "QLFT"),
                         selectInput('edgeR.normalize', 'edgeR.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "TMM"),

                         selectInput('limma.method.fit', 'limma.method.fit', c("ls", "robust"), selected = "ls"),
                         checkboxInput("limma.trend", "limma.trend", value = TRUE),
                         checkboxInput("limma.robust", "limma.robust", value = TRUE),
                         selectInput('limma.normalize', 'limma.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),

                         selectInput('Seurat.normalize', 'Seurat.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         selectInput('Seurat.method', 'Seurat.method', c("bimod", "negbinom", "poisson", "LR"), selected = "bimod"),

                         selectInput('MAST.method', 'MAST.method', c("bayesglm", "glm", "glmer"), selected = "bayesglm"),
                         selectInput('MAST.normalize', 'MAST.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         checkboxInput("MAST.parallel", "MAST.parallel", value = TRUE),

                         selectInput('monocle.normalize', 'monocle.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         # checkboxInput("monocle.", "monocle.parallel", value = TRUE),
                         numericInput("monocle.cores", label = "monocle.cores", value = 1, min = 1, max = 8, step = 1),

                         selectInput('scDD.normalize', 'scDD.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         numericInput("scDD.alpha1", label = "scDD.alpha1", value = 0.01, step = 0.1),
                         numericInput("scDD.mu0", label = "scDD.mu0", value = 0, step = 0.01),
                         numericInput("scDD.s0", label = "scDD.s0", value = 0.01, step = 0.01),
                         numericInput("scDD.a0", label = "scDD.a0", value = 0.01, step = 0.01),
                         numericInput("scDD.b0", label = "scDD.b0", value = 0.01, step = 0.01),
                         numericInput("scDD.permutation", label = "scDD.permutation", value = 0, step = 1),

                         selectInput('Ttest.normalize', 'Ttest.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),

                         selectInput('Wilcoxon.normalize', 'Wilcoxon.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),

                         selectInput('zingeR.edgeR.normalize', 'zingeR.edgeR.normalize', c("TMM", "RLE", "CPM", "TPM"), selected = "CPM"),
                         numericInput("zingeR.edgeR.maxit.EM", label = "zingeR.edgeR.maxit.EM", value = 100, min = 0, max = 200, step = 1),
                         ###Parameters of twelve individual methods

                         ### Parameters of ensemble methods



                         actionButton("Run_scDEA",
                                      label = "Run scDEA",
                                      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),

                         width = 3
                       ),
                       mainPanel(
                         fluidRow(
                           column(9,
                                  plotOutput("Heatmap_scDEA"),
                                  br(),
                                  textOutput("summary1"),
                                  br(),
                                  textOutput("summary2"),
                                  br(),
                                  textOutput("summary3"),
                                  br(),
                                  textOutput("summary4")
                           ),
                           column(3,
                                  numericInput("number_genes_show_scDEA",
                                               label = "number genes show",
                                               value = 10, min = 1, step = 1),
                                  radioButtons("fileformat",
                                               label = "Output File Format",
                                               choices = c(".txt (tab-delimited text)" = "txt",
                                                           ".csv (comma-separated values)" = "csv"),
                                               selected = "csv"),
                                  textInput("dlname",
                                            label = "Output File Name (Do not include file extension)"),
                                  downloadButton("download",
                                                 label = "Download")
                           )
                         ),
                         fluidRow(
                           DT::dataTableOutput("result_scDEA")),

                       )
                   )
)
)
)
