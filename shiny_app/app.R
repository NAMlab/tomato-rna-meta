library(shiny)
library(plotly)
library(DT)

source("./config.R")

pca.panel <- function(prefix) {
  return(
    sidebarLayout(
      sidebarPanel(
        numericInput(inputId = paste0(prefix, ".x.axis"),
                     label = "PC on X Axis",
                     min = 1,
                     max = 338,
                     value = 1),
        numericInput(inputId = paste0(prefix, ".y.axis"),
                     label = "PC on Y Axis",
                     min = 1,
                     max = 338,
                     value = 2),
        selectInput(inputId = paste0(prefix, ".col.factor"),
                    label = "Color points by...",
                    choices = c("Tissue", "Genotype", "Data Source", "Treatment"))
      ),
      mainPanel(
        plotlyOutput(outputId = paste0(prefix, ".plot"))
      )
    )
  )
}

table.panel <- function(title, id) {
  return(
    fluidPage(
      titlePanel(title),
      DT::dataTableOutput(id)
    )
  )
}

intersect.panel <- function(prefix, displaying="genes") {
  return(
    sidebarLayout(
      sidebarPanel(
        p(paste0("Only display ", displaying ," which are part of all following sets:")),
        h4("Union Sets (Figure 3)"),
        checkboxInput(inputId = paste0(prefix, ".union.up.heat"),
                      label = "Heat (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".union.down.heat"),
                      label = "Heat (Downregulated)"),
        checkboxInput(inputId = paste0(prefix, ".union.up.drought"),
                      label = "Drought (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".union.down.drought"),
                      label = "Drought (Downregulated)"),
        checkboxInput(inputId = paste0(prefix, ".union.up.salt"),
                      label = "Salt (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".union.down.salt"),
                      label = "Salt (Downregulated)"),
        h4("Core Sets (Figure 4)"),
        checkboxInput(inputId = paste0(prefix, ".core.up.heat"),
                      label = "Heat (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".core.down.heat"),
                      label = "Heat (Downregulated)"),
        checkboxInput(inputId = paste0(prefix, ".core.up.drought"),
                      label = "Drought (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".core.down.drought"),
                      label = "Drought (Downregulated)"),
        checkboxInput(inputId = paste0(prefix, ".core.up.salt"),
                      label = "Salt (Upregulated)"),
        checkboxInput(inputId = paste0(prefix, ".core.down.salt"),
                      label = "Salt (Downregulated)"),
      ),
      mainPanel(
        DT::dataTableOutput(paste0(prefix, ".table"))
      )
    )
  )
}

ui_home <- function() {
  return(
    fluidPage(
      titlePanel("Welcome!"),
      p("This Shiny app accompanies our RNA-Seq Meta-Analysis on Heat Stress in Tomato manuscript (DOI forthcoming).
        You can use it to review the samples used in more detail, inspect and/or download the raw mapped and diffexp data yourself,
        (through the 'Tables' menu), explore the PCAs in all their dimensions with interactive hover-over info for each data point ('PCA' menu),
        or find out which genes and GO terms are common to the different stress response sets ('Set Intersects' menu),
        i.e. which genes and GO terms are behind the bars in the upset plots in Figures 3 and 4."),
      p("The code for this Shiny app can be found in our repository along the other analyses at", a(href='https://github.com/NAMlab/tomato-rna-meta', "https://github.com/NAMlab/tomato-rna-meta"),
        "under ", code("shiny_app/"), ".")
    )
  )
}

ui <- navbarPage(
  "Tomato RNA-Seq Meta Analysis",
  tabPanel("Home", ui_home()),
  navbarMenu("Tables",
    tabPanel("Samples", value="table.samples", table.panel("Samples Overview", "table.samples")),
    tabPanel("Raw Data", value="table.raw.data", table.panel("Raw Data (big - may take a moment to load)", "table.raw.data")),
    ),
  navbarMenu("PCA", 
    tabPanel("logTPMs", value = "pca.logtpms", pca.panel("pca.logtpms")),
    tabPanel("logFCs", value = "pca.logfcs", pca.panel("pca.logfcs"))
  ),
  navbarMenu("Set Intersects",
    tabPanel("Genes", value="intersects.genes", intersect.panel("intersect.genes")),
    tabPanel("GO Terms", value="intersects.go", intersect.panel("intersect.go", "GO Terms"))),
)

server <- function(input, output) {
  load("data/scores_tpms.lfs.Rda")
  output$pca.logtpms.plot <- renderPlotly({
    col.factor = switch(input$pca.logtpms.col.factor,
                        "Tissue" = "tissue",
                        "Genotype" = "genotype.name",
                        "Data Source" = "datasource_id",
                        "Treatment" = "stress_type")
    col.palette = switch(input$pca.logtpms.col.factor,
                         "Tissue" = "tissue",
                         "Genotype" = "genotype",
                         "Data Source" = "datasource",
                         "Treatment" = "stress")

    plot1 <- plot_ly(scores.tpms, 
                     x=~get(paste0("PC", input$pca.logtpms.x.axis)), 
                     y=~get(paste0("PC", input$pca.logtpms.y.axis)),
                     color = ~get(col.factor),
                     colors = colors[[col.palette]],
                     mode = "markers", type = "scatter",
                     text = ~paste0("Sample: ",sample,'<br>Data Source: ', datasource_id,'<br>Genotype: ', genotype.name,
                                   '<br>Tissue: ', tissue, " (", tissue.comment, ")",'<br>Stress Type: ', stress_type,
                                   '<br>Temperature: ', temperature,'<br>Stress Duration: ', stress.duration,
                                   '<br>Sample Group: ', sample.group)) %>%
      layout(xaxis = list(title = paste("PC", input$pca.logtpms.x.axis)), 
            yaxis = list(title = paste("PC", input$pca.logtpms.y.axis)))
  })
  
  load("data/scores_logfc.lfs.Rda")
  output$pca.logfcs.plot <- renderPlotly({
    col.factor = switch(input$pca.logfcs.col.factor,
                        "Tissue" = "tissue",
                        "Genotype" = "genotype.name",
                        "Data Source" = "datasource_id",
                        "Treatment" = "stress_type")
    col.palette = switch(input$pca.logfcs.col.factor,
                         "Tissue" = "tissue",
                         "Genotype" = "genotype",
                         "Data Source" = "datasource",
                         "Treatment" = "stress")
    
    plot1 <- plot_ly(scores.logfc, 
                     x=~get(paste0("PC", input$pca.logfcs.x.axis)), 
                     y=~get(paste0("PC", input$pca.logfcs.y.axis)),
                     color = ~get(col.factor),
                     colors = colors[[col.palette]],
                     mode = "markers", type = "scatter",
                     text = ~paste0("Contrast: ",contrast,'<br>Data Source: ', datasource_id,'<br>Genotype: ', genotype.name,
                                    '<br>Tissue: ', tissue,'<br>Stress Type: ', stress_type,
                                    '<br>Temperature: ', temperature,'<br>Stress Duration: ', stress.duration)) %>%
      layout(xaxis = list(title = paste("PC", input$pca.logfcs.x.axis)), 
             yaxis = list(title = paste("PC", input$pca.logfcs.y.axis)))
  })
  
  samples.annotation = read.csv("data/samples_annotation.csv")
  output$table.samples <- DT::renderDataTable(samples.annotation,
                      options = list(buttons = c('csv', 'excel'), dom = 'Blfrtip', lengthMenu=c(10,25,50,100,250,1000,-1)),
                      filter = "top",
                      extensions = 'Buttons',
                      rownames = F)
  
  load("data/raw_data.lfs.Rda")
  output$table.raw.data <- DT::renderDataTable(raw.data,
                                              options = list(buttons = c('csv', 'excel'), dom = 'Blfrtip', lengthMenu=c(10,25,50,100,250,1000,-1)),
                                              filter = "top",
                                              extensions = 'Buttons',
                                              rownames = F)
  
  load("data/intersects_genes.lfs.Rda")
  output$intersect.genes.table <- renderDT({
    ret = intersects.genes
    print(input)
    for(f in apply(expand.grid(c("union", "core"), c("up", "down"), c("heat", "drought", "salt")), 1, paste, collapse=".")) {
      if(input[[paste0("intersect.genes.", f)]]) {
        print(f)
        ret = ret[ret[[f]],]
      }
    }
    ret = ret[1:4]
    datatable(ret, options = list(buttons = c('csv', 'excel'), dom = 'Blfrtip', lengthMenu=c(10,25,50,100,250,1000,-1)),
              filter = "top",
              extensions = 'Buttons',
              rownames = F)
  })
  
  load("data/intersects_gos.lfs.Rda")
  output$intersect.go.table <- renderDT({
    ret = intersects_gos
    print(input)
    for(f in apply(expand.grid(c("union", "core"), c("up", "down"), c("heat", "drought", "salt")), 1, paste, collapse=".")) {
      if(input[[paste0("intersect.go.", f)]]) {
        print(f)
        ret = ret[ret[[f]],]
      }
    }
    ret = ret[1:4]
    datatable(ret, options = list(buttons = c('csv', 'excel'), dom = 'Blfrtip', lengthMenu=c(10,25,50,100,250,1000,-1)),
              filter = "top",
              extensions = 'Buttons',
              rownames = F)
  })
}

shinyApp(ui = ui, server = server)