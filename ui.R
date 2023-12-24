library(plotly)
library(shiny)
library(ggplot2)

################################################################################
########### Samplemeta #########################################################
################################################################################

tab_samplemeta <- #sidebarLayout(
  
#  sidebarPanel(
    
  fluidPage(
    selectInput(
      inputId = "samplemeta_pool",
      label = "Pool:",
      selectize = FALSE,
      multiple = FALSE,
      choices = names(all_samplemeta), 
      selected = names(all_samplemeta)[1]
    ),

 # ),
#  mainPanel(
    #h3("UMAPs pool overview"),
    plotOutput(outputId = "plotSamplemetaUmap", height = "700px")
  #)
)


################################################################################
########### Gene viewer ########################################################
################################################################################

tab_grstats <- fluidPage(
    column(6,
      wellPanel(
        
        plotlyOutput("plot_grstats_volcano", height = "400px"),
        
        selectInput(
          inputId = "grstats_pool",
          label = "Pool:",
          selectize = FALSE,
          multiple = FALSE,
          choices = names(all_grstats), 
          selected = names(all_grstats)[1]
        ),
        
        selectInput(
          inputId = "grstats_volcano",
          label = "Volcano:",
          selectize = TRUE,
          multiple = FALSE,
          choices = c(""), 
          selected = NULL
        )
        

      )
    ),
    column(6,
       wellPanel(
         
         plotlyOutput("plot_grstats_scatterplot", height = "400px"),
  
         selectInput(
           inputId = "grstats_scatter",
           label = "Compare:",
           selectize = TRUE,
           multiple = FALSE,
           choices = c(""), 
           selected = NULL
         ),

         selectInput(
           inputId = "grstats_scatter_type",
           label = "Representation:",
           selectize = FALSE,
           multiple = FALSE,
           choices = c("Volcano plot","FC scatter plot"), 
           selected = "Volcano plot"
         ),
         
         
       )
    ),
    
    column(12,
       fluidPage(
         column(6,
            wellPanel(
              plotlyOutput("plot_grstats_tcplot", height = "400px")
            )
         ),
         column(6,
            wellPanel(
              selectInput(
                inputId = "grstats_gene",
                label = "Gene:",
                selectize = TRUE,
                multiple = FALSE,
                choices = c(""), 
                selected = NULL
              ),
              
              checkboxInput(
                inputId = "grstats_avg_mouse",
                label = "Average across mice",
                value = TRUE
              ),
              
              checkboxInput(
                inputId = "grstats_avg_grna",
                label = "Average across genetic constructs",
                value = TRUE
              ),

              checkboxInput(
                inputId = "grstats_avg_genotype",
                label = "Average across genotypes",
                value = FALSE
              ),

              checkboxInput(
                inputId = "grstats_avg_treatment",
                label = "Average across treatments",
                value = FALSE
              ),
              
              selectInput(
                inputId = "grstats_colorby",
                label = "Color by:",
                multiple = FALSE,
                choices = c("Gene","Mouse","Genotype","Treatment","Genotype+Treatment","Genotype+Treatment+Genetic construct"), 
                selected = c("Gene")
              ),
              
              selectInput(
                inputId = "grstats_units",
                label = "Units:",
                multiple = FALSE,
                selectize = TRUE,
                choices = c("Count/AllCount","Count/ControlCount"), 
                selected = c("Count/AllCount")
              ),
            )
         )
       )
    ),

)

################################################################################
########### About page #########################################################
################################################################################

tab_about <- verbatimTextOutput("todo description of project here")




################################################################################
########### Total page #########################################################
################################################################################


ui <- fluidPage(
  
  titlePanel("Bushell lab malaria screen viewer"),

  tabsetPanel(type = "tabs",
              tabPanel("Gene RGRs", tab_grstats),
              tabPanel("Pool overviews", tab_samplemeta),
              tabPanel("About", tab_about),
  )
  
)



