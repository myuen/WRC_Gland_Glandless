library(DT)
library(shiny)


fluidPage(
  fluidRow(
    column(4, 
           (h4('Aggregated differential expression analysis results of UBC glanded genotype 5314 against glandless genotype 5038 and JGI glanded genotype 5309 against glandless genotype 5038')),
           br(),
           br(),
           
           DT::DTOutput('DE'),
    ),
    
    column(4, offset = 4,
           # DT::DTOutput('annot'),
           tableOutput('annot'),
           br(),
           br(),
           # verbatimTextOutput('seq')
           textOutput('seq')
    )
  )
)
