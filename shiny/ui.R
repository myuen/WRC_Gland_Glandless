library(DT)
library(shiny)


fluidPage(
  tabsetPanel(
    tabPanel('Differential Expression',
             fluidRow(
               column(4, 
                      (h4('Aggregated differential expression analysis results of UBC glanded genotype 5314 against glandless genotype 5038 and JGI glanded genotype 5309 against glandless genotype 5038')),
                      br(),
                      br(),
                      
                      DT::DTOutput('DE'),
               ),
               
               column(4, offset = 3,
                      # DT::DTOutput('annot'),
                      tableOutput('annot'),
                      br(),
                      br(),
                      # verbatimTextOutput('seq')
                      textOutput('seq')
               )
             )
    ),
    
    tabPanel('Differential Expression by set',
             fluidRow(
               column(5,
                      (h4('Set presentataion on differential expression results on UBC glanded genotype 5314 and JGI glanded genotype 5309 against glandless genotype 5038')),
                      br(),
                      column(4,

                      selectInput(input = 'UBC_exp',
                                  label = 'UBC expression relative to glandless genotype',
                                  choices = c('Up' = 'up',
                                              'Down' = 'down')),
                      ),

                      column(4,
                      selectInput(input = 'JGI_exp',
                                  label = 'JGI expression relative to glandless genotype',
                                  choices = c('Up' = 'up',
                                              'Down' = 'down')),
                      ),

                      DT::DTOutput('DE_by_set'),
               ),

               column(1, offset = 1,
                      tableOutput('annot_by_set')
                      )
             ),
    )   
  )
)