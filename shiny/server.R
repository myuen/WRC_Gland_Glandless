library(dplyr)
library(DT)
library(magrittr)
library(readr)
library(reticulate)
library(shiny)
library(tibble)


seqio <- import('Bio.SeqIO')

### Read differential expression statistics
# ubc <-
#   read.table("results/UBC/ubc_dea_results.sigDE.txt", header = TRUE)
ubc <- 
  read_table("results/UBC/ubc_dea_results.sigDE.txt", show_col_types = FALSE)
# jgi <-
#   read.table("results/UBC/jgi_dea_on_ubc_transcriptome.sigDE.txt", header=TRUE)
jgi <-
  read_table("results/UBC/jgi_dea_on_ubc_transcriptome.sigDE.txt", 
             show_col_types = FALSE)


### Join the results so we can compare how the expression of the UBC gland 
### genotype differ from JGI gland genotype
ubc_jgi <- full_join(ubc, jgi, by='cds',
                     suffix = c(' UBC_5309', ' JGI_5314'))

ubc_jgi <- ubc_jgi |>
  select('adj.P.Val JGI_5314', 'logFC JGI_5314', 
         'cds', 
         'logFC UBC_5309', 'adj.P.Val UBC_5309') 

ubc_jgi <- ubc_jgi[order(ubc_jgi$`logFC UBC_5309`, decreasing = TRUE),]

ubc_jgi_cds_id <- unique(ubc_jgi$cds)


### Read fasta file through reticulate
nt_dict <- seqio$index('data/ubc_gland_glandless_trinity_23Aug2021.transdecoder.cds.nr.fasta', 'fasta')
aa_dict <- seqio$index('data/ubc_gland_glandless_trinity_23Aug2021.transdecoder.pep.fasta', 'fasta')


### Read annotations
refseq <- read_delim('data/all.sigDE.blastpRefSeqPlant213.tophit.txt', 
                     delim = "\t", show_col_types = FALSE)
refseq$source <- 'refseq'

swissprot <- read_delim('data/all.sigDE.blastpSwissProt.tophit.txt', 
                        delim = "\t", show_col_types = FALSE)
swissprot$source <- 'swissprot'

  
function(input, output) {
  
  # Output for the first tab
  output$DE <- DT::renderDT(ubc_jgi, server = TRUE,
                            selection = 'single', rownames = FALSE)

  # Annotation on the side for the first tab
  output$annot <- renderTable({
    id <- ubc_jgi[input$DE_rows_selected, 'cds']
    selected = input$DE_rows_selected
    
    if (length(selected)) {
      refseq_selected <- refseq |> 
        filter(qseqid == id)
      
      swissprot_selected <- swissprot |> 
        filter(qseqid == id)

      refseq_count <- as.logical(length(refseq_selected))
      swissprot_count <- as.logical(length(swissprot_selected))
      
      annot <- rbind(refseq_selected, swissprot_selected)
      
      annot <- annot |> 
        select(source, qseqid, sseqid, evalue, salltitles)

    }
  }, rownames = FALSE, colnames = FALSE, align = 'c')
  
  
  # Sequence output on the side for the first tab
  output$seq = renderText({
    selected = input$DE_rows_selected
    
    id <- ubc_jgi[input$DE_rows_selected, 'cds']
    
    if (length(selected)) {
      seqs <- paste(nt_dict[id]$format('fasta'), aa_dict[id]$format('fasta'), sep = '')
    }
  })
  
  output$Annot <- DT::renderDT(annot, server = TRUE, 
                               selection = 'single', rownames = FALSE)

  output$DE <- DT::renderDT(ubc_jgi, server = TRUE,
                            selection = 'single', rownames = FALSE)
  
  output$annot <- renderTable({
    id <- ubc_jgi[input$DE_rows_selected, 'cds']
    selected = input$DE_rows_selected
    
    if (length(selected)) {
      refseq_selected <- refseq |> 
        filter(qseqid == id)
      
      swissprot_selected <- swissprot |> 
        filter(qseqid == id)
      
      refseq_count <- as.logical(length(refseq_selected))
      swissprot_count <- as.logical(length(swissprot_selected))
      
      annot <- rbind(refseq_selected, swissprot_selected)
      
      annot <- annot |> 
        select(source, qseqid, sseqid, evalue, salltitles)
      
    }
  }, rownames = FALSE, colnames = FALSE, align = 'c')
  # End first tab

  
  # Output for second tab
  output$DE_by_set <- DT::renderDT({
    
    # Remove all rows with NA
    ubc_jgi_noNA <- ubc_jgi |> 
      filter(!is.na(`logFC UBC_5309`) & !is.na(`logFC JGI_5314`))

    if(input$UBC_exp == 'up' & input$JGI_exp == 'up') {
      ubc_jgi_noNA |> 
        filter(`logFC UBC_5309` >=2 & `logFC JGI_5314` >= 2)
      
      } else if (input$UBC_exp == 'up' & input$JGI_exp == 'down') {
        ubc_jgi_noNA |> 
          filter(`logFC UBC_5309` >= 2 & `logFC JGI_5314` <= -2) 
        
      } else if (input$UBC_exp == 'down' & input$JGI_exp == 'up') {
        ubc_jgi_noNA |> 
          filter(`logFC UBC_5309` <= -2 & `logFC JGI_5314` >= 2)
        
      } else if (input$UBC_exp == 'down' & input$JGI_exp == 'down') {
        ubc_jgi_noNA |> 
          filter(`logFC UBC_5309` <= -2 & `logFC JGI_5314` <= -2)
      }
  }, server = TRUE, selection = 'single', rownames = FALSE)
}
