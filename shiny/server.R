library(dplyr)
library(DT)
library(magrittr)
library(readr)
library(reticulate)
library(shiny)
library(tibble)

# Import SeqIO from Bioconductor
seqio <- import('Bio.SeqIO')

### Read differential expression statistics
ubc <-
  read.delim("results/UBC/ubc_dea_results.sigDE.txt", sep = "\t")
jgi <-
  read.delim("results/UBC/jgi_dea_on_ubc_transcriptome.sigDE.txt", sep = "\t")



### Join the results so we can compare how the expression of the UBC gland 
### genotype differ from JGI gland genotype
ubc_jgi <- full_join(ubc, jgi, by='cds',
                     suffix = c('_UBC_5309', '_JGI_5314'))

ubc_jgi <- ubc_jgi |>
  select(adj.P.Val_JGI_5314, logFC_JGI_5314, 
         cds, 
         logFC_UBC_5309, adj.P.Val_UBC_5309) 

ubc_jgi <- ubc_jgi[order(ubc_jgi$logFC_UBC_5309, decreasing = TRUE),]

ubc_jgi_cds_id <- unique(ubc_jgi$cds)


### Read fasta file through reticulate
nt_dict <- seqio$index('data/ubc_gland_glandless_trinity_23Aug2021.transdecoder.cds.nr.fasta', 'fasta')
aa_dict <- seqio$index('data/ubc_gland_glandless_trinity_23Aug2021.transdecoder.pep.fasta', 'fasta')


### Read RefSeq annotations
refseq <- read_delim('data/all.sigDE.blastpRefSeqPlant213.tophit.txt', 
                     delim = "\t", show_col_types = FALSE)
refseq$source <- 'refseq'


# Read SwissProt annotations
swissprot <- read_delim('data/all.sigDE.blastpSwissProt.tophit.txt', 
                        delim = "\t", show_col_types = FALSE)
swissprot$source <- 'swissprot'

  
function(input, output) {
  
  # Output for the first tab
  
  # Output DEA statistics
  output$DE <- DT::renderDT(ubc_jgi, server = TRUE,
                            selection = 'single', rownames = FALSE)

  # Put annotation and sequences on the right side
  output$annot <- renderTable({
    
    selected <- input$DE_rows_selected
    id <- ubc_jgi[selected, 'cds']
    
    
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
  }, rownames = FALSE, colnames = FALSE, align = 'l')
  
  
  # Sequence output on the side for the first tab
  output$seq = renderText({

    selected <- input$DE_rows_selected
    id <- ubc_jgi[selected, 'cds']

    if (length(selected)) {
      seqs <- paste(nt_dict[id]$format('fasta'), aa_dict[id]$format('fasta'), sep = '')
    }
  })
  # End first tab

  
  # Output for second tab
  output$DE_by_set <- DT::renderDT({

    # Remove all rows with NA
    ubc_jgi_noNA <- ubc_jgi |>
      filter(!is.na(logFC_UBC_5309) & !is.na(logFC_JGI_5314))

    if(input$UBC_exp == 'up' & input$JGI_exp == 'up') {
      ubc_jgi_noNA |>
        filter(logFC_UBC_5309 >=2 & logFC_JGI_5314 >= 2)

      } else if (input$UBC_exp == 'up' & input$JGI_exp == 'down') {
        ubc_jgi_noNA |>
          filter(logFC_UBC_5309 >= 2 & logFC_JGI_5314 <= -2)

      } else if (input$UBC_exp == 'down' & input$JGI_exp == 'up') {
        ubc_jgi_noNA |>
          filter(logFC_UBC_5309 <= -2 & logFC_JGI_5314 >= 2)

      } else if (input$UBC_exp == 'down' & input$JGI_exp == 'down') {
        ubc_jgi_noNA |>
          filter(logFC_UBC_5309 <= -2 & logFC_JGI_5314 <= -2)
      }
  }, server = TRUE, selection = 'single', rownames = FALSE)


  output$annot_by_set <- renderTable({
    selected <- input$DE_by_set_rows_selected
    id <- ubc_jgi[selected, 'cds']


    if (length(selected)) {
      refseq_selected <- refseq |>
        filter(qseqid == id)

      swissprot_selected <- swissprot |>
        filter(qseqid == id)

      refseq_count <- as.logical(length(refseq_selected))
      swissprot_count <- as.logical(length(swissprot_selected))

      annot <- rbind(refseq_selected, swissprot_selected)
      
      annot |>
        select(source, qseqid, sseqid, evalue, salltitles)

    }
  }, rownames = FALSE, colnames = FALSE, align = 'l')
}
