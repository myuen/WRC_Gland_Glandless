library(dplyr)

#####
# BLAST column header
blast_header <- 
  c("ctg_cds", "NCBI_accession", "evalue", "bitscore", "length", "pident", "ppos", 
    "qcovs", "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "staxid", 
    "sskingdom", "sscinames", "scomnames", "NCBI_annotation")
#####

# Read BLAST results  
ubc_ncbi_blast_results <- 
  read.table("results/ubc_dea_sig.blastpRefSeqPlant212.txt", sep = "\t", 
             quote = "", stringsAsFactors = FALSE, col.names = blast_header)

str(ubc_ncbi_blast_results)
# 'data.frame':	12031 obs. of  20 variables:

# Only keep relevant columns
ubc_ncbi_blast_results <- ubc_ncbi_blast_results %>% 
  select(c(ctg_cds, NCBI_accession, NCBI_annotation))

# Take the top blast hit
ubc_top_blast_result <- ubc_ncbi_blast_results %>% 
  group_by(ctg_cds) %>% 
  slice_head(n = 1)

write.table(ubc_top_blast_result, "results/ubc_dea_sig_ncbi-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#####

# Read BLAST results  
jgi_blast_results <- 
  read.table("results/jgi_dea_sig.blastpRefSeqPlant212.txt", sep = "\t", 
             quote = "", stringsAsFactors = FALSE, col.names = blast_header)

str(jgi_blast_results)
# 'data.frame':	159030 obs. of  20 variables:

# Only keep relevant columns
jgi_blast_results <- jgi_blast_results %>% 
  select(c(ctg_cds, NCBI_accession, NCBI_annotation))

# Take the top blast hit
jgi_top_blast_result <- jgi_blast_results %>% 
  group_by(ctg_cds) %>% 
  slice_head(n = 1)

write.table(jgi_top_blast_result, "results/jgi_dea_sig_ncbi-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
