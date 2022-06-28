library(dplyr)
library(purrr)

#####
# BLAST column header
blast_header <- 
  c("ctg_cds", "NCBI_accession", "evalue", "bitscore", "length", "pident", "ppos", 
    "qcovs", "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "staxid", 
    "sskingdom", "sscinames", "scomnames", "NCBI_annotation")
#####

# Read the statistics from UBC DEA
ubc_sigDE <- read.table("results/ubc_dea_sig_results.txt", 
                        header = TRUE, stringsAsFactors = FALSE)

colnames(ubc_sigDE)[1] <- "ctg_cds"

str(ubc_sigDE)
# 'data.frame':	2399 obs. of  3 variables:


# Read BLAST results  
ubc_blast_results <- 
  read.table("results/ubc_dea_sig.blastpRefSeqPlant212.txt", sep = "\t", 
             quote = "", stringsAsFactors = FALSE, col.names = blast_header)

str(ubc_blast_results)
# 'data.frame':	12031 obs. of  20 variables:

# Only keep relevant columns
ubc_blast_results <- ubc_blast_results %>% 
  select(c(ctg_cds, NCBI_accession, NCBI_annotation))

# Take the top blast hit
ubc_top_blast_result <- ubc_blast_results %>% 
  group_by(ctg_cds) %>% 
  slice_head(n = 1)

# Join statistics with annotation
ubc_sigDE_annotated <- 
  left_join(ubc_sigDE, ubc_top_blast_result, by = "ctg_cds")
# 'data.frame':	2399 obs. of  5 variables:


table(is.na(ubc_sigDE_annotated$salltitles))
# FALSE  TRUE
#  1164  1235 
# 1164 sequences have no annotation

write.table(ubc_sigDE_annotated, "results/ubc_dea_sig_ncbi-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#####

# Read the statistics from JGI DEA
jgi_sigDE <- read.table("results/jgi_dea_sig_results.txt", 
                        header = TRUE, stringsAsFactors = FALSE)

colnames(jgi_sigDE)[1] <- "ctg_cds"

str(jgi_sigDE)
# 'data.frame':	22404 obs. of  4 variables:


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

# Join statistics with annotation
jgi_sigDE_annotated <- 
  left_join(jgi_sigDE, jgi_top_blast_result, by = c("ctg_cds"))

str(jgi_sigDE_annotated)
# 'data.frame':	22404 obs. of  19 variables:


table(is.na(jgi_sigDE_annotated$salltitles))
# FALSE  TRUE 
# 16063  6341 
# 6341 sequences have no annotation

write.table(jgi_sigDE_annotated, "results/jgi_dea_sig_ncbi-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
