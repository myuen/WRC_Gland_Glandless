library(dplyr)
library(purrr)
library(stringr)


#Read condensed gff3 file
gff3_mRNA <- read.delim("data/gff3_condensed.txt", 
                        sep = "\t", stringsAsFactors = FALSE)

str(gff3_mRNA)
# 'data.frame':	39659 obs. of  4 variables:


### Read DE statistics
UBC_de_stats <- read.delim("results/ubc_dea_sig_results.txt", 
                           stringsAsFactors = FALSE)

colnames(UBC_de_stats)[1] <- "ctg_cds"

str(UBC_de_stats)
# 'data.frame':	2399 obs. of  3 variables:


### Read UBC DE CDS with NCBI annotations
UBC_de_ncbi_annot <- read.delim("results/ubc_dea_sig_ncbi-annotated.txt")

str(UBC_de_ncbi_annot)
# 'data.frame':	1164 obs. of  3 variables:


### Read UBC DE CDS with genome annotations
UBC_de_genome_annot <- read.delim("results/ubc_dea_sig_genome-annotated.txt",
                                  stringsAsFactors = FALSE)

str(UBC_de_genome_annot)
# 'data.frame':	7850 obs. of  6 variables:


### Aggregating annotations
UBC_annots <- 
  left_join(UBC_de_stats, UBC_de_genome_annot, by = "ctg_cds")

UBC_annots <- 
  left_join(UBC_annots, gff3_mRNA, by = "mRNA")

UBC_annots <- UBC_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "mRNA", "start", "end",
  "domain_id", "domain_source", "domain_desc")
  
UBC_annots <- 
  left_join(UBC_annots, UBC_de_ncbi_annot, by = "ctg_cds")

write.table(UBC_annots, "results/ubc_dea_sig_aggregated-annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


### Repeat workflow on JGI dataset
JGI_de_stats <- read.delim("results/jgi_dea_sig_results.txt", 
                           stringsAsFactors = FALSE)

colnames(JGI_de_stats)[1] <- "ctg_cds"

JGI_young_de <- JGI_de_stats %>% 
  filter(focus == "gYoung_glYoung") %>% 
  select(-focus)

JGI_mature_de <- JGI_de_stats %>% 
  filter(focus == "gMature_glMature") %>% 
  select(-focus)


JGI_de_ncbi_annot <- read.delim("results/jgi_dea_sig_ncbi-annotated.txt")

str(JGI_de_ncbi_annot)
# 'data.frame':	15848 obs. of  3 variables:


JGI_de_genome_annot <- read.delim("results/jgi_dea_sig_genome-annotated.txt")

str(JGI_de_genome_annot)
# 'data.frame':	88069 obs. of  6 variables:


### Aggregating annotations
JGI_young_annots <- 
  left_join(JGI_de_stats, JGI_de_genome_annot, by = "ctg_cds")

JGI_young_annots <- 
  left_join(JGI_young_annots, gff3_mRNA, by = "mRNA")

JGI_young_annots <- JGI_young_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "start", "end",
  "domain_id", "domain_source", "domain_desc")

JGI_young_annots <- 
  left_join(JGI_young_annots, JGI_de_ncbi_annot, by = "ctg_cds")

write.table(JGI_young_annots, "results/jgi_young_dea_sig_aggregated-annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

###

JGI_mature_annots <- 
  left_join(JGI_mature_de, JGI_de_genome_annot, by = "ctg_cds")

JGI_mature_annots <- 
  left_join(JGI_mature_annots, gff3_mRNA, by = "mRNA")

JGI_mature_annots <- JGI_mature_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "start", "end",
  "domain_id", "domain_source", "domain_desc")

JGI_mature_annots <- 
  left_join(JGI_mature_annots, JGI_de_ncbi_annot, by = "ctg_cds")

str(JGI_mature_annots)
# 'data.frame':	2413 obs. of  12 variables:

write.table(JGI_mature_annots, "results/jgi_mature_dea_sig_aggregated-annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
