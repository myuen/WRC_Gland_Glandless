library(dplyr)
library(purrr)
library(stringr)



#Read genome annotations from GFF3 file
gff3 <- read.delim("data/Tplicatav3.1c.primaryTrs.gff3", header = FALSE,
                   comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

str(gff3)
# 'data.frame':	293312 obs. of  9 variables:


#Only keeping genes from GFF3 file
scaffold_genes <- gff3 %>%
  filter(V3 == "gene") %>%
  select(4,5,9)

colnames(scaffold_genes) <- 
  c("locus_start", "locus_end", "attributes")

str(scaffold_genes)
# 'data.frame':	39659 obs. of  3 variables:


# Extracting gene names and respective locations
scaffold_genes$locus <-
  map(str_split(scaffold_genes[,"attributes"], ";"), function(a){
    b = split(a[2], "=")
    c = str_split(b, "=")
    map(c, function(d) {
      return(d[2])
    })
  }) %>% unlist()

scaffold_genes <- scaffold_genes %>% 
  select(locus, locus_start, locus_end)


### Read DE statistics
UBC_de_stats <- read.delim("results/ubc_dea_sig_results.txt", 
                           stringsAsFactors = FALSE)

colnames(UBC_de_stats)[1] <- "ctg_cds"

str(UBC_de_stats)
# 'data.frame':	2399 obs. of  3 variables:


# Read annotations
### Read UBC DE CDS with NCBI annotations
UBC_de_ncbi_annot <- read.delim("results/ubc_dea_sig_ncbi-annotated.txt")

str(UBC_de_ncbi_annot)
# 'data.frame':	1164 obs. of  3 variables:


### Read UBC DE CDS with genome annotations
UBC_de_genome_annot <- read.delim("results/ubc_dea_sig_genome-annotated.txt",
                                  stringsAsFactors = FALSE)

str(UBC_de_genome_annot)
# 'data.frame':	9484 obs. of  6 variables:


### Aggregating annotations
UBC_annots <- 
  left_join(UBC_de_stats, UBC_de_genome_annot, by = "ctg_cds")

UBC_annots <- 
  left_join(UBC_annots, scaffold_genes, by = "locus")

UBC_annots <- UBC_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "locus_start", "locus_end",
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
UBC_annots <- 
  left_join(UBC_de_stats, UBC_de_genome_annot, by = "ctg_cds")

UBC_annots <- 
  left_join(UBC_annots, scaffold_genes, by = "locus")

UBC_annots <- 
  left_join(UBC_annots, UBC_de_ncbi_annot, by = "ctg_cds")



###
JGI_young_annots <- 
  left_join(JGI_young_de, JGI_de_genome_annot, by = "ctg_cds")

JGI_young_annots <- 
  left_join(JGI_young_annots, scaffold_genes, by = "locus")

JGI_young_annots <- JGI_young_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "locus_start", "locus_end",
  "domain_id", "domain_source", "domain_desc")

JGI_young_annots <- 
  left_join(JGI_young_annots, JGI_de_ncbi_annot, by = "ctg_cds")

str(JGI_young_annots)
# 'data.frame':	110240 obs. of  12 variables:

write.table(JGI_young_annots, "results/jgi_young_dea_sig_aggregated-annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

###

JGI_mature_annots <- 
  left_join(JGI_mature_de, JGI_de_genome_annot, by = "ctg_cds")

JGI_mature_annots <- 
  left_join(JGI_mature_annots, scaffold_genes, by = "locus")

JGI_mature_annots <- JGI_mature_annots %>% select(
  "ctg_cds", "logFC", "adj.P.Val",
  "scaffold", "locus", "locus_start", "locus_end",
  "domain_id", "domain_source", "domain_desc")

JGI_mature_annots <- 
  left_join(JGI_mature_annots, JGI_de_ncbi_annot, by = "ctg_cds")

str(JGI_mature_annots)
# 'data.frame':	2413 obs. of  12 variables:

write.table(JGI_mature_annots, "results/jgi_mature_dea_sig_aggregated-annotation.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
