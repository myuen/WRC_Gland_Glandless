library(dplyr)
library(purrr)
library(stringr)


#Read annotation GFF3
gff3 <- read.delim("data/Tplicatav3.1c.primaryTrs.gff3", header = FALSE,
                   comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

str(gff3)
# 'data.frame':	293312 obs. of  9 variables:


# Only keeping genes
scaffold_genes <- gff3 %>%
  filter(V3 == "gene") %>%
  select(4,5,9)

colnames(scaffold_genes) <- 
  c("gene_start", "gene_end", "attributes")

str(scaffold_genes)
# 'data.frame':	39659 obs. of  3 variables:

scaffold_genes$gene <-
  map(str_split(scaffold_genes[,"attributes"], ";"), function(a){
    b = split(a[2], "=")
    c = str_split(b, "=")
    map(c, function(d) {
      return(d[2])
    })
  }) %>% unlist()

scaffold_genes <- scaffold_genes %>% 
  select(gene, gene_start, gene_end)

#####
# Read UBC annotated DE contigs CDS annotations
UBC_de_ncbi_annot <- read.delim("results/ubc_dea_sig_ncbi-annotated.txt")

str(UBC_de_ncbi_annot)
# 'data.frame':	2399 obs. of  5 variables:

UBC_de_genome_annot <- read.delim("results/ubc_dea_sig_genome-annotated.txt")

UBC_de_genome_annot <- UBC_de_genome_annot %>% select(-c(logFC, adj.P.Val))

str(UBC_de_genome_annot)
# 'data.frame':	2399 obs. of  4 variables:

UBC_annots <- left_join(UBC_de_ncbi_annot, UBC_de_genome_annot, by = "ctg_cds")

UBC_annots <- left_join(UBC_annots, scaffold_genes, by = "gene")

next_closest_DE_ctg <- UBC_annots %>%
  select(scaffold, gene, gene_start, gene_end) %>% 
  unique() %>% 
  group_by(scaffold) %>%
  arrange(gene_start) %>% 
  mutate(closest_DE_gene = lead(gene_start) - gene_start) %>% 
  ungroup() %>% 
  select(gene, closest_DE_gene)

next_closest_DE_ctg$closest_DE_gene <- 
  abs(next_closest_DE_ctg$closest_DE_gene)

UBC_annots <- 
  left_join(UBC_annots, next_closest_DE_ctg, by = "gene")

DE_count_by_scaffold <- UBC_annots[!is.na(UBC_annots$scaffold), ] %>% 
  group_by(scaffold) %>% 
  count()

max(DE_count_by_scaffold$n)
# [1] 10

###

JGI_DE_CDS_on_scaffold <- 
  read.delim("data/jgi_dea_sig.blastWRCv3Annotation.txt",
             stringsAsFactors = FALSE, sep = "\t", col.names = blast_header) %>% 
  select(-c("staxid", "sskingdom", "sscinames", "scomnames"))

str(JGI_DE_CDS_on_scaffold)
# 'data.frame':	103175 obs. of  16 variables:
