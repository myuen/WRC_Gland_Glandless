library(dplyr)
library(purrr)
library(stringr)


#Read gene function
gene_function <- 
  read.delim("data/gene.functions.txt", comment.char = "#",
             stringsAsFactors = FALSE,
             col.names = c("locus", "domain_id", "domain_source", "domain_desc"))

#Remove emtpy function description
gene_function <- 
  gene_function %>% filter(!description == "")

str(gene_function)
# 'data.frame':	132271 obs. of  4 variables:
  
#####
  
#BLAST header
blast_header <- 
  c("ctg_cds", "mRNA", "evalue", "bitscore", "length", "pident", "ppos", 
    "qcovs", "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "staxid", "sskingdom", "sscinames", "scomnames", "salltitles")

#####

# Read BLAST results
ubc_genome_blast <- 
  read.delim("data/ubc_dea_sig.blastWRCv3Annotation.txt", sep = "\t",
             stringsAsFactors = FALSE, col.names = blast_header) %>% 
  select(ctg_cds, mRNA, ppos, qcovs, salltitles)

str(ubc_genome_blast)
# 'data.frame':	13669 obs. of  16 variables:


ubc_genome_blast <- ubc_genome_blast %>%
  filter(ppos >= 95 & qcovs >= 50)

str(ubc_genome_blast)
# 'data.frame':	2150 obs. of  5 variables:


#Take the top best hit from BLAST
ubc_top_genome_blast <- ubc_genome_blast %>% 
  group_by(ctg_cds) %>% 
  arrange(ppos, qcovs) %>% 
  slice_head(n = 1) %>% 
  ungroup()

ubc_top_genome_blast$locus <- 
  str_replace(ubc_top_genome_blast$mRNA, "\\.\\d$", "")

ubc_top_genome_blast$scaffold <- ubc_top_genome_blast$locus %>% 
  str_replace("Thupl.", "") %>% str_replace("s\\d+", "")

ubc_top_genome_blast <- ubc_top_genome_blast %>%
  select(ctg_cds, scaffold, locus)

#Add gene function
UBC_annots <- left_join(ubc_top_genome_blast, gene_function, by = "locus")

write.table(UBC_annots, "results/ubc_dea_sig_genome-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#####

# Read BLAST results
jgi_genome_blast <- 
  read.delim("data/jgi_dea_sig.blastWRCv3Annotation.txt", sep = "\t", 
             stringsAsFactors = FALSE, col.names = blast_header) %>% 
  select(ctg_cds, mRNA, ppos, qcovs, salltitles)

str(jgi_genome_blast)
# 'data.frame':	103175 obs. of  16 variables:

jgi_genome_blast <- jgi_genome_blast %>%
  filter(ppos >= 95 & qcovs >= 50)

str(jgi_genome_blast)
# 'data.frame':	23479 obs. of  5 variables:


#Take the top best hit from BLAST
jgi_top_genome_blast <- jgi_genome_blast %>% 
  group_by(ctg_cds) %>% 
  arrange(ppos, qcovs) %>% 
  slice_head(n = 1) %>% 
  ungroup()


jgi_top_genome_blast$locus <- 
  str_replace(jgi_top_genome_blast$mRNA, "\\.\\d$", "")

jgi_top_genome_blast$scaffold <- jgi_top_genome_blast$locus %>%
  str_replace("Thupl.", "") %>% str_replace("s\\d+", "")

jgi_top_genome_blast <- jgi_top_genome_blast %>%
  select(ctg_cds, scaffold, locus)

#Add gene function
JGI_annots <- 
  left_join(jgi_top_genome_blast, gene_function, by = "locus")

write.table(JGI_annots, "results/jgi_dea_sig_genome-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
