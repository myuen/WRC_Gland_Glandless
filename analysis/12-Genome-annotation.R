library(dplyr)
library(purrr)
library(stringr)

#BLAST header
blast_header <- 
  c("ctg_cds", "mRNA", "evalue", "bitscore", "length", "pident", "ppos", 
    "qcovs", "qlen", "qstart", "qend", "qframe",
    "slen", "sstart", "send", "staxid", "sskingdom", "sscinames", "scomnames", "salltitles")

#####

# Read the statistics from UBC DEA
ubc_sigDE <- read.table("results/ubc_dea_sig_results.txt", 
                        header = TRUE, stringsAsFactors = FALSE)

colnames(ubc_sigDE)[1] <- "ctg_cds"

str(ubc_sigDE)


# Read BLAST results
ubc_genome_blast <- 
  read.delim("data/ubc_dea_sig.blastWRCv3Annotation.txt",
             stringsAsFactors = FALSE, sep = "\t", col.names = blast_header) %>% 
  select(-c("staxid", "sskingdom", "sscinames", "scomnames"))

str(ubc_genome_blast)
# 'data.frame':	13669 obs. of  16 variables:

ubc_genome_blast <- ubc_genome_blast %>%
  filter(ppos >= 95 & qcovs >= 50)

str(ubc_genome_blast)
# 'data.frame':	2150 obs. of  16 variables:


#Take the top best hit from BLAST
ubc_genome_blast <- ubc_genome_blast %>% 
  group_by(ctg_cds) %>% 
  arrange(length, ppos) %>% 
  slice_head(n = 1) %>% 
  ungroup()

ubc_genome_blast <- ubc_genome_blast %>%
  select(ctg_cds, mRNA)

ubc_genome_blast$gene <- 
  str_replace(ubc_genome_blast$mRNA, "\\.\\d$", "")

ubc_genome_blast$scaffold <- ubc_genome_blast$gene %>% 
  str_replace("Thupl.", "") %>% str_replace("s\\d+", "")

UBC_annots <- left_join(ubc_sigDE, ubc_genome_blast, by = "ctg_cds")

UBC_annots <- UBC_annots %>% 
  select(1,2,3,6,5,4)

table(is.na(UBC_annots$scaffold))
# FALSE  TRUE 
#  1610   789 

write.table(UBC_annots, "results/ubc_dea_sig_genome-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

#####

# Read the statistics from UBC DEA
jgi_sigDE <- read.table("results/jgi_dea_sig_results.txt", 
                        header = TRUE, stringsAsFactors = FALSE)

colnames(jgi_sigDE)[1] <- "ctg_cds"

str(jgi_sigDE)
# 'data.frame':	22404 obs. of  4 variables:


# Read BLAST results
jgi_genome_blast <- 
  read.delim("data/jgi_dea_sig.blastWRCv3Annotation.txt",
             stringsAsFactors = FALSE, sep = "\t", col.names = blast_header) %>% 
  select(-c("staxid", "sskingdom", "sscinames", "scomnames"))

str(jgi_genome_blast)
# 'data.frame':	103175 obs. of  16 variables:

jgi_genome_blast <- jgi_genome_blast %>%
  filter(ppos >= 95 & qcovs >= 50)

str(jgi_genome_blast)
# 'data.frame':	23479 obs. of  16 variables:


#Take the top best hit from BLAST
jgi_genome_blast <- jgi_genome_blast %>% 
  group_by(ctg_cds) %>% 
  arrange(length, ppos) %>% 
  slice_head(n = 1) %>% 
  ungroup()

jgi_genome_blast <- jgi_genome_blast %>%
  select(ctg_cds, mRNA)

jgi_genome_blast$gene <- 
  str_replace(jgi_genome_blast$mRNA, "\\.\\d$", "")

jgi_genome_blast$scaffold <- jgi_genome_blast$gene %>% 
  str_replace("Thupl.", "") %>% str_replace("s\\d+", "")

JGI_annots <- left_join(jgi_sigDE, jgi_genome_blast, by = "ctg_cds")

JGI_annots <- JGI_annots %>% 
  select(1,2,3,4,7,6,5)

table(JGI_annots$focus)
# gMature_glMature   gYoung_glYoung 
#              669            21735 

table(is.na(JGI_annots$scaffold))
# FALSE  TRUE 
#  1610   789 

write.table(JGI_annots, "results/jgi_dea_sig_genome-annotated.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
