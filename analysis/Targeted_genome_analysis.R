library(dplyr)

library(stringr)



#####

tbb_blast_results <- 
  read.delim("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/UBC_TBB_orthologs.blastWRCv3Annotations.txt",
             stringsAsFactors = FALSE, header = FALSE) %>% 
  select(1,2)
colnames(tbb_blast_results) = c("ctg_cds", "mRNA")
tbb_blast_results$type <- "tbb"

tps_blast_results <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/putative_FL-TPS.UBC_assembly.blastWRCv3Annotations.txt",
             stringsAsFactors = FALSE, header = FALSE) %>% 
  select(1,2)
colnames(tps_blast_results) = c("ctg_cds", "mRNA")
tps_blast_results$type <- "tps"

p450_blast_results <- 
  read.delim("data/targeted-pathway-annotation/03-P450/putative_P450.UBC_assembly.blastWRCv3Annotations.txt",
             stringsAsFactors = FALSE, header = FALSE) %>% 
  select(1,2)
colnames(p450_blast_results) = c("ctg_cds", "mRNA")
p450_blast_results$type <- "p450"

sdr_blast_results <- 
  read.delim("data/targeted-pathway-annotation/04-SDR_ADH/putative_SDR.UBC_assembly.size_filtered.blastWRCv3Annotations.txt",
             stringsAsFactors = FALSE, header = FALSE) %>% 
  select(1,2)
colnames(sdr_blast_results) = c("ctg_cds", "mRNA")
sdr_blast_results$type <- "sdr"

akr_blast_results <- 
  read.delim("data/targeted-pathway-annotation/05-Aldo-keto-reductase/putative_AKR.UBC_assembly.blastWRCv3Annotations.txt",
             stringsAsFactors = FALSE, header = FALSE) %>% 
  select(1,2)
colnames(akr_blast_results) = c("ctg_cds", "mRNA")
akr_blast_results$type <- "akr"


##### TBP = Thujone Bionsynthesis Pathway
TBP_genome_orthologs <- 
  rbind(tbb_blast_results, tps_blast_results, p450_blast_results, 
        sdr_blast_results, akr_blast_results)


gff3_mRNA <- read.delim("data/gff3_condensed.txt")

TBP_genome_orthologs <- 
  left_join(TBP_genome_orthologs, gff3_mRNA, by = "mRNA")


##### Number of targeted genes on the same scaffold
table(duplicated(TBP_genome_orthologs$scaffold))
# FALSE  TRUE 
#   426  1449 

gt_one_DE <- TBP_genome_orthologs %>%
  select(scaffold, start, end) %>% 
  unique() %>% 
  group_by(scaffold) %>% 
  count(name = "count") %>% 
  filter(count > 1) %>%
  ungroup()

str(gt_one_DE)
# tibble [106 × 2] (S3: tbl_df/tbl/data.frame)

next_closest_DE <-
  TBP_genome_orthologs %>%
  filter(TBP_genome_orthologs$scaffold %in% gt_one_DE$scaffold) %>% 
  select(-ctg_cds) %>% 
  unique() %>% 
  group_by(scaffold) %>%
  arrange(start) %>% 
  mutate(closest_DE_locus = start-lag(start)) %>% 
  ungroup() %>% 
  arrange(scaffold)


# For illustrative purpose
next_closest_DE %>% filter(scaffold == 29381898)

next_closest_DE %>% filter(scaffold == 29382445)

next_closest_DE %>% filter(scaffold == 29377589)

next_closest_DE %>% filter(scaffold == 29380088)
