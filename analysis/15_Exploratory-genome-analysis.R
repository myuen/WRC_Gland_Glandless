library(dplyr)
library(ggplot2)
library(gt)
library(patchwork)


UBC_annots <- 
  read.delim("results/ubc_dea_sig_aggregated-annotation.txt")


# Note scaffold is read as character instead of integer.  Some scaffold
# name contain character '-'
UBC_annots$scaffold <- as.character(UBC_annots$scaffold)


# 1. Number of loci with higher expression in gland
UBC_annots %>% 
  filter(logFC >= 0) %>% 
  select(mRNA) %>% 
  unique() %>% count()
# 1,726 up-regulated contigs matched 955 genome loci


# 2. Number of loci with higher expression in glandless
UBC_annots %>% 
  filter(logFC <= 0) %>% 
  select(mRNA) %>% 
  unique() %>% count()
# 673 up-regulated contigs matched 364 genome loci


scaffold_length <- 
  read.delim("data/genome_scaffold_length.tsv", header = FALSE, 
             col.names = c("scaffold", "length"), 
             stringsAsFactors = FALSE)


# 3. How many number of scaffolds contain DE loci?
(UBC_annots %>% 
    filter(!is.na(scaffold)) %>% 
    select(scaffold) %>% 
    unique() %>% 
    count())
# n = 909


# Enumerate number of DE locus *** per scaffold ***?
(num_DE_locus_by_scaffold <- UBC_annots %>% 
  filter(!is.na(scaffold)) %>% 
  group_by(scaffold) %>% 
  select(scaffold, locus) %>%
  unique() %>% 
  count(name = "num_DE_locus") %>% 
  ungroup()
)


# Add scaffold length 
num_DE_locus_by_scaffold <- 
  left_join(num_DE_locus_by_scaffold, scaffold_length, by = "scaffold")


# 4a. Minimum # of DE loci on ***a single*** scaffold
min(num_DE_locus_by_scaffold$num_DE_locus)
# n = 1

# 4b. Maximum # of DE loci on ***a single*** scaffold
max(num_DE_locus_by_scaffold$num_DE_locus)
# n = 8


# Tabulate counts of scaffolds by num of DE locus
(count_of_scaffold_by_DE_locus <- 
    num_DE_locus_by_scaffold %>% 
    group_by(num_DE_locus) %>% 
    count(name = "num_scaffolds") %>% 
    ungroup()
)
#   num_DE_locus num_scaffolds
# 1            1           709
# 2            2           144
# 3            3            39
# 4            4            11
# 5            5             4
# 6            6             1
# 7            8             1

max_x <- max(count_of_scaffold_by_DE_locus$num_DE_locus)
max_y <- max(count_of_scaffold_by_DE_locus$num_scaffolds)

# labs(caption = expression(paste0("Differential expression = >caption here", beta))) +


(p <- ggplot(count_of_scaffold_by_DE_locus, aes(x = num_DE_locus, y = num_scaffolds)) + 
    geom_point(size = 2) +
    
    ggtitle("Exploring the Genome",
            subtitle = "Number of DE locus per scaffolds") +
    
    labs(caption = "|logFC| \u2265 2,  Adj. p-value \u2264 0.05")+
    
    scale_x_continuous("Number of DE locus",
                       breaks = c(1:max_x), labels = c(1:max_x)) +
    
    scale_y_continuous("Number of scaffolds", breaks = seq(1, max_y, 100)) +
    
    theme_bw() +

    theme(panel.grid = element_line(colour = "grey97"))
)


(q <- ggplot(num_DE_locus_by_scaffold, 
             aes(x = num_DE_locus, y = length, group = num_DE_locus)) + 
    geom_jitter(colour = "grey60", alpha = 0.75, width = 0.2) + 
    geom_boxplot(colour = "black", fill = "white", width = 0.35, alpha = 0.5) +
    
    ggtitle("",
            subtitle = "Scaffold Length") +

    scale_x_continuous("Number of DE locus",
                       breaks = c(1:max_x), labels = c(1:max_x)) +
    ylab("Scaffold Length") +
    
    theme_bw() +
    
    theme(panel.grid = element_line(colour = "grey97"))
)

p + q

ggsave("results/figures/DE_locus_with_length.svg", p+q)

ggsave("results/figures/DE_locus_with_length.jpg", p+q)


#####


# 5.For scaffold with more than 1 DE locus, calculate the closest 
# distance to the next DE locus

# Get scaffold ID for those with more than 1 DE locus
gt_one_DE_locus <- num_DE_locus_by_scaffold %>% 
  filter(num_DE_locus > 1) %>% 
  select(scaffold)

gt_one_DE_locus <- as.character(gt_one_DE_locus$scaffold)


# Find the next closest DE locus
next_closest_DE_by_locus <-
  UBC_annots %>%
  filter(scaffold %in% gt_one_DE_locus) %>%
  select(scaffold, locus, start, end) %>%
  unique() %>%
  group_by(scaffold) %>%
  arrange(start) %>%
  mutate(closest_DE_locus = start-lag(start)) %>%
  ungroup() %>%
  arrange(scaffold)


(min_x <- next_closest_DE_by_locus %>% 
  filter(closest_DE_locus > 0) %>% 
  select(closest_DE_locus) %>% 
  min() / 1000
)
# [1] 2.996

(max_x <- next_closest_DE_by_locus %>% 
  filter(closest_DE_locus > 0) %>% 
  select(closest_DE_locus) %>% 
  max() / 1000
) 
# [1] 8088.799


# Plot next closest DE locus
(g <- ggplot(next_closest_DE_by_locus, aes(closest_DE_locus/1000)) + 
    geom_density() +
    ggtitle("Exploring the Genome", 
            subtitle = expression(paste("Next closest differentially ",
                                        "expressed locus on the ", bold("same scaffold")))) +
    scale_x_continuous("Distance (kbp)",
                       breaks = c(min_x, seq(1000, max_x, 1000)),
                       labels = c(min_x, seq(1000, max_x, 1000))) +
    scale_y_continuous("", labels = NULL) +
    
    theme_bw() +
    
    theme(panel.grid = element_line(colour = "grey97"))
)

ggsave("results/figures/Next_closest_DE.svg", g)
ggsave("results/figures/Next_closest_DE.jpg", g)


# 6. Find the next closest DE by gene of interest (goi)
tbb <- 
  read.delim("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/UBC_TBB_orthologs.txt")

tbb <- tbb %>% mutate(type = "tbb") %>% 
  select(cds, type)

tps <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/putative_TPS.UBC_assembly.cdsID.txt",
             col.names = c("cds")) %>% 
  mutate(type = "tps")

p450 <- 
  read.delim("data/targeted-pathway-annotation/03-P450/putative_P450.UBC_assembly.cdsID.txt", 
             col.names = c("cds")) %>% 
  mutate(type = "p450")

sdr <- 
  read.delim("data/targeted-pathway-annotation/04-SDR_ADH/putative_SDR.UBC_assembly.cdsID.txt",
             col.names = c("cds")) %>%
  mutate(type = "sdr")

akr <- 
  read.delim("data/targeted-pathway-annotation/05-Aldo-keto-reductase/putative_AKR.UBC_assembly.cdsID.txt",
             col.names = c("cds")) %>%
  mutate(type = "akr")

goi <- rbind(tbb, tps,p450, sdr, akr) %>% 
  unique()
# 'data.frame':	2229 obs. of  2 variables:

next_closest_DE_by_goi <-
  inner_join(UBC_annots, goi, by = c("ctg_cds" = "cds")) %>% 
  filter(scaffold %in% gt_one_DE_locus) %>% 
  select(scaffold, locus, start, end, type) %>%
  unique() %>%
  group_by(scaffold) %>%
  arrange(start) %>%
  mutate(closest_DE_locus = start-lag(start)) %>%
  ungroup() %>%
  arrange(scaffold)

(min_x <- next_closest_DE_by_goi %>% 
    filter(closest_DE_locus > 0) %>% 
    select(closest_DE_locus) %>% 
    min() / 1000
)
# [1] 8.218

(max_x <- next_closest_DE_by_goi %>% 
    filter(closest_DE_locus > 0) %>% 
    select(closest_DE_locus) %>% 
    max() / 1000
) 
# [1] 4068.358

(goi_plot <- ggplot(next_closest_DE_by_goi, aes(closest_DE_locus/1000)) + 
    geom_density() +
    ggtitle("Exploring the Genome", 
            subtitle = expression(paste("Next closest differentially ",
                                        "expressed locus by ", bold("gene of interests "),
                                        "on the ", bold("same scaffold")))) +
    scale_x_continuous("Distance (kbp)",
                       breaks = c(min_x, seq(1000, max_x, 1000)),
                       labels = c(min_x, seq(1000, max_x, 1000))) +
    scale_y_continuous("", labels = NULL) +
    
    theme_bw() +
    
    theme(panel.grid = element_line(colour = "grey97"))
)

ggsave("results/figures/Next_closest_DE_goi.svg", goi_plot)
ggsave("results/figures/Next_closest_DE_goi.jpg", goi_plot)


# For illustration
next_closest_DE_by_goi$type <- toupper(next_closest_DE_by_goi$type)

for_illustration <- next_closest_DE_by_goi %>%
  group_by(scaffold) %>% 
  arrange(desc(closest_DE_locus), .by_group = TRUE) %>% 
  count(name = "count") %>% 
  filter(count > 1) %>% 
  ungroup()

for_illustration <- next_closest_DE_by_goi %>% 
  filter(scaffold %in% for_illustration$scaffold)

colnames(for_illustration) <- 
  c("Scaffold", "Locus", "Locus Start", "Locus Stop",
    "Class", "Next Closest Locus")

tbl_for_illustration <- for_illustration %>% 
  gt(rowname_col = "Locus",
     groupname_col = "Scaffold") %>% 
  fmt_integer(c("Locus Start", "Locus Stop", "Next Closest Locus")) %>% 

  tab_style(style = list(cell_text(weight = "bold")), 
            locations = cells_row_groups()) %>% 

  tab_header(title = "Next Closest DE by Gene-of-Interest (GOI)",
             subtitle = "Grouped by Scaffold") %>% 
  tab_style(style = list(cell_text(weight = "bold"),
                         cell_fill()),
          locations = list(cells_column_labels(), cells_title()))

gtsave(tbl_for_illustration, "
       results/undated_lab_meeting/tbl_for_illustration.html", inline_css = TRUE)



# Capitalized descriptions to remove redundancy
UBC_annots$domain_desc <- toupper(UBC_annots$domain_desc)


# Most abundant PFAM domains by count
pfam_count <- UBC_annots %>% 
  filter(domain_source == "PFAM") %>% 
  select(domain_id, domain_desc) %>% 
  count(domain_id, domain_desc, name = "count", sort = TRUE)

(top20_pfam <- pfam_count %>% slice_max(count, n = 20))

tbl_for_top20_pfam <-
  top20_pfam %>% 
    rename("Domain ID" = domain_id, 
           "Domain Description" = domain_desc,
           "Count" = count) %>% 
  gt() %>% 
  fmt_integer(Count) %>% 
  tab_header(title = "Top 20 Pfam domain by count") %>% 
  tab_style(style = list(cell_text(weight = "bold"),
                         cell_fill()),
            locations = list(cells_column_labels(), cells_title()))

gtsave(tbl_for_top20_pfam, 
       "results/undated_lab_meeting/tbl_for_top20_pfam.html", inline_css = TRUE)

top20_pfam_stats <-
  UBC_annots %>% 
  filter(domain_id %in% top20_pfam$domain_id)

(top20_pfam_stats_plot <- ggplot(top20_pfam_stats, aes(x = logFC, y = domain_id)) + 
    geom_boxplot() +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    ggtitle("Top 20 most abundant PFAM domain with gene expression") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("Domain Description") +
    scale_y_discrete(
      labels = c(
        "PF00201" = "UDP-GLUCORONOSYL AND UDP-GLUCOSYL TRANSFERASE",
        "PF00560" = "LEUCINE RICH REPEAT",
        "PF01397" = "TERPENE SYNTHASE, N-TERMINAL DOMAIN",
        "PF03936" = "TERPENE SYNTHASE FAMILY, METAL BINDING DOMAIN",
        "PF07693" = "KAP FAMILY P-LOOP DOMAIN",
        "PF07714" = "PROTEIN TYROSINE KINASE")) +
    theme_bw() + 
    
    theme(panel.grid = element_line(colour = "grey97"))
  )

ggsave("results/figures/top20_pfam_exp_plot.svg", top20_pfam_stats_plot)

ggsave("results/figures/top20_pfam_exp_plot.jpg", top20_pfam_stats_plot)
