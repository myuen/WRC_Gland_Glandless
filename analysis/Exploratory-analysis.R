library(dplyr)
library(ggplot2)


UBC_annots <- read.delim("results/ubc_dea_sig_aggregated-annotation.txt")

# How many number of scaffolds contain DE contigs?
(num_scaffolds_with_DE_ctgs <- UBC_annots %>% 
  filter(!is.na(scaffold)) %>% 
  select(scaffold) %>% 
  unique())

(UBC_annots %>% 
    filter(!is.na(scaffold))) %>% 
  select(scaffold) %>% 
  unique() %>% 
  count()

#     n
# 1 909


# How many DE locus on ***a single*** scaffolds?
(num_DE_locus_per_scaffold <- UBC_annots %>% 
  filter(!is.na(scaffold)) %>% 
  group_by(scaffold) %>% 
  select(scaffold, locus) %>%
  unique() %>% 
  count() %>% 
  ungroup())

str(num_DE_locus_per_scaffold)
# tibble [909 × 2] (S3: tbl_df/tbl/data.frame)

colnames(num_DE_locus_per_scaffold)[2] <- "count"

# Minimum # of DE contig on ***a single*** scaffold
min(num_DE_locus_per_scaffold$count)
# [1] 1

# Maximum # of DE contig on ***a single*** scaffold
max(num_DE_locus_per_scaffold$count)
# [1] 8


# Plot for counts of # of DE contigs per scaffold 
count_num_DE_locus <- num_DE_locus_per_scaffold %>% 
  group_by(count) %>% 
  count() %>% 
  ungroup()

#   count     n
# 1     1   709
# 2     2   144
# 3     3    39
# 4     4    11
# 5     5     4
# 6     6     1
# 7     8     1

max_x <- max(count_num_DE_locus$count)
max_y <- max(count_num_DE_locus$n)

ggplot(count_num_DE_locus, aes(x = count, y = n)) + 
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous("Number of DE contigs",
                     breaks = c(1:max_x), labels = c(1:max_x)) + 
  scale_y_continuous("Number of scaffolds", breaks = seq(0, max_y, 100)) +
  ggtitle("Exploring the Genome",
    subtitle = "Counts of number of differential expressed locus/loci per scaffolds") +
  theme_bw()


# For scaffold with more than 1 DE locus, calculate the closest distance to 
# the next DE locus
gt_one_DE_locus <- num_DE_locus_per_scaffold %>% 
  filter(count > 1) %>% 
  select(scaffold)

gt_one_DE_locus <- as.integer(gt_one_DE_locus$scaffold)


next_closest_DE_by_locus <-
  UBC_annots %>%
  filter(scaffold %in% gt_one_DE_locus) %>% 
  select(scaffold, locus, locus_start, locus_end) %>% 
  unique() %>% 
  group_by(scaffold) %>%
  arrange(locus_start) %>% 
  mutate(closest_DE_locus = lead(locus_start)-locus_start) %>% 
  filter(!is.na(closest_DE_locus)) %>%
  ungroup()


# Plot next closest DE locus
min_x <- ceiling(min(next_closest_DE_by_locus$closest_DE_locus)/1000)
max_x <- max(next_closest_DE_by_locus$closest_DE_locus)/1000

ggplot(next_closest_DE_by_locus, aes(closest_DE_locus/1000)) + 
  geom_density() +
  ggtitle("Exploring the Genome", 
          subtitle = paste0("Next closest differentially expressed locus ",
          "on the same scaffold")) + 
  scale_x_continuous("Distance (kbp)",
                     breaks = c(min_x, seq(1000, max_x, 1000)),
                     labels = c(min_x, seq(1000, max_x, 1000))) +
  scale_y_continuous("", labels = NULL) +
  theme_bw()


# Most happened PFAM domain
UBC_annots$description <- toupper(UBC_annots$description)

pfam_count <- UBC_annots %>% 
  filter(source == "PFAM") %>% 
  select(ID, description) %>% 
  group_by(ID, description) %>% 
  count() %>% 
  ungroup()

top20_pfam <- slice_max(pfam_count, order_by = n, n = 20)
