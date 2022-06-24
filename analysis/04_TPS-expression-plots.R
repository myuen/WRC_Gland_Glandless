library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

#####

TPS_exp_plot <- function(stats, title) {

  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
  geom_point(size = 2.5, aes(shape = type)) +
  scale_x_continuous(limits = c(-bound, bound), 
                     breaks = c(-bound, -10, -5, -2, 0, 2, 5, 10, bound)) + 
  scale_y_discrete(labels = NULL) +
  scale_shape_manual(name = "",
                     # labels = c("Putative", "Annotated", "Sabinene Synthase"),
                     # labels = character(stats$type),
                     values = c(1, 19, 15)) +
  
  ggtitle(title) +
  labs(caption = "p-value <= 0.05") +
  xlab(expression(paste("log"[2],"FC"))) +
  ylab("") +
  
  geom_vline(xintercept = c(-2, 2), colour = "black", linetype = 2) +
  
  theme(legend.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, "black"), 
        panel.grid = element_line(colour = "grey95"),
        axis.ticks = element_line(size = 0))
}

######

UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


# Read putative TPS CDS IDs from functional domain scan with hmmscan and rpsblast.
UBC_putative_TPSs <- 
  scan("data/targeted-pathway-annotation/02-TPS/putative_FL-TPS.UBC_assembly.cdsID.txt",
              what = "character")
# Read 85 items


# Construct dataframe.  Easier to remove duplicates with CDS ID as rownames
UBC_TPSs <- data.frame(row.names = c(unique(UBC_putative_TPSs)))

# Default annotation as "Putative"
UBC_TPSs$type <- "Putative"

str(UBC_TPSs)
# 'data.frame':	85 obs. of  1 variable:


# Read annotated TPS in assembly
UBC_annotated_TPSs <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/WRC_TPS.blastp.UBC_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)


# Filter out low percent positive (V7) from BLASTp
UBC_annotated_TPSs <- UBC_annotated_TPSs %>% 
  filter(V7 >= 90) %>% 
  select(V2) %>% unique()

str(UBC_annotated_TPSs)
# 'data.frame':	12 obs. of  1 variable:

UBC_TPSs[UBC_annotated_TPSs$V2, "type"] <- "Annotated"


# Catalogue Sabinene Synthase in dataframe
UBC_sabinene_syn <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/Sabinene_synthase.blastp.UBC_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)

UBC_sabinene_syn <- UBC_sabinene_syn$V2

UBC_TPSs[UBC_sabinene_syn, "type"] <- "Sabinene Synthase"


# Turn rowname to column
UBC_TPSs <- UBC_TPSs %>% rownames_to_column("cds")

str(UBC_TPSs)
# 'data.frame':	304 obs. of  2 variables:


UBC_TPS_stats <- inner_join(UBC_dea_stats, UBC_TPSs)

str(UBC_TPS_stats)  
# 'data.frame':	107 obs. of  4 variables:

UBC_TPS_stats$type <- factor(UBC_TPS_stats$type,
                             levels = c("Putative", "Annotated", "Sabinene Synthase"))

(UBC_TPS_exp_plot <- 
  TPS_exp_plot(UBC_TPS_stats, "UBC Full Length Terpene Synthase Expression Plot"))

ggsave("results/figures/UBC_TPS_expression.svg", plot = UBC_TPS_exp_plot)

######

JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(JGI_dea_stats)
# 'data.frame':	84192 obs. of  4 variables:


# Read putative TPS CDS IDs from functional domain scan with hmmscan and rpsblast.
JGI_putative_TPSs <- 
  scan("data/targeted-pathway-annotation/02-TPS/putative_FL-TPS.JGI_assembly.cdsID.txt",
       what = "character")
# Read 62 items

JGI_TPSs <- data.frame(row.names = c(unique(JGI_putative_TPSs)))

JGI_TPSs$type <- "Putative"

str(JGI_TPSs)
# 'data.frame':	62 obs. of  1 variable:


JGI_annotated_TPSs <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/WRC_TPS.blastp.JGI_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)

#Filter out low percent positive (V7) from BLASTp
JGI_annotated_TPSs <- JGI_annotated_TPSs %>% 
  filter(V7 >= 90) %>% 
  select(V2) %>% unique()

str(JGI_annotated_TPSs)
# 'data.frame':	11 obs. of  1 variable:

JGI_TPSs[JGI_annotated_TPSs$V2, "type"] <- "Annotated"


JGI_sabinene_syn <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/Sabinene_synthase.blastp.JGI_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)

JGI_sabinene_syn <- JGI_sabinene_syn$V2

JGI_TPSs[JGI_sabinene_syn, "type"] <- "Sabinene Synthase"

JGI_TPSs <- JGI_TPSs %>% rownames_to_column("cds")

str(JGI_TPSs)
# 'data.frame':	65 obs. of  2 variables:

JGI_TPS_stats <- inner_join(JGI_dea_stats, JGI_TPSs)

str(JGI_TPS_stats)  
# 'data.frame':	348 obs. of  5 variables:

JGI_TPS_stats$type <- factor(JGI_TPS_stats$type,
                             levels = c("Putative", "Annotated", "Sabinene Synthase"))

JGI_TPS_stats$focus <- factor(JGI_TPS_stats$focus, 
                              levels = c("gYoung_glYoung", "gMature_glMature"))

JGI_TPS_exp_plot <-
  TPS_exp_plot(JGI_TPS_stats, "JGI Full Length Terpene Synthase Expression Plot")

(JGI_TPS_exp_plot <- 
    JGI_TPS_exp_plot + facet_wrap(~focus, nrow = 2, strip.position = "right"))

ggsave("results/figures/JGI_TPS_expression.svg", plot = JGI_TPS_exp_plot)
