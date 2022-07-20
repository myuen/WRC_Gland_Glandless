library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

#####

TPS_exp_plot <- function(stats) {

  legend_for_annotated <- 
    substitute(paste(italic('in silico'), ' predicted mono-TPS'^1))
  
  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 2.5, aes(shape = type)) +
    
    scale_x_continuous(limits = c(-bound, bound), 
                       breaks = c(-bound, -10, -5, -2, 0, 2, 5, 10, bound)) + 
    scale_y_discrete(labels = NULL) +

    scale_shape_manual(name = "",
                       labels = c("Predicted" = legend_for_annotated),
                       values = c(1, 5, 16, 17, 15)) +

    labs(caption =
           expression(paste("p-value <= 0.05\t", ""^1, "Shalev T et al. (2018)"))) +

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


# Read differential expressed statistics results
UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


### Read functionally characterized orthologs
UBC_cloned_TPSs <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/UBC_TPS_orthologs.cdsID.txt",
             header = FALSE)
             
colnames(UBC_cloned_TPSs) <- c("type", "accession", "cds")

UBC_cloned_TPSs <-
  UBC_cloned_TPSs %>% select(cds, type)

str(UBC_cloned_TPSs)
# 'data.frame':	4 obs. of  2 variables:


### Read in silico predicted mono-TPS in assembly
UBC_predicted_monoTPSs <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/WRC_monoTPS.blastp.UBC_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)

# Remove rows if predicted monoTPS already funcationally characterized
UBC_predicted_monoTPSs <- 
  UBC_predicted_monoTPSs[!UBC_predicted_monoTPSs$V2 %in% UBC_cloned_TPSs$cds,]

# Filter out low percent positive (V7) from BLASTp
UBC_predicted_monoTPSs <- UBC_predicted_monoTPSs %>% 
  filter(V7 >= 90) %>% 
  select(V2) %>% 
  unique()

UBC_predicted_monoTPSs$type <- "Predicted"

colnames(UBC_predicted_monoTPSs)[1] <- "cds"

str(UBC_predicted_monoTPSs)
# 'data.frame':	9 obs. of  2 variables:


# Concatenate all TPSs collected so far (cloned + predicted)
UBC_TPSs <- rbind(UBC_cloned_TPSs, UBC_predicted_monoTPSs)


# Read putative TPS CDS IDs from functional domain scan with hmmscan and rpsblast.
UBC_putative_TPSs <- 
  read.csv("data/targeted-pathway-annotation/02-TPS/putative_FL-TPS.UBC_assembly.cdsID.txt",
           header = FALSE, col.names = c("cds")) %>% mutate("type" = "Putative")

UBC_putative_TPSs <- 
  UBC_putative_TPSs[!UBC_putative_TPSs$cds %in% UBC_TPSs$cds,]

str(UBC_putative_TPSs)
# 'data.frame':	73 obs. of  2 variables:

UBC_TPSs <- rbind(UBC_TPSs, UBC_putative_TPSs)

dim(UBC_TPSs)
# [1] 86  2

UBC_TPSs$type <- factor(UBC_TPSs$type, 
                        levels = c("Putative", "Predicted", 
                                  "Diterpene_synthase", "Sabinene_synthase",
                                  "Alpha_pinene_synthase"))

UBC_TPS_stats <- inner_join(UBC_TPSs, UBC_dea_stats, by = c("cds"))


UBC_TPS_exp_plot <-
  TPS_exp_plot(UBC_TPS_stats)

(UBC_TPS_exp_plot <- 
  UBC_TPS_exp_plot + 
  ggtitle("UBC Putative Full Length Terpene Synthase Expression Plot",
          subtitle = "Size filtered >= 500 a.a"))

ggsave("results/figures/UBC_TPS_expression.svg", plot = UBC_TPS_exp_plot)

######
# 
# JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
#                             header = TRUE, stringsAsFactors = FALSE)
# str(JGI_dea_stats)
# # 'data.frame':	84192 obs. of  4 variables:
# 
# 
# # Read putative TPS CDS IDs from functional domain scan with hmmscan and rpsblast.
# JGI_putative_TPSs <- 
#   scan("data/targeted-pathway-annotation/02-TPS/putative_FL-TPS.JGI_assembly.cdsID.txt",
#        what = "character")
# # Read 62 items
# 
# JGI_TPSs <- data.frame(row.names = c(unique(JGI_putative_TPSs)))
# 
# JGI_TPSs$type <- "Putative"
# 
# str(JGI_TPSs)
# # 'data.frame':	62 obs. of  1 variable:
# 
# 
# JGI_annotated_TPSs <- 
#   read.delim("data/targeted-pathway-annotation/02-TPS/WRC_TPS.blastp.JGI_assembly.txt",
#              stringsAsFactors = FALSE, header = FALSE)
# 
# #Filter out low percent positive (V7) from BLASTp
# JGI_annotated_TPSs <- JGI_annotated_TPSs %>% 
#   filter(V7 >= 90) %>% 
#   select(V2) %>% unique()
# 
# str(JGI_annotated_TPSs)
# # 'data.frame':	11 obs. of  1 variable:
# 
# JGI_TPSs[JGI_annotated_TPSs$V2, "type"] <- "Annotated"
# 
# 
# JGI_sabinene_syn <- 
#   read.delim("data/targeted-pathway-annotation/02-TPS/Sabinene_synthase.blastp.JGI_assembly.txt",
#              stringsAsFactors = FALSE, header = FALSE)
# 
# JGI_sabinene_syn <- JGI_sabinene_syn$V2
# 
# JGI_TPSs[JGI_sabinene_syn, "type"] <- "Sabinene Synthase"
# 
# JGI_TPSs <- JGI_TPSs %>% rownames_to_column("cds")
# 
# str(JGI_TPSs)
# # 'data.frame':	65 obs. of  2 variables:
# 
# JGI_TPS_stats <- inner_join(JGI_dea_stats, JGI_TPSs)
# 
# str(JGI_TPS_stats)  
# # 'data.frame':	348 obs. of  5 variables:
# 
# JGI_TPS_stats$type <- factor(JGI_TPS_stats$type,
#                              levels = c("Putative", "Annotated", "Sabinene Synthase"))
# 
# JGI_TPS_stats$focus <- factor(JGI_TPS_stats$focus, 
#                               levels = c("gYoung_glYoung", "gMature_glMature"))
# 
# JGI_TPS_exp_plot <- TPS_exp_plot(JGI_TPS_stats)
# 
# (JGI_TPS_exp_plot <- JGI_TPS_exp_plot + 
#   ggtitle("JGI Putative Full Length Terpene Synthase Expression Plot",
#           subtitle = "Size filtered >= 500 a.a") + 
#   facet_wrap(~focus, nrow = 2, strip.position = "right"))
# 
# ggsave("results/figures/JGI_TPS_expression.svg", plot = JGI_TPS_exp_plot)
# 
