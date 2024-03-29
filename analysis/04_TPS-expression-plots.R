library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

#####

TPS_exp_plot <- function(stats) {

  legend_for_annotated <-
    substitute(paste(italic('in silico'), 'predicted mono-TPS'^1, sep = ''))

  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 4, aes(shape = type, colour = type)) +
    
    scale_x_continuous(limits = c(-bound, bound), 
                       breaks = c(-bound, -10, -5, -2, 0, 2, 5, 10, bound)) + 
    scale_y_discrete(labels = NULL) +

    scale_shape_manual(
      name = "",
      labels = c("Putative" = "Putative TPS",
                 "Predicted" = legend_for_annotated,
                 "Diterpene_synthase" = "Diterpene Synthase",
                 "Sabinene_synthase" = "Sabinene Synthase",
                 "Alpha_pinene_synthase" = "Alpha Pinene Synthase"),
      values = c(5, 1, 16, 17, 15)) +
    
    scale_color_manual(
      # guide = NULL,
      name = "",
      labels = c("Putative" = "Putative TPS",
                 "Predicted" = legend_for_annotated,
                 "Diterpene_synthase" = "Diterpene Synthase",
                 "Sabinene_synthase" = "Sabinene Synthase",
                 "Alpha_pinene_synthase" = "Alpha Pinene Synthase"),
      values = c("black", "black", "black", "red", "black")) +

    labs(caption =
           expression(paste("p-value <= 0.05 ", "\t"^1, "Shalev T et al. (2018)"))) +

    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +

    geom_vline(xintercept = c(-2, 2), colour = "red", linetype = 2) +

    theme(
      title = element_text(size = 20),
      
      legend.background = element_blank(),
      legend.text = element_text(size = 20),
      
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, "black"),
      panel.grid = element_line(colour = "grey98"),
      
      axis.text = element_text(size = 14),
      axis.ticks = element_line(size = 0))
}


######


# Read differential expressed statistics results
UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)


### Read functionally characterized orthologs
UBC_cloned_TPSs <- 
  read.delim("data/targeted-pathway-annotation/02-TPS/UBC_TPS_orthologs.cdsID.txt",
             header = FALSE)
             
colnames(UBC_cloned_TPSs) <- c("type", "accession", "cds")

UBC_cloned_TPSs <-
  UBC_cloned_TPSs %>% select(cds, type)

str(UBC_cloned_TPSs)
# 'data.frame':	5 obs. of  2 variables:


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
# 'data.frame':	8 obs. of  2 variables:


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

UBC_TPS_stats <- UBC_TPS_stats %>% 
  filter(adj.P.Val < 0.05)
# 'data.frame':	16 obs. of  4 variables:

UBC_TPS_exp_plot <-
  TPS_exp_plot(UBC_TPS_stats)

(UBC_TPS_exp_plot <- 
  UBC_TPS_exp_plot + 
  ggtitle("UBC Putative Full Length Terpene Synthase Expression Plot",
          subtitle = "Size filtered >= 500 a.a"))

ggsave("results/figures/UBC_TPS_expression.svg", plot = UBC_TPS_exp_plot, 
       width = 4000, height = 2500, units = 'px')

ggsave("results/figures/UBC_TPS_expression.jpg", plot = UBC_TPS_exp_plot, 
       width = 4000, height = 2500, units = 'px')
