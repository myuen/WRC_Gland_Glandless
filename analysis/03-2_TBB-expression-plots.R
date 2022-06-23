library(dplyr)
library(ggplot2)
library(stringr)

#####

TBB_exp_plot <- function(stats, title){
  
  ggplot(stats, aes(y = enzyme, x = logFC)) +   
    geom_point(size = 2, colour = "black") + 
    
    facet_wrap(~ pathway, drop = TRUE, 
               scale = "free_y", ncol = 1,
               strip.position = "right") +
    
    geom_vline(xintercept = c(-2, 2), colour = "black", linetype = 2) + 
    
    ggtitle(title) +
    labs(caption = "p-value <= 0.05") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +
    
    geom_vline(xintercept = c(-2, 2), colour = "black", linetype = 2) +
    
    theme(strip.text.y = element_text(face = "bold", size = 10),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, "black"),
          panel.grid = element_line(colour = "grey90"))
}

#####

UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


# Read putative TBB CDS IDs from functional domain scan with hmmscan and rpsblast.
UBC_tbb_orthologs <- 
  read.table("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/UBC_tbb_orthologs.txt",
             stringsAsFactors = FALSE, header = TRUE)
str(UBC_tbb_orthologs)
# 'data.frame':	60 obs. of  4 variables:


UBC_tbb_ortholog_stats <- 
  inner_join(UBC_dea_stats, UBC_tbb_orthologs)
# Joining, by = "cds"

str(UBC_tbb_ortholog_stats)
# 'data.frame':	54 obs. of  6 variables:

UBC_tbb_ortholog_stats$pathway <- 
  factor(UBC_tbb_ortholog_stats$pathway, 
         levels = c("MEP", "MEV", "C1020"))

(UBC_tbb_plot <-
  TBB_exp_plot(UBC_tbb_ortholog_stats, "UBC Terpenoid Backbone Biosynthesis"))

ggsave("results/figures/UBC_TBB_expression.svg", plot = UBC_tbb_plot)

#####

JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(JGI_dea_stats)
# 'data.frame':	84192 obs. of  4 variables:


JGI_tbb_orthologs <- 
  read.table("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/JGI_tbb_orthologs.txt",
             stringsAsFactors = FALSE, header = TRUE)
str(JGI_tbb_orthologs)
# 'data.frame':	90 obs. of  4 variables:


JGI_tbb_ortholog_stats <- 
  inner_join(JGI_dea_stats, JGI_tbb_orthologs)
# Joining, by = "cds"

str(JGI_tbb_ortholog_stats)
# 'data.frame':	128 obs. of  7 variables:


JGI_tbb_ortholog_stats$pathway <- 
  factor(JGI_tbb_ortholog_stats$pathway, 
         levels = c("MEP", "MEV", "C1020"))


(JGI_tbb_plot <-
    TBB_exp_plot(JGI_tbb_ortholog_stats, "JGI Terpenoid Backbone Biosynthesis"))

ggsave("results/figures/JGI_TBB_expression.svg", plot = JGI_tbb_plot)
