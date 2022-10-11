library(dplyr)
library(ggplot2)
library(stringr)


# Function to create plot
TBB_exp_plot <- function(stats){

  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(y = enzyme, x = logFC)) +   
    geom_point(size = 3, colour = "black") + 
    
    facet_wrap(~ pathway, drop = TRUE, 
               scale = "free_y", ncol = 1,
               strip.position = "right") +
  
    scale_x_continuous(limits = c(-bound, bound)) + 
    
    geom_vline(xintercept = c(-2, 2), colour = "red", linetype = 2) + 
    
    labs(caption = "p-value <= 0.05") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +
    
    theme(
      title = element_text(size = 20),
      
      strip.text.y = element_text(face = "bold", size = 16),
      
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),

      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, "black"),
      panel.grid = element_line(colour = "grey95"))
}


# Read differential expression results
UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)


# Read putative TBB CDS IDs from functional domain scan with hmmscan and rpsblast.
UBC_TBB_orthologs <- 
  read.table("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/UBC_TBB_orthologs.txt",
             stringsAsFactors = FALSE, header = TRUE)

str(UBC_TBB_orthologs)
# 'data.frame':	60 obs. of  4 variables:


# Inner join to keep only enzyme involved in terpenoid backbone biosynthesis
UBC_TBB_ortholog_stats <- 
  inner_join(UBC_dea_stats, UBC_TBB_orthologs, by = 'cds')


# Filter gene by statistics significance
UBC_TBB_ortholog_stats <- UBC_TBB_ortholog_stats %>%
  filter(adj.P.Val <= 0.05)

str(UBC_TBB_ortholog_stats)
# 'data.frame':	35 obs. of  6 variables:

# Relevel pathway factor
UBC_TBB_ortholog_stats$pathway <- 
  factor(UBC_TBB_ortholog_stats$pathway, 
         levels = c("MEP", "MEV", "C1020"))

UBC_TBB_plot <-
  TBB_exp_plot(UBC_TBB_ortholog_stats)

(UBC_TBB_plot <-
    UBC_TBB_plot + 
      ggtitle("Terpenoid Backbone Biosynthesis"))

ggsave("results/figures/UBC_TBB_expression.svg", plot = UBC_TBB_plot,
       width = 4000, height = 2500, units = 'px')

ggsave("results/figures/UBC_TBB_expression.jpg", plot = UBC_TBB_plot,
       width = 4000, height = 2500, units = 'px')
