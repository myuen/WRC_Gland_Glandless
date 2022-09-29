library(dplyr)
library(ggplot2)
library(stringr)

#####

AKR_exp_plot <- function(stats) {
  
  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 4, aes(shape = factor(type))) +
    
    scale_x_continuous(limits = c(-bound, bound), 
                       breaks = c(-bound, seq(-10, 10, 5), bound)) + 
    scale_y_discrete(labels = NULL) +
    
    scale_shape_manual(name = "", label = NULL, values = 1) +
    
    labs(caption = "p-value <= 0.05") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +
    
    geom_vline(xintercept = c(-2, 2), colour = "red", linetype = 2) +
    
    annotate(geom = "text", x = 8, y = 20, size = 8,
             label = "4 AKRs with higher\nexpression in gland\n") +
    
    annotate(geom = "text", x = -8, y = 18, size = 8,
             label = "1 AKR with higher\nexpression in glandess") +

    theme(
      title = element_text(size = 20),
      
      legend.position = "none",
      
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, "black"), 
      panel.grid = element_line(colour = "grey99"),
      
      axis.text = element_text(size = 14),
      axis.ticks = element_line(size = 0))
}

#####


UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)


UBC_putative_AKRs <- 
  scan("data/targeted-pathway-annotation/05-Aldo-keto-reductase/putative_AKR.UBC_assembly.cdsID.txt",
       what = "character")

UBC_putative_AKRs <- unique(UBC_putative_AKRs)

length(UBC_putative_AKRs)
# [1] 217


UBC_AKR_stats <- UBC_dea_stats %>% 
  filter(cds %in% UBC_putative_AKRs)

str(UBC_AKR_stats)
# 'data.frame':	100 obs. of  3 variables:

UBC_AKR_stats$type <- "Putative"

UBC_AKR_stats <- UBC_AKR_stats %>% 
  filter(adj.P.Val <= 0.05)

str(UBC_AKR_stats)
# 'data.frame':	23 obs. of  4 variables:

UBC_AKR_exp_plot <- AKR_exp_plot(UBC_AKR_stats)

(UBC_AKR_exp_plot <- UBC_AKR_exp_plot + 
    ggtitle("UBC Putative Aldo-Keto Reductase (AKR) Expression"))

ggsave("results/figures/UBC_AKR_expression.svg", plot = UBC_AKR_exp_plot,
       width = 4000, height = 2500, units = 'px')
