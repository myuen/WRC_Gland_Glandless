library(dplyr)
library(ggplot2)
library(stringr)

#####

AKR_exp_plot <- function(stats) {
  
  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 2, aes(shape = factor(type))) +
    scale_x_continuous(limits = c(-bound, bound), 
                       breaks = c(-bound, -10, -5, -2, 0, 2, 5, 10, bound)) + 
    scale_y_discrete(labels = NULL) +
    scale_shape_manual("", values = 1) +

  labs(caption = "p-value <= 0.05") +
  xlab(expression(paste("log"[2],"FC"))) +
  ylab("") +
  
  geom_vline(xintercept = c(-2, 2), colour = "black", linetype = 2) +
  
  theme(
    legend.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, "black"), 
    panel.grid = element_line(colour = "grey95"),
    axis.ticks = element_line(size = 0))
}

#####


UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


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

UBC_AKR_exp_plot <- AKR_exp_plot(UBC_AKR_stats)

(UBC_AKR_exp_plot <- UBC_AKR_exp_plot + 
    ggtitle("UBC Putative Aldo-Keto Reductase (AKR) Expression"))

ggsave("results/figures/UBC_AKR_expression.svg", plot = UBC_AKR_exp_plot)


#####

JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(JGI_dea_stats)
# 'data.frame':	84192 obs. of  4 variables:


JGI_putative_AKRs <- 
  scan("data/targeted-pathway-annotation/05-Aldo-keto-reductase/putative_AKR.JGI_assembly.cdsID.txt",
       what = "character")

JGI_putative_AKRs <- unique(JGI_putative_AKRs)

length(JGI_putative_AKRs)
# [1] 770


JGI_AKR_stats <- JGI_dea_stats %>% 
  filter(cds %in% JGI_putative_AKRs)

str(JGI_AKR_stats)
# 'data.frame':	248 obs. of  4 variables:

JGI_AKR_stats$type <- "Putative"

JGI_AKR_stats$focus <- factor(JGI_AKR_stats$focus, 
                              levels = c("gYoung_glYoung", "gMature_glMature"))

JGI_AKR_exp_plot <- AKR_exp_plot(JGI_AKR_stats)

(JGI_AKR_exp_plot <- JGI_AKR_exp_plot + 
    ggtitle("JGI Putative Aldo-Keto Reductase (AKR) Expression"))

(JGI_AKR_exp_plot <- 
    JGI_AKR_exp_plot + facet_wrap(~ focus, nrow = 2, strip.position = "right"))

ggsave("results/figures/JGI_AKR_expression.svg", plot = JGI_AKR_exp_plot)
