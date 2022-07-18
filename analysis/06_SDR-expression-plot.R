library(dplyr)
library(ggplot2)
library(stringr)

#####

SDR_exp_plot <- function(stats) {

  bound <- ceiling(max(abs(stats$logFC)))

  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 2.5, aes(shape = type)) +
  scale_x_continuous(limits = c(-bound, bound),
                     breaks = c(-bound, 10, -5, -2, 0, 2, 5, 10, bound)) +
  scale_y_discrete(labels = NULL) +
    scale_shape_manual(name = "",
                       labels = c("Putative", "Cloned"),
                       values = c(1,19)) +

  labs(caption = "p-value <= 0.05") +
  xlab(expression(paste("log"[2],"FC"))) +
  ylab("") +

  geom_vline(xintercept = c(-2, 2), colour = "black", linetype = 2) +

  theme(
    legend.background = element_blank(),
    # strip.text.y = element_text(face = "bold", size = 12),
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


UBC_putative_SDRs <- 
  scan("data/targeted-pathway-annotation/04-SDR_ADH/putative_SDR.UBC_assembly.size_filtered.cdsID.txt",
       what = "character")

UBC_putative_SDRs <- unique(UBC_putative_SDRs)

length(UBC_putative_SDRs)
# Read 357 items


UBC_cloned_SDRs <- 
  scan("data/targeted-pathway-annotation/04-SDR_ADH/cloned_SDR_UBC_cdsID.txt", 
       what = "character")

UBC_cloned_SDRs <- unique(UBC_cloned_SDRs)

length(UBC_cloned_SDRs)
# [1] 5


UBC_SDRs <- 
  data.frame(rownames = c(union(UBC_putative_SDRs, UBC_cloned_SDRs)))
                       
UBC_SDRs$type <- "Putative"

colnames(UBC_SDRs)[1] <- 'cds'

str(UBC_SDRs)
# 'data.frame':	359 obs. of  2 variables:


UBC_SDRs[UBC_SDRs$cds %in% UBC_cloned_SDRs,][2] <- "Cloned"

str(UBC_SDRs)
# 'data.frame':	357 obs. of  2 variables:


UBC_SDR_stats <- inner_join(UBC_dea_stats, UBC_SDRs)
# Joining, by = "cds"

str(UBC_SDR_stats)
# 'data.frame':	159 obs. of  4 variables:

UBC_SDR_stats$type <- factor(UBC_SDR_stats$type,
                             levels = c("Putative", "Cloned"))

UBC_SDR_exp_plot <- SDR_exp_plot(UBC_SDR_stats)

(UBC_SDR_exp_plot <- UBC_SDR_exp_plot + 
  ggtitle("UBC Short Chained Alcohol Dehydrogenase (SDR) Expression", 
          subtitle = "Size filtered between 250 to 350 a.a"))

ggsave("results/figures/UBC_SDR_expression.svg", plot = UBC_SDR_exp_plot)


#####
# 
# JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
#                             header = TRUE, stringsAsFactors = FALSE)
# str(JGI_dea_stats)
# # 'data.frame':	84192 obs. of  4 variables:
# 
# 
# JGI_putative_SDRs <- 
#   scan("data/targeted-pathway-annotation/04-SDR_ADH/putative_SDR.JGI_assembly.size_filtered.cdsID.txt",
#        what = "character")
# 
# JGI_putative_SDRs <- unique(JGI_putative_SDRs)
# 
# length(JGI_putative_SDRs)
# # [1] 383
# 
# 
# JGI_cloned_SDRs <- 
#   scan("data/targeted-pathway-annotation/04-SDR_ADH/cloned_SDR_JGI_cdsID.txt", 
#        what = "character")
# 
# JGI_cloned_SDRs <- unique(JGI_cloned_SDRs)
# 
# length(JGI_cloned_SDRs)
# # [1] 3
# 
# 
# JGI_SDRs <- data.frame(rownames = c(unique(JGI_putative_SDRs, JGI_cloned_SDRs)))
# 
# JGI_SDRs$type <- "Putative"
# 
# colnames(JGI_SDRs)[1] <- 'cds'
# 
# str(JGI_SDRs)
# # 'data.frame':	383 obs. of  2 variables:
# 
# 
# JGI_SDRs[JGI_SDRs$cds %in% JGI_cloned_SDRs,][2] <- "Cloned"
# 
# str(JGI_SDRs)
# # 'data.frame':	699 obs. of  2 variables:
# 
# 
# JGI_SDR_stats <- inner_join(JGI_dea_stats, JGI_SDRs)
# # Joining, by = "cds"
# 
# str(JGI_SDR_stats)
# # 'data.frame':	114 obs. of  5 variables:
# 
# JGI_SDR_stats$type <- factor(JGI_SDR_stats$type,
#                              levels = c("Putative", "Cloned"))
# 
# JGI_SDR_stats$focus <- factor(JGI_SDR_stats$focus,
#                              levels = c("gYoung_glYoung", "gMature_glMature"))
# 
# JGI_SDR_exp_plot <- SDR_exp_plot(JGI_SDR_stats)
# 
# (JGI_SDR_exp_plot <- JGI_SDR_exp_plot +
#     ggtitle("JGI Short Chained Alcohol Dehydrogenase (SDR) Expression", 
#           subtitle = "Size filtered between 250 to 350 a.a") +
#     facet_wrap(~ focus, nrow = 2, strip.position = "right"))
# 
# ggsave("results/figures/JGI_SDR_expression.svg", plot = JGI_SDR_exp_plot)
