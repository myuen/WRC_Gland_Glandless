library(dplyr)
library(ggplot2)
library(stringr)

#####

SDR_exp_plot <- function(stats) {
  
  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(size = 4, aes(shape = type)) +
    
    scale_x_continuous(limits = c(-bound, bound),
                       breaks = c(-bound, seq(-10, 10, 5), bound)) +
    scale_y_discrete(labels = NULL) +
    scale_shape_manual(name = "",
                       labels = c("Putative", "Cloned"),
                       values = c(1,19)) +
    
    labs(caption = "p-value <= 0.05") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +
    
    geom_vline(xintercept = c(-2, 2), colour = "red", linetype = 2) +
    
    annotate(geom = "text", x = 6, y = 50, size = 6,
             label = "17 SDRs with higher\nexpression in gland\n and 4 cloned") +
    
    annotate(geom = "text", x = -6, y = 50, size = 6,
             label = "9 SDRs with higher\nexpression in glandess") +

    theme(
      text = element_text(size = 16),
      
      title = element_text(size = 20),
      
      legend.background = element_blank(),
      legend.text = element_text(size = 16),
      
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, "black"),
      panel.grid = element_line(colour = "grey99"),
      
      axis.text = element_text(size = 14),
      axis.ticks = element_line(size = 0))
}

#####

UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                        header = TRUE, stringsAsFactors = FALSE)


UBC_putative_SDRs <- 
  scan("data/targeted-pathway-annotation/04-SDR_ADH/putative_SDR.UBC_assembly.size_filtered.cdsID.txt",
       what = "character")
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


UBC_SDR_stats <- inner_join(UBC_dea_stats, UBC_SDRs, by = 'cds')

str(UBC_SDR_stats)
# 'data.frame':	159 obs. of  4 variables:

UBC_SDR_stats$type <- factor(UBC_SDR_stats$type,
                             levels = c("Putative", "Cloned"))

UBC_SDR_stats <- UBC_SDR_stats %>% 
  filter(adj.P.Val <= 0.05)

UBC_SDR_exp_plot <- SDR_exp_plot(UBC_SDR_stats)

(UBC_SDR_exp_plot <- UBC_SDR_exp_plot + 
  ggtitle("UBC Short Chained Alcohol Dehydrogenase (SDR) Expression", 
          subtitle = "Size filtered between 250 to 350 a.a"))

ggsave("results/figures/UBC_SDR_expression.svg", plot = UBC_SDR_exp_plot,
       width = 4000, height = 2500, units = 'px')

