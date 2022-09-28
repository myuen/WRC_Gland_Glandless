library(dplyr)
library(ggplot2)


# Read differential expressed statistics results
UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


x_bound <- max(ceiling(abs(UBC_dea_stats$logFC)))
y_bound <- max(ceiling(-log10(UBC_dea_stats$adj.P.Val)))


(volcano <- ggplot(UBC_dea_stats, aes(x = logFC, y = -log10(adj.P.Val))) + 
    geom_point(shape = 21, colour = 'black', fill = 'black', alpha = 0.3) +
    
    # Add break points on x-axis to show min, max, log2FC threshold
    scale_x_continuous(breaks = c(seq(-x_bound, x_bound, 5))) +
    
    scale_y_continuous(breaks = c(seq(0, y_bound, 2))) +
    
    theme_classic() +
    
    # Define text size and remove grid lines
    theme(text = element_text(size = 16),
          legend.text = element_text(size = 16),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    
    # Add line to show log2FC and adj p-value threshold
    geom_hline(yintercept = -log10(0.05), 
               colour = "blue", linetype = 2, size = 0.8) +
    geom_vline(xintercept = c(-2, 2), colour = "red", 
               linetype = 2, size = 0.8) + 
    
    # Add title for plot, x and y-axis
    labs(title = "Differential Expression Results",
         x = expression(paste('log'[2], ' Fold Change')), 
         y = expression(paste('-log'[10], ' Adj. p-Value'))) + 
    
    annotate("text", x = 8, y = 6.5, size = 6,
             label = expression(bold("33,990 contigs in total"))) +
    
    annotate("text", x = 10, y = 4, size = 5,
             label = expression(bold("1,726 with higher\nexpression in Gland"))) +
    
    annotate("rect", xmin = 2, xmax = x_bound, 
             ymin = 1.3, ymax = y_bound, alpha = .1) +
    
    annotate("text", x = -10, y = 4, size = 5,
             label = expression(bold("673 with higher\nexpressionin Glandless"))) +
    
    annotate("rect", xmin = -2, xmax = -x_bound,
             ymin = 1.3, ymax = y_bound, alpha = .1)
)

ggsave("results/figures/UBC_volcano.svg", plot = volcano,
       width = 4000, height = 2500, units = 'px')
