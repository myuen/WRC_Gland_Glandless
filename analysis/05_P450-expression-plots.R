library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)


#####

# geom_points shapes in preferencial order
# shapes <- c(0, 21, 23, 17)

P450_exp_plot <- function(stats, pt_shapes) {
  
  bound <- ceiling(max(abs(stats$logFC)))
  
  ggplot(stats, aes(x = logFC, y = cds)) +
    geom_point(aes(size = type, shape = type, colour = type, fill = type)) +
    
    scale_x_continuous(limits = c(-bound, bound), 
                       breaks = c(-bound, -10, -5, -2, 0, 2, 5, 10, bound)) + 
    scale_y_discrete(labels = NULL) +

    scale_alpha_manual(name = "",
                       values = c("Putative" = 0.1, 
                                  "Annotated" = 1, 
                                  "CYP750B1" = 1, 
                                  "CYP76AA25" = 1)) +
    scale_colour_manual(name = "",
                        values = c("black", "deepskyblue4", 
                                   "red", "brown4")) +

    scale_fill_manual(name = "",
                      values = c("black", "deepskyblue4", 
                                 "red", "brown4")) +
    
    scale_shape_manual(name = "",
                       values = c(1, 21, 24, 22)) +

    scale_size_manual(name = "",
                      values = c("Putative" = 3, 
                                 "Annotated" = 5, 
                                 "CYP750B1" = 5, 
                                 "CYP76AA25" = 5)) +
    
    
    

    labs(caption = "p-value <= 0.05") +
    xlab(expression(paste("log"[2],"FC"))) +
    ylab("") +
    
    geom_vline(xintercept = c(-2, 2), colour = "red", linetype = 2) +
    
    theme(
      title = element_text(size = 20),
      
      legend.background = element_blank(),
      legend.text = element_text(size = 16),
          
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, "black"), 
      panel.grid = element_line(colour = "white"),

      axis.text = element_text(size = 14),
      axis.ticks = element_line(size = 0))
}


######


UBC_dea_stats <- read.delim("results/ubc_dea_results.txt",
                            header = TRUE, stringsAsFactors = FALSE)
str(UBC_dea_stats)
# 'data.frame':	33989 obs. of  3 variables:


# Read putative P450s CDS IDs from functional domain scan with hmmscan and rpsblast.
UBC_putative_P450s <- 
  scan("data/targeted-pathway-annotation/03-P450/putative_P450.UBC_assembly.cdsID.txt",
       what = "character")
# Read 1066 items


# Construct dataframe.  Easier to remove duplicates with CDS ID as rownames
UBC_P450s <- data.frame(row.names = c(unique(UBC_putative_P450s)))

# Default annotation as "Putative"
UBC_P450s$type <- "Putative"

str(UBC_P450s)
# 'data.frame':	1066 obs. of  1 variable:


# Read annotated TPS in assembly
UBC_annotated_P450s <- 
  read.delim("data/targeted-pathway-annotation/03-P450/P450s.blastp.UBC_assembly.txt",
             stringsAsFactors = FALSE, header = FALSE)


# Get CDS ID for cloned P450s
# AKH41019.1 = CYP750B1
# AKH41025.1 = CYP76AA25
CYP750B1 <- 
  UBC_annotated_P450s[UBC_annotated_P450s$V1 == "AKH41019.1",]$V2 

CYP76AA25 <- 
  UBC_annotated_P450s[UBC_annotated_P450s$V1 == "AKH41025.1",]$V2


# Filter out low percent positive (V7) from BLASTp
UBC_annotated_P450s <- UBC_annotated_P450s %>%  
  filter(V7 >= 90) %>%  
  select(V2) %>% unique()

str(UBC_annotated_P450s)
# 'data.frame':	8 obs. of  1 variable:

UBC_P450s[UBC_annotated_P450s$V2, "type"] <- "Annotated"

UBC_P450s[CYP750B1,] <- "CYP750B1"
UBC_P450s[CYP76AA25,] <- "CYP76AA25"

# Turn rowname to column
UBC_P450s <- UBC_P450s %>% rownames_to_column("cds")

str(UBC_P450s)
# 'data.frame':	1066 obs. of  2 variables:


UBC_P450_stats <- inner_join(UBC_dea_stats, UBC_P450s)
# Joining, by = "cds"

str(UBC_P450_stats)  
# 'data.frame':	472 obs. of  4 variables:

UBC_P450_stats$type <-
  factor(UBC_P450_stats$type,
         levels = c("Putative", "Annotated", "CYP750B1", "CYP76AA25"))

# Number of types
num_types <- length(unique(UBC_P450s$type))

UBC_P450_exp_plot <- 
  P450_exp_plot(UBC_P450_stats, shapes[1:num_types])

(UBC_P450_exp_plot <- UBC_P450_exp_plot + 
    ggtitle("UBC P450s Expression Plot"))

ggsave("results/figures/UBC_P450_expression.svg", plot = UBC_P450_exp_plot,
       width = 4000, height = 2500, units = 'px')


######

# 
# JGI_dea_stats <- read.delim("results/jgi_dea_results.txt",
#                             header = TRUE, stringsAsFactors = FALSE)
# str(JGI_dea_stats)
# # 'data.frame':	84192 obs. of  4 variables:
# 
# 
# # Read putative P450s CDS IDs from functional domain scan with hmmscan and rpsblast.
# JGI_putative_P450s <- 
#   scan("data/targeted-pathway-annotation/03-P450/putative_P450.JGI_assembly.cdsID.txt",
#        what = "character")
# # Read 1763 items
# 
# 
# # Construct dataframe.  Easier to remove duplicates with CDS ID as rownames
# JGI_P450s <- data.frame(row.names = c(unique(JGI_putative_P450s)))
# 
# # Default annotation as "Putative"
# JGI_P450s$type <- "Putative"
# 
# str(JGI_P450s)
# # 'data.frame':	1763 obs. of  1 variable:
# 
# 
# # Read annotated P450s in assembly
# JGI_annotated_P450s <- 
#   read.delim("data/targeted-pathway-annotation/03-P450/P450s.blastp.JGI_assembly.txt",
#              stringsAsFactors = FALSE, header = FALSE)
# 
# 
# # Get CDS ID for cloned P450s
# # AKH41019.1 = CYP750B1
# # AKH41025.1 = CYP76AA25
# CYP750B1 <- 
#   JGI_annotated_P450s[JGI_annotated_P450s$V1 == "AKH41019.1",]$V2 
# 
# CYP76AA25 <- 
#   JGI_annotated_P450s[JGI_annotated_P450s$V1 == "AKH41025.1",]$V2
# 
# 
# # Filter out low percent positive (V7) from BLASTp
# JGI_annotated_P450s <- JGI_annotated_P450s %>% 
#   filter(V7 >= 90) %>% 
#   select(V2) %>% unique()
# 
# str(JGI_annotated_P450s)
# # 'data.frame':	9 obs. of  1 variable:
# 
# JGI_P450s[JGI_annotated_P450s$V2, "type"] <- "Annotated"
# 
# JGI_P450s[CYP750B1,] <- "CYP750B1"
# JGI_P450s[CYP76AA25,] <- "CYP76AA25"
# 
# # Turn rowname to column
# JGI_P450s <- JGI_P450s %>% rownames_to_column("cds")
# 
# str(JGI_P450s)
# # 'data.frame':	1763 obs. of  2 variables:
# 
# 
# JGI_P450_stats <- inner_join(JGI_dea_stats, JGI_P450s)
# # Joining, by = "cds"
# 
# str(JGI_P450_stats)  
# # 'data.frame':	1452 obs. of  5 variables:
# 
# JGI_P450_stats$type <-
#   factor(JGI_P450_stats$type,
#          levels = c("Putative", "Annotated", "CYP750B1", "CYP76AA25"))
# 
# JGI_P450_stats$focus <- 
#   factor(JGI_P450_stats$focus,
#          levels = c("gYoung_glYoung", "gMature_glMature"))
# 
# # Number of types
# num_types <- length(unique(JGI_P450s$type))
# 
# JGI_P450_exp_plot <- 
#   P450_exp_plot(JGI_P450_stats, shapes[1:num_types])
# 
# (JGI_P450_exp_plot <- JGI_P450_exp_plot + 
#     ggtitle("JGI P450s Expression Plot"))
# 
# (JGI_P450_exp_plot <- JGI_P450_exp_plot + 
#   facet_wrap(~ focus, nrow = 2, strip.position = "right"))
# 
# ggsave("results/figures/JGI_P450_expression.svg", plot = JGI_P450_exp_plot)
