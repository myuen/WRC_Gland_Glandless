# Function to create PCA plot

# Argument 1 = experimental design
# Argument 2 = EList object from limma

PCA_maker <- function(expDes, v) {
  require(ggplot2)
  
  pca <- prcomp(t(v$E), scale. = TRUE)
  data <- cbind(pca$x[,1:2], expDes)

  # Add 25% margin on both axis to the max value
  limit <- max(abs(c(data$PC1, data$PC2))) * 1.25
  
  pca_summary <- summary(pca)
  
  pc1 <- pca_summary$importance["Proportion of Variance", "PC1"] * 100
  pc2 <- pca_summary$importance["Proportion of Variance", "PC2"] * 100
  

  p <-
    ggplot(
    data,
    aes_string(x = "PC1", y = "PC2", shape = data$phenotype, color = data$age)) +
    geom_point(size = 2) +
    
    # Naming title and axis
    ggtitle("Principal Component Analysis") + 
    xlab(paste("PC1 (", pca_summary$importance["Proportion of Variance", "PC1"]*100, "%)")) +
    ylab(paste("PC2 (", pca_summary$importance["Proportion of Variance", "PC2"]*100, "%)")) +
    
    # Text on labels
    geom_text(aes(label = sample),
              size = 2, vjust = 2., hjust = 0.25) +

    scale_x_continuous(limits = c(-limit, limit)) + 
    scale_y_continuous(limits = c(-limit, limit)) + 
    
    theme_bw(base_size = 6) + 
    
    
    # Setting legend colour and shape
    guides(shape = guide_legend("Phenotype"),
           color = guide_legend("Age")) +
    
    scale_colour_manual(values = c("firebrick1", "dodgerblue3"))
  
  return(p)
}
