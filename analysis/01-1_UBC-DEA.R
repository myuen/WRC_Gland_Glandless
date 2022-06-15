library(edgeR)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tximport)


source("analysis/helper01_PCA-maker.R")


### Differential Expression Analysis on gland and glandless genotype 
### in Western Redcedar

# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


# Import file with tximport
quant.files <-
  dir("data",
      pattern = "UBC\\S+quant.sf",
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE)


# Read Salmon quantification with tximport
txi <-
  tximport(quant.files,
           type = "salmon",
           txOut = TRUE,
           varReduce = TRUE,
           countsFromAbundance = "lengthScaledTPM")


samples <-
  quant.files %>% 
  str_extract("\\w{3}_gland\\S*_quant.sf") %>% 
  str_replace("_quant.sf", "")

names(quant.files) <- samples

metadata <-
  str_split(samples, "_")

metadata <-
  map_df(metadata, function(l){
    df = data.frame("phenotype" = l[2],
                    "replicate" = l[3]
    )
    df
})


metadata$phenotype <- 
  factor(metadata$phenotype, levels = c("glandless", "gland"))


# Create experimental design
expDes <- 
  data.frame(
    sample = samples,
    phenotype = metadata$phenotype,
    replicate = metadata$replicate)

#            sample phenotype replicate
# 1     UBC_gland_A     gland         A
# 2     UBC_gland_B     gland         B
# 3     UBC_gland_C     gland         C
# 4     UBC_gland_D     gland         D
# 5 UBC_glandless_A glandless         A
# 6 UBC_glandless_B glandless         B
# 7 UBC_glandless_C glandless         C
# 8 UBC_glandless_D glandless         D


x <- DGEList(counts = txi$counts, group = expDes$phenotype)

dim(x)
# [1] 97018     8


# Filter low-expression contigs.  Keep only genes with at 
# least 1 count-per-million reads (cpm) in at least 2 libraries (i.e. half of 
# the number of libraries per condition)
y <- x[(rowSums(cpm(x) > 1) >= 2), ]

dim(y)
# [1] 33989     8


# Reset depth
y$samples$lib.size <- colSums(y$counts)


# TMM Normalization by Depth
y <- calcNormFactors(y)

CPM <- cpm(y)

colnames(CPM) <- expDes$sample

write.table(CPM, "results/ubc_normalized_cpm.txt", 
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)


# make model matrix
modMat <-
  model.matrix(~ 0 + phenotype, expDes)

colnames(modMat) <- 
  colnames(modMat) %>% str_replace("phenotype", "")

#   glandless gland
# 1         0     1
# 2         0     1
# 3         0     1
# 4         0     1
# 5         1     0
# 6         1     0
# 7         1     0
# 8         1     0


# voom transformation
v <- voom(y, modMat, plot = TRUE)


# Create PCA plot
(p <- PCA_maker(expDes, v))

ggsave("results/figures/ubc_pca.svg", plot = p)


cont_matrix <- makeContrasts(
  g_gl = gland - glandless,
  levels = modMat)


# Linear modeling
fit <- lmFit(v, modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

# Up-regulation = higher expression in gland in reference to glandes
summary(decideTests(fit3, method = "separate", 
                    adjust.method = "fdr", p.value = pCutoff, 
                    lfc = lfcCutoff))

#         g_gl
# Down     673
# NotSig 31590
# Up      1726


results <- topTable(fit3, sort.by = "logFC", number = Inf)

results <- results %>% rownames_to_column("cds")

str(results)
# 'data.frame':	33989 obs. of  7 variables:


# Reorganize columns
results <- results %>% 
  select(cds, logFC, adj.P.Val)


### Write out stats of all contigs
write.table(
  results, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/ubc_dea_results.txt"
)


### Write out all sig DE stats
sig_results <- results %>% 
  filter(adj.P.Val <= pCutoff & abs(logFC) >= lfcCutoff)

str(sig_results)
# 'data.frame':	2399 obs. of  3 variables:

write.table(
  sig_results, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/ubc_dea_sig_results.txt"
)

write(sig_results$cds, "results/ubc_dea_sig.cdsID.txt")
