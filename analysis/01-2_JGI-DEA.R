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
      pattern = "JGI\\S+quant.sf",
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
  str_extract("gland\\S*_quant.sf") %>% 
  str_replace("_foliage", "") |> 
  str_replace("_quant.sf", "")

names(quant.files) <- samples


metadata <-
  str_split(samples, "_")


metadata <-
  map_df(metadata, function(l){
    df = data.frame("phenotype" = l[1],
                    "age" = l[2],
                    "replicate" = l[3]
    )
    df
  })


metadata$phenotype <- 
  factor(metadata$phenotype, levels = c("glandless", "gland"))

metadata$age <- 
  factor(metadata$age, levels = c("young", "mature"))


# Create experimental design
expDes <- 
  data.frame(
    sample = samples,
    phenotype = metadata$phenotype,
    age = metadata$age,
    replicate = metadata$replicate)

#                   sample phenotype    age replicate
# 1      gland_mature_rep1     gland mature      rep1
# 2      gland_mature_rep2     gland mature      rep2
# 3      gland_mature_rep3     gland mature      rep3
# 4      gland_mature_rep4     gland mature      rep4
# 5      gland_mature_rep5     gland mature      rep5
# 6       gland_young_rep1     gland  young      rep1
# 7       gland_young_rep2     gland  young      rep2
# 8       gland_young_rep3     gland  young      rep3
# 9       gland_young_rep4     gland  young      rep4
# 10      gland_young_rep5     gland  young      rep5
# 11      gland_young_rep6     gland  young      rep6
# 12 glandless_mature_rep1 glandless mature      rep1
# 13 glandless_mature_rep2 glandless mature      rep2
# 14 glandless_mature_rep3 glandless mature      rep3
# 15 glandless_mature_rep4 glandless mature      rep4
# 16 glandless_mature_rep5 glandless mature      rep5
# 17  glandless_young_rep1 glandless  young      rep1
# 18  glandless_young_rep2 glandless  young      rep2
# 19  glandless_young_rep3 glandless  young      rep3
# 20  glandless_young_rep4 glandless  young      rep4
# 21  glandless_young_rep5 glandless  young      rep5



x <- DGEList(counts = txi$counts)

dim(x)
# [1] 248844     21


# Filter low-expression contigs.  Keep only genes with at 
# least 1 count-per-million reads (cpm) in at least 3 libraries (i.e. half of 
# the number of libraries per condition)
y <- x[(rowSums(cpm(x) > 1) >= 3), ]

dim(y)
# [1] 42096    21


# Reset depth
y$samples$lib.size <- colSums(y$counts)


# TMM Normalization by Depth
y <- calcNormFactors(y)

CPM <- cpm(y)

colnames(CPM) <- expDes$sample

write.table(CPM, "results/jgi_normalized_cpm.txt",
            sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)


# Create model matrix
modMat <-
  model.matrix(~ phenotype * age, expDes)

colnames(modMat)[1] <- "glandless"

colnames(modMat) <-
  colnames(modMat) %>% 
  str_replace("phenotype", "") |> 
  str_replace("age", "") |> 
  str_replace(":", "_")

#    glandless gland mature gland_mature
# 1          1     1      1            1
# 2          1     1      1            1
# 3          1     1      1            1
# 4          1     1      1            1
# 5          1     1      1            1
# 6          1     1      0            0
# 7          1     1      0            0
# 8          1     1      0            0
# 9          1     1      0            0
# 10         1     1      0            0
# 11         1     1      0            0
# 12         1     0      1            0
# 13         1     0      1            0
# 14         1     0      1            0
# 15         1     0      1            0
# 16         1     0      1            0
# 17         1     0      0            0
# 18         1     0      0            0
# 19         1     0      0            0
# 20         1     0      0            0
# 21         1     0      0            0



# voom transformation
v <- voom(y, modMat, plot = TRUE)


# Create PCA plot
(p <- PCA_maker(expDes, v))

ggsave("results/figures/jgi_pca.svg", plot = p)


cont_matrix <- makeContrasts(
  # young gland vs young glandless, complementing 
  # the comparison in UBC experiment
  gYoung_glYoung = gland - glandless,
  # glandless mature vs glandless young
  glMature_glYoung = mature,
  # gland mature vs gland young
  gMature_gYoung = gland_mature,
  levels = modMat)


# Linear modeling
fit <- lmFit(v, modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)


# Up-regulation = higher expression in gland in reference to glandless
summary(decideTests(fit3, method = "separate", 
                    adjust.method = "fdr", p.value = pCutoff, 
                    lfc = lfcCutoff))

#        gYoung_glYoung glMature_glYoung gMature_gYoung
# Down            19034                2              0
# NotSig          20361            42094          42096
# Up               2701                0              0


focus <- colnames(cont_matrix)

all_results <-
  map_df(focus, function(f) {
    ret <- topTable(fit3, coef = f, sort.by = "logFC", number = Inf)
    ret <- ret |> add_column("focus" = f)
    ret <- ret |> rownames_to_column("cds")
    ret <- ret |> select(cds, focus, logFC, adj.P.Val)
})

str(all_results)
# 'data.frame':	126288 obs. of  4 variables:

table(all_results$focus)
# glMature_glYoung   gMature_gYoung   gYoung_glYoung 
#            42096            42096            42096 


### Write out stats of all contigs
write.table(
  all_results, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/jgi_dea_results.txt"
)


### Write out all sig DE stats
sig_results <- all_results %>% 
  filter(adj.P.Val <= pCutoff & abs(logFC) >= lfcCutoff)

str(sig_results)
# 'data.frame':	21737 obs. of  4 variables:


write.table(
  sig_results, quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE,
  "results/jgi_dea_sig_results.txt"
)

write(sig_results$cds, "results/jgi_dea_sig.cdsID.txt")


sig_upReg <- sig_results %>% 
  filter(logFC >= lfcCutoff)

str(sig_upReg)
# 'data.frame':	2701 obs. of  4 variables:

write(sig_upReg$cds, "results/jgi_dea_upReg.cdsID.txt")


sig_downReg <- sig_results %>% 
  filter(logFC <= lfcCutoff)

str(sig_downReg)
# 'data.frame':	19036 obs. of  4 variables:

write(unique(sig_downReg$cds), "results/jgi_dea_downReg.cdsID.txt")
