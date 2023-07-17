library(dplyr)
library(readr)
library(VennDiagram)


ubc <- read_delim('results/UBC/ubc_dea_results.sigDE.txt', delim = "\t")
jgi <- read_delim('results/UBC/jgi_dea_on_ubc_transcriptome.sigDE.txt', delim = "\t")


### Up-regulated
ubc_up <- ubc |> 
  filter(adj.P.Val <= 0.05 & logFC >=2) |> 
  select(cds)

dim(ubc_up)
# [1] 2003    1

jgi_up <- jgi |> 
  filter(adj.P.Val <= 0.05 & logFC >=2) |> 
  select(cds)

dim(jgi_up)
# [1] 4714    1

(up_intersect <- intersect(ubc_up, jgi_up) |> 
  count())
#       n
# 1  1425

### Down-regulated
### Up-regulated
ubc_down <- ubc |> 
  filter(adj.P.Val <= 0.05 & logFC <=2) |> 
  select(cds)

dim(ubc_down)
# [1] 798   1

jgi_down <- jgi |> 
  filter(adj.P.Val <= 0.05 & logFC <=2) |> 
  select(cds)

dim(jgi_down)
# [1] 3523    1

(down_intersect <- intersect(ubc_down, jgi_down) |> 
    count())
#       n
# 1   508



drawVennDiagram <- function(a1, a2, i) {
  VennDiagram::draw.pairwise.venn(area1 = a1,
                                  area2 = a2,
                                  cross.area = i,
                                  category = c('UBC', 'JGI'),
                                  lwd = c(0, 0),

                                  fill = c('#944C06', '#446C98'),
                                  label.col = c('white', 'white', 'white'),
                                  alpha = c(0.4, 0.4),
                                  cex = c(5, 5, 5),
                                  fontface = 'bold',
                                  fontfamily = 'Helvetica',
                                  
                                  cat.cex = c(5, 5),
                                  cat.col = c('#944C06', '#446C98'),
                                  cat.fontface = 'bold',
                                  cat.fontfamily = 'Helvetica',
  )}

jpeg('results/figures/fig2b_up-regulated.jpg', 
     width = 1600, height = 1600,
     units = 'px', quality = 100)

drawVennDiagram(length(ubc_up$cds), length(jgi_up$cds), up_intersect$n)

dev.off()


jpeg('results/figures/fig2b_down-regulated.jpg', 
     width = 1600, height = 1600,
     units = 'px', quality = 100)

drawVennDiagram(length(ubc_down$cds), length(jgi_down$cds), down_intersect$n)

dev.off()
