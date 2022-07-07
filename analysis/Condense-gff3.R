library(dplyr)
library(purrr)
library(stringr)


### Read GFF3 file for genome annotation
gff3 <- read.delim("data/Tplicatav3.1c.primaryTrs.gff3", header = FALSE,
                   comment.char = "#", sep = "\t", stringsAsFactors = FALSE)

#Only keeping genes from GFF3 file
gff3_mRNA <- gff3 %>%
  filter(V3 == "mRNA") %>%
  select(1,4,5,9)

colnames(gff3_mRNA) <- 
  c("scaffold", "start", "end", "attributes")

str(gff3_mRNA)
# 'data.frame':	39659 obs. of  4 variables:

gff3_mRNA$mRNA <- map(gff3_mRNA$attributes, function(a){
  m <- str_split(a, ";") %>% unlist()
  str_replace(m[2], "Name=", "")
}) %>% unlist()

gff3_mRNA <- gff3_mRNA %>% 
  select(scaffold, mRNA, start, end)

write.table(gff3_mRNA, "data/gff3_condensed.txt", sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)
