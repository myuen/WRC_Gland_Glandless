library(dplyr)
library(purrr)
library(stringr)

#####
#Terpenoid backbone biosynthesis (TBB) from UBC data

UBC_blast_result_files <- 
  list.files("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/",
             pattern = "*UBC_assembly.txt", full.names = TRUE, recursive = TRUE)

length(UBC_blast_result_files)
# [1] 18


#BLAST was run in batch and some result files can be empty
#Remove empty result files
UBC_blast_result_files <- 
  UBC_blast_result_files[file.info(UBC_blast_result_files)$size != 0]

length(UBC_blast_result_files)
# [1] 18


ubc_tbb_blast_results <-
  map_df(UBC_blast_result_files, function(f){
  df <- read.delim(f, sep = "\t", 
                    stringsAsFactors = FALSE, header = FALSE)

  metadata <- f |> 
    str_match("(\\w+)-(\\d)-(\\w+).fasta.blastp_UBC_assembly.txt") |> 
    as.data.frame()

  df$pathway = metadata[,2]
  df$step = metadata[,3]
  df$enzyme = metadata[,4]
  return(df)
})

colnames(ubc_tbb_blast_results)[1:2] <- c("query", "cds")

ubc_tbb_blast_results <- ubc_tbb_blast_results |> 
  select("cds", "pathway", "step", "enzyme") |> 
  unique()

str(ubc_tbb_blast_results)
# 'data.frame':	60 obs. of  4 variables:

write.table(ubc_tbb_blast_results, 
            "data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/UBC_tbb_orthologs.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#####

#Terpenoid backbone biosynthesis (TBB) from JGI data

JGI_blast_result_files <- 
  list.files("data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/",
             pattern = "*JGI_assembly.txt", full.names = TRUE, recursive = TRUE)

length(JGI_blast_result_files)
# [1] 18


jgi_tbb_blast_results <-
  map_df(JGI_blast_result_files, function(f){
    df <- read.delim(f, sep = "\t", 
                     stringsAsFactors = FALSE, header = FALSE)
    
    metadata <- f |> 
      str_match("(\\w+)-(\\d)-(\\w+).fasta.blastp_JGI_assembly.txt") |> 
      as.data.frame()
    
    df$pathway = metadata[,2]
    df$step = metadata[,3]
    df$enzyme = metadata[,4]
    return(df)
  })

colnames(jgi_tbb_blast_results)[1:2] <- c("query", "cds")

jgi_tbb_blast_results <- jgi_tbb_blast_results |> 
  select("cds", "pathway", "step", "enzyme") |> 
  unique()

str(jgi_tbb_blast_results)
# 'data.frame':	90 obs. of  4 variables:

write.table(jgi_tbb_blast_results, 
            "data/targeted-pathway-annotation/01-Terpenoid_backbone_biosynthesis/JGI_tbb_orthologs.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
