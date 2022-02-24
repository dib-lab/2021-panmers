library(readr)
library(dplyr)
library(pagoo)

# read in roary csv as pagoo 
roary_pg <- roary_2_pagoo(snakemake@input[['pa']])
# reset core level to genes present in 100% of samples
roary_pg$core_level <- 100
# extract core genes and get cluster names
core_df <- unlist(roary_pg$core_genes)
core_cluster_names <- core_df %>%
  as.data.frame() %>%
  select(cluster) %>%
  distinct()
# read in sequence names
seq_names <- read_delim(snakemake@input[['names']], delim = " ", 
                        col_names = c("gene", "cluster"), col_types = "cc") %>%
  mutate(gene = gsub(">", "", gene))
# find core gene sequence names
core_seq_names <- seq_names %>%
  filter(cluster %in% core_cluster_names$cluster)
# write core seq names for intake by seqtk
write.table(core_seq_names$gene, snakemake@output[["core_names"]], row.names = F, quote = F, col.names = F)
