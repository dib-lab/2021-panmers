library(dplyr)
library(readr)
library(janitor)
library(tibble)
library(arrow)
library(ggplot2)
library(vegan)
setwd("~/github/2021-panmers/")

roary <- read_csv("outputs/roary/s__Faecalibacterium-prausnitzii_D/gene_presence_absence.csv") %>%
  clean_names() %>%
  select(-non_unique_gene_name, -annotation, -no_isolates, -no_sequences, 
         -avg_sequences_per_isolate, -genome_fragment, -order_within_fragment, 
         -accessory_fragment, -accessory_order_with_fragment, -qc, 
         -min_group_size_nuc, -max_group_size_nuc, -avg_group_size_nuc) %>%
  as.data.frame() %>%
  column_to_rownames("gene")

# replace gene name with 1, NA with 0
roary[] <- lapply(roary, function(x) as.integer(!is.na(x)))
roary <- t(roary)
roary <- as.data.frame(roary)

mers <- arrow::read_feather("outputs/sourmash_sketch_tables/s__Faecalibacterium-prausnitzii_D_k10_scaled1_wide.feather")
mers$acc <- gsub("_k10_scaled1", "", mers$acc)
mers$acc <- make_clean_names(mers$acc)
mers <- column_to_rownames(mers, var = "acc")


# correlate number of genes with number of kmers --------------------------
total_genes_per_genome <- rowSums(roary)
total_mers_per_genome <- rowSums(mers)

total_per_genome <- data.frame(total_genes_per_genome, total_mers_per_genome)
ggplot(total_per_genome, aes(x = total_genes_per_genome, y = total_mers_per_genome)) +
  geom_point() +
  theme_minimal() +
  labs(x = "genes per genome", y = "amino acid k-mers per genome (k = 10)", 
       title = "s__Faecalibacterium-prausnitzii_D",
       subtitle = "R-squared:  0.9605 (p < 0.001)")
summary(lm(total_genes_per_genome ~ total_mers_per_genome, data = total_per_genome))

# mantel test between presence/absence matrices ---------------------------

roary_dist <- dist(roary, method = "binary")
mers_dist <- dist(mers, method = "binary")
mantel_res <- mantel(roary_dist, mers_dist, method = "spearman", 
                     permutations = 9999, na.rm = TRUE)
mantel_res
# Mantel statistic r: 0.9782 
# Significance: 1e-04

# correlate number of unique genes with number of unique kmers ------------

# filter roary to columns that sum to 1, aka genes that are unique
roary_unique <- roary[ , colSums(roary == 1) == 1]
roary_unique <- rowSums(roary_unique)
mers_unique <- mers[ , colSums(mers == 1) == 1]
mers_unique <- rowSums(mers_unique)

unique_per_genome <- data.frame(roary_unique, mers_unique)
ggplot(unique_per_genome, aes(x = roary_unique, y = mers_unique)) +
  geom_point() +
  theme_minimal() +
  labs(x = "unique genes per genome", y = "unique amino acid k-mers per genome (k = 10)", 
       title = "s__Faecalibacterium-prausnitzii_D",
       subtitle = "R-squared:  0.9872 (p < 0.001)")
summary(lm(roary_unique ~ mers_unique, data = unique_per_genome))

# do the rarefaction curves match? ----------------------------------------


sp_roary <- specaccum(roary, "exact")
sp_roary_df <- data.frame(genomes = sp_roary$sites, genes = sp_roary$richness, sd = sp_roary$sd) 

ggplot() +
  geom_ribbon(data = sp_roary_df, aes(x = genomes, y = genes, ymin = genes - sd, ymax = genes + sd), fill = "grey", alpha = 1/3) +
  geom_point(data = sp_roary_df, aes(x = genomes, y = genes), size = 1) + 
  theme_minimal() +
  labs(title = "s__Faecalibacterium-prausnitzii_D genes")

sp_mers <- specaccum(mers, "exact")
sp_mers_df <- data.frame(genomes = sp_mers$sites, mers = sp_mers$richness, sd = sp_mers$sd) 

ggplot() +
  geom_ribbon(data = sp_mers_df, aes(x = genomes, y = mers, ymin = mers - sd, ymax = mers + sd), fill = "grey", alpha = 1/3) +
  geom_point(data = sp_mers_df, aes(x = genomes, y = mers), size = 1) + 
  theme_minimal() +
  labs(title = "s__Faecalibacterium-prausnitzii_D amino acid kmers (k = 10)")

library(BiodiversityR)
tmp <- accumresult(roary, method = "exact")

tmp2 <- accumresult(mers, "exact")
