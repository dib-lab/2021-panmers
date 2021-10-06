library(dplyr) #
library(readr) #
library(janitor) 
library(tibble) #
library(arrow)
library(ggplot2) #
library(vegan)
library(broom)

# functions ---------------------------------------------------------------

read_roary_presence_absence <- function(path){
  roary <- read_csv(path) %>%
    clean_names() %>%
    select(-non_unique_gene_name, -annotation, -no_isolates, -no_sequences, 
           -avg_sequences_per_isolate, -genome_fragment, -order_within_fragment, 
           -accessory_fragment, -accessory_order_with_fragment, -qc, 
           -min_group_size_nuc, -max_group_size_nuc, -avg_group_size_nuc) %>%
    as.data.frame() %>%
    column_to_rownames("gene")
  # replace gene name with 1, NA with 0
  roary[] <- lapply(roary, function(x) as.integer(!is.na(x)))
  roary <- as.data.frame(t(roary)) # switch orientation genome (rows) x genes (columns)
}

read_mers_presence_absence <- function(path) {
  mers <- arrow::read_feather(path)
  mers$acc <- gsub("_k10_scaled1", "", mers$acc)
  mers$acc <- make_clean_names(mers$acc)
  mers <- as.data.frame(mers)
  rownames(mers) <- mers$acc
  mers <- mers[ , -ncol(mers)]
}

correlate_total_per_genome <- function(roary, mers, species){
  
  total_genes_per_genome <- rowSums(roary)
  total_mers_per_genome <- rowSums(mers)
  
  total_per_genome <- data.frame(total_genes_per_genome, total_mers_per_genome)
  
  lm_result <- lm(total_genes_per_genome ~ total_mers_per_genome, data = total_per_genome) %>%
    glance() %>%
    mutate(species = species)
  
  r_squared <- round(lm_result$r.squared, digits = 2)
  p_value <- ifelse(lm_result$p.value > 0.001, 
                    paste0("= ", round(lm_result$p.value, digits = 3)), 
                    "< 0.001")
  
  plt <- ggplot(total_per_genome, aes(x = total_genes_per_genome, y = total_mers_per_genome)) +
    geom_point() +
    theme_minimal() +
    labs(x = "genes per genome", y = "k-mers per genome (protein, k=10, scaled=1)", 
         title = species,
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  return(list(lm_result = lm_result, plt = plt))
}

correlate_unique_per_genome <- function(roary, mers, species) {
  # filter roary to columns that sum to 1, aka genes that are unique
  roary_unique <- roary[ , colSums(roary == 1) == 1]
  roary_unique <- rowSums(roary_unique)
  mers_unique <- mers[ , colSums(mers == 1) == 1]
  mers_unique <- rowSums(mers_unique)
  unique_per_genome <- data.frame(roary_unique, mers_unique)
  
  lm_result <- lm(roary_unique ~ mers_unique, data = unique_per_genome) %>%
    glance() %>%
    mutate(species = species)
  
  r_squared <- round(lm_result$r.squared, digits = 2)
  p_value <- ifelse(lm_result$p.value > 0.001, 
                    paste0("= ", round(lm_result$p.value, digits = 3)), 
                    "< 0.001")
  
  ggplot(unique_per_genome, aes(x = roary_unique, y = mers_unique)) +
    geom_point() +
    theme_minimal() +
    labs(x = "unique genes per genome", y = "unique k-mers per genome (protein, k=10, scaled=1)", 
         title = species,
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  return(list(lm_result = lm_result, plt = plt))
}

mantel_roary_and_mers <- function(roary, mers, species){
  roary_dist <- dist(roary, method = "binary")
  mers_dist <- dist(mers, method = "binary")
  mantel_res <- mantel(roary_dist, mers_dist, method = "spearman", 
                       permutations = 9999, na.rm = TRUE)
  mantel_res <- data.frame(species = species, 
                           statistic = mantel_res$statistic,
                           signif = mantel_res$signif)
  return(mantel_res)
}

# apply -------------------------------------------------------------------

roary <- read_roary_presence_absence(path = snakemake@input[["roary"]])
mers <- read_mers_presence_absence(path = snakemake@input[["mers"]])
species <- unlist(snakemake@wildcards[['species']])
# correlate number of genes with number of kmers --------------------------

total_per_genome <- correlate_total_per_genome(roary = roary, mers = mers, species = species)

write_tsv(total_per_genome$lm_result, snakemake@output[['genes']])
pdf(snakemake@output[['genes_pdf']], height = 3, width = 3)
total_per_genome$plt
dev.off()
# correlate number of unique genes with number of unique kmers ------------

unique_per_genome <- correlate_unique_per_genome(roary = roary, mers = mers, species = species)

write_tsv(unique_per_genome$lm_result, snakemake@output[['unique']])
pdf(snakemake@output[['unique_pdf']], height = 3, width = 3)
unique_per_genome$plt
dev.off()
# mantel test between presence/absence matrices ---------------------------

mantel_results <- mantel_roary_and_mers(roary = roary, mers = mers, species = species)
write_tsv(mantel_results, snakemake@output[['mantel']])