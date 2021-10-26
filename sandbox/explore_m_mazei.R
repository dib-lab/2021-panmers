library(dplyr)
library(readr)
library(tidyr)
library(broom)
library(tibble)
library(janitor)
library(ggplot2)
library(feather)
library(pagoo)
library(plotly)

# functions ---------------------------------------------------------------

read_long_sketch_table_as_pagoo <- function(path){
  sketch_table <- read_csv(path)
  sketch_table <- sketch_table %>%
    select(gene    = minhash, 
           org     = acc,
           cluster = minhash) %>%
    filter(org != "GCF_000979855.1")
  p <- pagoo(data = as.data.frame(sketch_table))
  return(p)
}

compare_summary_stats <- function(roary_pagoo, kmer_pagoo){
  all_summary <- roary$summary_stats %>%
    as.data.frame() %>%
    select(Category, roary = Number) %>%
    left_join(as.data.frame(kmer_pagoo$summary_stats), by = "Category") %>%
    select(Category, roary, kmers = Number)
  
  all_summary %>%
    pivot_longer(cols = -Category, names_to = "method", values_to = "units") %>%
    pivot_wider(id_cols = c(method), names_from = Category, values_from = units) %>%
    mutate(core_pct = Core/Total) %>%
    mutate(shell_pct = Shell/Total) %>%
    mutate(cloud_pct = Cloud/Total) %>%
    pivot_longer(cols = -method)
}

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
  #mers$acc <- gsub("_k10_scaled1", "", mers$acc)
  mers$acc <- make_clean_names(mers$acc)
  mers <- as.data.frame(mers)
  rownames(mers) <- mers$acc
  mers <- mers[ , -ncol(mers)]
}

correlate_total_per_genome <- function(roary, mers, species, alpha, ksize, scaled){
  roary = roary2
  mers = pk10_2
  species = "M. mazei"
  alpha = "protein" 
  ksize = 10
  scaled = 100
  plot_params <- paste0("(", alpha, ", k=", ksize, ", scaled=", scaled, ")")
  
  total_genes_per_genome <- rowSums(roary)
  total_mers_per_genome <- rowSums(mers)
  
  total_per_genome <- data.frame(total_genes_per_genome, total_mers_per_genome)
  
  total_per_genome <- total_per_genome %>%
    rownames_to_column("genome") %>%
    filter(genome != "gca_001315865_1")
  
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
    theme(axis.title = element_text(size = 8),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(x = "genes per genome", y = paste0("k-mers per genome ", plot_params), 
         title = species,
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  return(list(lm_result = lm_result, plt = plt))
}
# read data ---------------------------------------------------------------

roary <- roary_2_pagoo("outputs/roary/s__Methanosarcina-mazei/gene_presence_absence.csv")
roary$summary_stats

pk10 <- read_long_sketch_table_as_pagoo("outputs/sourmash_sketch_tables/protein-k10_scaled100/s__Methanosarcina-mazei_long.csv")
pk10$summary_stats

compare_summary_stats(roary, pk10)


# redo corr ---------------------------------------------------------------

pk10$pan_matrix
roary$pan_matrix

tmp <- correlate_total_per_genome(roary = roary$pan_matrix,
                                  mers = pk10$pan_matrix, 
                                  species = "M. mazei", 
                                  alpha = "protein", 
                                  ksize = 10, 
                                  scaled = 100)

roary2 <- read_roary_presence_absence("outputs/roary/s__Methanosarcina-mazei/gene_presence_absence.csv")
pk10_2 <- read_mers_presence_absence("outputs/sourmash_sketch_tables/protein-k10_scaled100/s__Methanosarcina-mazei_wide.feather")

tmp2 <- correlate_total_per_genome(roary = roary2,
                                  mers = pk10_2, 
                                  species = "M. mazei", 
                                  alpha = "protein", 
                                  ksize = 10, 
                                  scaled = 100)
ggplotly(tmp2$plt)
# remove gca_001315865_1 and try again
# Excluded from RefSeq: many frameshifted proteins
# Definition: the CDSs predicted by the NCBI Prokaryotic Genome Annotation Pipeline 
#             have a suspiciously high number of frameshifts. For any clade 
#             containing at least 10 good quality assemblies the cutoff is more 
#             than three standard deviations from average or 5% of annotated CDSs,
#             whichever is larger. For any clade containing less than 10 good 
#             quality assemblies the cutoff is more than 30% of total CDSs.


# other plots -------------------------------------------------------------

pk10$summary_stats
pk10$gg_barplot()
pk10$gg_curves()
pk10$gg_binmap()
pk10$gg_dist()
pk10$gg_pca()
pk10$gg_pie()
