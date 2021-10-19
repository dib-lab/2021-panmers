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
  #mers$acc <- gsub("_k10_scaled1", "", mers$acc)
  mers$acc <- make_clean_names(mers$acc)
  mers <- as.data.frame(mers)
  rownames(mers) <- mers$acc
  mers <- mers[ , -ncol(mers)]
}

correlate_total_per_genome <- function(roary, mers, species, alpha, ksize, scaled){
  plot_params <- paste0("(", alpha, ", k=", ksize, ", scaled=", scaled, ")")
  
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
    theme(axis.title = element_text(size = 8),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(x = "genes per genome", y = paste0("k-mers per genome ", plot_params), 
         title = species,
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  return(list(lm_result = lm_result, plt = plt))
}

correlate_unique_per_genome <- function(roary, mers, species, alpha, ksize, scaled){
  plot_params <- paste0("(", alpha, ", k=", ksize, ", scaled=", scaled, ")")

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
  
  plt <- ggplot(unique_per_genome, aes(x = roary_unique, y = mers_unique)) +
    geom_point() +
    theme_minimal() +
    theme(axis.title = element_text(size = 8),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(x = "unique genes per genome", y = paste0("unique k-mers per genome ", plot_params), 
         title = species,
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  return(list(lm_result = lm_result, plt = plt))
}

correlate_specaccum <- function(roary, mers, speciesalpha, ksize, scaled){
  plot_params <- paste0("(", alpha, ", k=", ksize, ", scaled=", scaled, ")")
  
  sp_roary <- specaccum(roary, "exact")
  sp_roary_df <- data.frame(genomes = sp_roary$sites, genes = sp_roary$richness, sd = sp_roary$sd) 
  
  sp_mers <- specaccum(mers, "exact")
  sp_mers_df <- data.frame(genomes = sp_mers$sites, mers = sp_mers$richness, sd = sp_mers$sd) 
  
  tmp1 <- data.frame(genomes = sp_roary$sites, units = sp_roary$richness, sd = sp_roary$sd, measure = "genes")
  tmp2 <-  data.frame(genomes = sp_mers$sites, units = sp_mers$richness, sd = sp_mers$sd, measure = "kmers")
  both_sp_long <- rbind(tmp1, tmp2)
  
  specaccum_plt <- ggplot() +
    geom_ribbon(data = both_sp_long, aes(x = genomes, y = units, ymin = units - sd, ymax = units + sd, fill = measure), alpha = 1/3) +
    geom_point(data = both_sp_long, aes(x = genomes, y = units, color = measure), size = 1) + 
    theme_minimal() +
    theme(plot.title.position = "plot",
          legend.position = "top",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(title = species, 
         subtitle = plot_params)
  
  both_sp_wide <- left_join(sp_mers_df, sp_roary_df, by = "genomes")
  lm_result <- lm(mers ~ genes, both_sp_wide) %>%
    glance()
  
  r_squared <- round(lm_result$r.squared, digits = 2)
  p_value <- ifelse(lm_result$p.value > 0.001, 
                    paste0("= ", round(lm_result$p.value, digits = 3)), 
                    "< 0.001")
  
  specaccum_cor_plt <- ggplot(both_sp_wide, aes(x = genes, y = mers)) +
    geom_point() +
    theme_minimal() +
    theme(plot.title.position = "plot", 
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(x = "genes", y = "k-mers", 
         title = "Accumulation curve values per genome added",
         subtitle = paste0("R-squared:  ", r_squared, " (p ", p_value, ")"))
  
  ggarrange(specaccum_plt, specaccum_cor_plt)
  
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

total_per_genome <- correlate_total_per_genome(roary = roary, mers = mers, species = species,
                                               alpha = snakemake@wildcards[['alpha']],
                                               ksize = snakemake@wildcards[['ksize']],
                                               scaled = snakemake@wildcards[['scaled']])

write_tsv(total_per_genome$lm_result, snakemake@output[['genes']])
pdf(snakemake@output[['genes_pdf']], height = 3, width = 3)
total_per_genome$plt
dev.off()

# correlate number of unique genes with number of unique kmers ------------

unique_per_genome <- correlate_unique_per_genome(roary = roary, mers = mers, species = species,
                                                 alpha = snakemake@wildcards[['alpha']],
                                                 ksize = snakemake@wildcards[['ksize']],
                                                 scaled = snakemake@wildcards[['scaled']])

write_tsv(unique_per_genome$lm_result, snakemake@output[['unique']])
pdf(snakemake@output[['unique_pdf']], height = 3, width = 3)
unique_per_genome$plt
dev.off()


# correlate specaccum -----------------------------------------------------

correlate_specaccum_plt <- correlate_specaccum(roary = roary, mers = mers, species = species,
                                               alpha = snakemake@wildcards[['alpha']],
                                               ksize = snakemake@wildcards[['ksize']],
                                               scaled = snakemake@wildcards[['scaled']])
pdf(snakemake@output[['specaccum_pdf']], height = 3, width = 6)
correlate_specaccum_plt
dev.off()

# mantel test between presence/absence matrices ---------------------------

mantel_results <- mantel_roary_and_mers(roary = roary, mers = mers, species = species)
write_tsv(mantel_results, snakemake@output[['mantel']])
