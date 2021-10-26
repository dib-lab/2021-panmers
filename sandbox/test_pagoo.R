library(dplyr)
library(ggplot2)
library(pagoo)
library(readr)
library(tidyr)
# pagoo -------------------------------------------------------------------

# SEE MICROPAN FLUIDITY micropan::fluidity

df <- read_csv("outputs/sourmash_sketch_tables/s__Cuniculiplasma-divulgatum_k10_scaled100_long.csv")
df <- df %>%
  select(gene    = minhash, 
         org     = acc,
         cluster = minhash)
p2 <- pagoo(data = as.data.frame(df))
p2$summary_stats
p2$gg_barplot()
p2$gg_curves()
p2$gg_binmap()
p2$gg_dist()
p2$gg_pca()
p2$gg_pie()


# import from roary -------------------------------------------------------


roary <- roary_2_pagoo(gene_presence_absence_csv = "outputs/roary/s__Cuniculiplasma-divulgatum/gene_presence_absence.csv")
roary$gg_curves()
roary$summary_stats


# compare panmers to roary with pagoo -------------------------------------

all_summary <- roary$summary_stats %>%
  as.data.frame() %>%
  select(Category, roary = Number) %>%
  left_join(as.data.frame(p2$summary_stats), by = "Category") %>%
  select(Category, roary, kmers = Number)

all_summary %>%
  pivot_longer(cols = -Category, names_to = "method", values_to = "units") %>%
  pivot_wider(id_cols = c(method), names_from = Category, values_from = units) %>%
  mutate(core_pct = Core/Total) %>%
  mutate(shell_pct = Shell/Total) %>%
  mutate(cloud_pct = Cloud/Total) %>%
  pivot_longer(cols = -method)

