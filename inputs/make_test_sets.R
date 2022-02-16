library(dplyr)
library(readr)
setwd("~/github/2021-panmers/")

# This file was output by the 2020-ibd pipeline. 
species <- read_csv("inputs/gather_vita_vars_gtdb_shared_assemblies.x.genbank.gather.csv") %>%
  select(name) %>%
  mutate(accession = gsub(" .*", "", name))

lineages <- read_csv("https://osf.io/p6z3w/download") 

species <- left_join(species, lineages, by = c("accession" = "ident"))
ibd_species <- species$species

# define controlled species sets for panmers ------------------------------

# preferentially select species that are relevant to IBD.
# also choose species sets with many representatives (>=20)

ibd_pangenomes <- lineages %>%
  filter(species %in% ibd_species) %>%
  group_by(superkingdom, phylum, class, order, family, genus, species) %>%
  tally() %>%
  filter(n >= 20)

other_pangenomes <- lineages %>%
  group_by(superkingdom, phylum, class, order, family, genus, species) %>%
  tally() %>%
  filter(n >= 20) %>%
  filter(n < 1000) %>% # temporarily remove very large pangenomes for testing.
  arrange(desc(n))

other_pangenomes_slice <- other_pangenomes %>%
  group_by(phylum) %>%
  dplyr::slice(1)

test_pangenomes <- bind_rows(ibd_pangenomes, other_pangenomes_slice) %>%
  mutate(species_no_space = gsub(" ", "-", species))

write_tsv(test_pangenomes, "inputs/metadata.tsv")

# add pangenomes that were explored in the metapangenomes section of the paper
# B. uniformis and E. bolteae are already included; add B. fragilis, P. distasonis, P. merdae, and P. vulgatus

other_pangenomes_metap <- other_pangenomes %>%
  filter(species %in% c("s__Bacteroides fragilis", "s__Parabacteroides distasonis", 
                        "s__Parabacteroides merdae", "s__Phocaeicola vulgatus"))

test_pangenomes <- bind_rows(ibd_pangenomes, other_pangenomes_slice, other_pangenomes_metap) %>%
  mutate(species_no_space = gsub(" ", "-", species))

write_tsv(test_pangenomes, "inputs/metadata.tsv")

# create a test set to test the pipeline on. This will ensure the pipeline
# is reproducible even if I link in the desired genomes from the ones
# Tessa already has downloaded on farm (/home/ntpierce/2021-rank-compare/genbank/genomes)
test_pangenomes_small <- test_pangenomes %>%
  arrange(n) %>%
  ungroup() %>%
  group_by(superkingdom) %>%
  dplyr::slice(1)

write_tsv(test_pangenomes_small, "inputs/metadata_small.tsv")
