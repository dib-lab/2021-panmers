library(readr)
library(dplyr)

# metadata <- read_tsv("inputs/metadata_small.tsv")
metadata <- read_tsv(snakemake@input[['metadata']])
lineages <- read_csv(snakemake@input[['lineages']])

# grab the species of interest
species_no_space_w <- unlist(snakemake@wildcards[['species']])
#species_no_space_w <- "s__Cuniculiplasma-divulgatum"
metadata_w <- metadata %>%
  filter(species_no_space == species_no_space_w)

# write a tsv file that contains the lineage information for each set of species
# specified in the metadata file

lineages_w <- lineages %>% 
  filter(species %in% metadata_w$species) %>%
  select(name = ident, superkingdom, phylum, class, order, family, genus, species)

write_csv(lineages_w, snakemake@output[["csv"]], col_names = F)