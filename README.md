# Panmers

Pangenomes are comprised of genes, panproteomes are comprised of proteins, panmers are comprised of kmers.
This repository explores whether amino acid k-mers can be used as a substitute for gene or protein sequences for the construction and interpretation of microbial "pangenomes."

A pangenome is the entire set of genes from all strains (genomes) in a clade/group/species.
The core genome is shared by all (>90%) of genomes and accounts for general ecological and phenotypic properties of species.
In contrast, the accessory genes are present in one or more genomes. 
The pangenome reflects metabolic and ecological plasticity of a population.
Accessory genes facilitate adaptation to environmental changes.
Pangenomes can be open, where adding additional genomes adds more genes to the pangenome, or closed, where adding additional genomes does not add more genes to the pangenome.

To infer a pangenome, coding domain sequences that represent genes are predicted from genome assemblies from many representatives in a clade/group/species.
Traditionally, these genomes are derived from isolate sequencing, as isolate sequencing libraries should contain sequences from a single genome.
This practice works well when there are sufficient isolate genomes for the species of interest.
However, reference- and assembly-free pangenomic analysis methods are desirable for cases in which isolate sequences do not exist.
In particular, such methods would facilite analysis of genome bins (or genome query neighborhoods) derived from metagenome sequecing.
Pangenome analysis is not traditionally applied to metagenome-assembled genomes given that these assemblies are likely species composites derived from all closely-related organisms present in a microbial community at the time of sampling and sequencing. 
However, metapangenomics has recently risen in popularity as an analysis technique to reveal the comprehensive complement of genes present in a species in a certain environment, as well as the genes that are most commonly present in those environments.
While genome assembly of isolate sequences often produces genomes that are ~100% complete, therefore lending these genomes to successful pangenome analysis and inference, metagenome assembly often fails on a large percentage of the sequencing library. 
This complicates metapangenome analysis, as genes that are present in a sample may not be assembled, and therefore would be left unanalyzed in the metapangenome. 
If amino acid k-mers can be used as a substitue for gene sequences for the construction and interpretation of metapangenomes, assembly would no longer be a necessary prerequisite for metapangenome analysis.
This would improve metapangenome recall and inference.

## Getting started with this repo

```
conda env create --name panmers --file environment.yml
conda activate panmers

snakemake -j 16 --use-conda --rerun-incomplete --latency-wait 15 --resources mem_mb=200000 --cluster "sbatch -t 10080 -J bsub -p bmm -n 1 -N 1 -c {threads} --mem={resources.mem_mb}" -k -n
```
