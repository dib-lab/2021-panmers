import pandas as pd
import csv
import re

metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SPECIES = metadata['species_no_space'].tolist()

SCALED = [1, 100]

class Checkpoint_GatherResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self, species):
        gather_csv = f'outputs/genbank/{species}.x.genbank.gather.csv'
        assert os.path.exists(gather_csv)

        genome_accs = []
        with open(gather_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name'].split(' ')[0]
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {gather_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of rule 'format_bsub_accessions'; 
        # this will trigger exception until that rule has been run.
        checkpoints.grab_species_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.species)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p

rule all:
    input:
         expand('outputs/roary/{species}/pan_genome_reference.fa', species = SPECIES),
         expand("outputs/sourmash_sketch_tables/{species}_k10_scaled{scaled}_wide.feather", species = SPECIES, scaled = SCALED)

checkpoint grab_species_accessions:
    input: 
        lineages="/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv",
        metadata="inputs/metadata.tsv"
    output: csv="outputs/genbank/{species}.x.genbank.gather.csv",
    conda: "envs/tidyverse.yml"
    resources:
        mem_mb = 4000
    threads: 1
    script: "scripts/grab_species_accessions.R"

#checkpoint format_species_accessions:
#    input: "outputs/genbank/{species}_acc.csv",
#    output: "outputs/genbank/{species}.x.genbank.gather.csv",
#    resources:
#        mem_mb = 4000
#    threads: 1
#    shell:"""
#    sed '1iname,lineage' {input} > {output}
#    """

# I don't think we'll need the lineages for anything -- these were
# generated for charcoal decontamination, but that's not a necessary
# step in this pipeline. Leave rule here for now just in case.
#rule generate_species_lineages:
#    input: "outputs/genbank/bsub_acc.csv"
#    output: "outputs/genbank/bsub_assemblies.x.genbank.lineages.csv",
#    resources:
#        mem_mb = 4000
#    threads: 1
#    shell:"""
#    sed 's/,/_genomic.fna.gz,/1' {input} > {output}
#    """

rule make_genome_grist_conf_file:
    input: 
        species=expand("outputs/genbank/{species}.x.genbank.gather.csv", species = SPECIES)
    output:
        conf="conf/genome-grist-conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        species_list = [re.sub("outputs/genbank/", "", x) for x in input.species]
        species_list = [re.sub(".x.genbank.gather.csv", "", x) for x in species_list]
        species_list = "\n- ".join(species_list)
        with open(output.conf, 'wt') as fp:
           print(f"""\
sample:
- {species_list}
outdir: outputs
metagenome_trim_memory: 1e9
""", file=fp)

# specifying make_sgc_conf as the target downloads the genomes of interest,
# circumventing a bug in genome grist that prevents using the target
# download gather genomes. The sgc conf file is a dummy file -- it will be
# written to outputs/sgc, but the conf file has the wrong catlas bases.
rule download_species_assemblies:
    input: 
        conf = "conf/genome-grist-conf.yml"
    output: "genbank_genomes/download_done.txt"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    touch {output}
    '''

rule ls_download_species_assemblies:
    input: "genbank_genomes/download_done.txt"
    output: "genbank_genomes/{acc}_genomic.fna.gz"
    resources: mem_mb = 500
    threads: 1
    shell:'''
    ls {output}
    '''

rule gunzip_species_genomes:
    input: "genbank_genomes/{acc}_genomic.fna.gz"
    output: "genbank_genomes/{acc}_genomic.fna"
    resources:
        mem_mb = 1000
    threads: 1
    shell:'''
    gunzip -c {input} > {output}
    '''

rule prokka_species_genomes:
    output: 
        ffn = 'outputs/prokka/{species}/{acc}/{acc}.ffn',
        faa = 'outputs/prokka/{species}/{acc}/{acc}.faa',
        gff = 'outputs/prokka/{species}/{acc}/{acc}.gff',
    input: 'genbank_genomes/{acc}_genomic.fna'
    conda: 'envs/prokka.yml'
    resources:
        mem_mb = 8000
    threads: 2
    params: 
        outdir = lambda wildcards: 'outputs/prokka/' + wildcards.species + "/" +  wildcards.acc
        #prefix = lambda wildcards: wildcards.acc,
    shell:'''
    prokka {input} --outdir {params.outdir} --prefix {wildcards.acc} --metagenome --force --locustag {wildcards.acc} --cpus {threads} --centre X --compliant
    '''

rule roary_species_genomes:
    input: Checkpoint_GatherResults('outputs/prokka/{{species}}/{acc}/{acc}.gff')
    output: 
        'outputs/roary/{species}/pan_genome_reference.fa',
        'outputs/roary/{species}/gene_presence_absence.csv' 
    conda: 'envs/roary.yml'
    resources:
        mem_mb = 64000
    threads: 8
    benchmark: "benchmarks/roary/{species}.txt"
    params: 
        outdir = lambda wildcards: 'outputs/roary/' + wildcards.species 
    shell:'''
    roary -e -n -f {params.outdir} -p {threads} -z {input}
    mv {params.outdir}_*/* {params.outdir}/
    rmdir {params.outdir}_*
    '''

rule sourmash_sketch_species_genomes:
    input: 'outputs/prokka/{species}/{acc}/{acc}.faa'
    output: 'outputs/sourmash_sketch/{species}/{acc}_k10_scaled{scaled}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000
    threads: 1
    benchmark: "benchmarks/sourmash_sketch/{species}_{acc}_scaled{scaled}.txt"
    shell:"""
    sourmash sketch protein -p k=10,scaled={wildcards.scaled} -o {output} --name {wildcards.acc} {input}
    """

rule convert_signature_to_csv:
    input: 'outputs/sourmash_sketch/{species}/{acc}_k10_scaled{scaled}.sig'
    output: 'outputs/sourmash_sketch/{species}/{acc}_k10_scaled{scaled}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    resources:
        mem_mb=2000
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long:
    input: 
        Checkpoint_GatherResults("outputs/sourmash_sketch/{{species}}/{acc}_k10_scaled{scaled}.csv")
    output: csv = "outputs/sourmash_sketch_tables/{species}_k10_scaled{scaled}_long.csv"
    conda: 'envs/r.yml'
    threads: 1
    resources:
        mem_mb=64000
    script: "scripts/sketch_csv_to_long.R"

rule make_hash_table_wide:
    input:  "outputs/sourmash_sketch_tables/{species}_k10_scaled{scaled}_long.csv"
    output: "outputs/sourmash_sketch_tables/{species}_k10_scaled{scaled}_wide.feather"
    threads: 1
    resources:
        mem_mb=300000
    run:
        import pandas as pd
        import feather
        
        tab = pd.read_csv(str(input), dtype = {"minhash" : "float64", "acc" : "object", "present" : "int64"})
        tab_wide=tab.pivot(index='acc', columns='minhash', values='present')
        tab_wide = tab_wide.fillna(0)
        tab_wide['acc'] = tab_wide.index
        tab_wide = tab_wide.reset_index(drop=True)
        tab_wide.columns = tab_wide.columns.astype(str)
        tab_wide.to_feather(str(output)) 

rule correlate_pan_units:
    input:
        roary="outputs/roary/{species}/gene_presence_absence.csv",
        mers="outputs/sourmash_sketch_tables/{species}_k10_scaled1_wide.feather"
    output:
        genes="outputs/correlate_pan_units/{species}_genes.tsv",
        genes_pdf="outputs/correlate_pan_units/{species}_genes.pdf",
        unique="outputs/correlate_pan_units/{species}_unique.tsv",
        unique_pdf="outputs/correlate_pan_units/{species}_unique.pdf",
        mantel="outputs/correlate_pan_units/{species}_mantel.tsv",
    threads: 1
    resources: mem_mb=16000
    conda: "envs/r_cor_pan.yml"
    script: "scripts/correlate_pan_units.R"
