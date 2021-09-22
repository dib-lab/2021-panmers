import pandas as pd
import csv
metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SPECIES = metadata['species_no_space'].tolist()


class Checkpoint_GatherResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self):
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
        checkpoints.format_species_accessions.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs(w.species)

        p = expand(self.pattern, acc=genome_accs, **w)
        return p

rule all:
    input:

# TODO grep/python from metadata
# or write R script to do it
rule grab_species_accessions:
    input: "/group/ctbrowngrp/gtdb/gtdb-rs202.taxonomy.v2.csv"
    output: "outputs/genbank/{species}_acc.csv",
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    grep -i s__Bacillus {input} | grep -i subtilis > {output}
    """

checkpoint format_species_accessions:
    input: "outputs/genbank/{species}_acc.csv",
    output: "outputs/genbank/{species}.x.genbank.gather.csv",
    resources:
        mem_mb = 4000
    threads: 1
    shell:"""
    sed '1iname,lineage' {input} > {output}
    """

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

rule make_genome_grist_conf_files:
    input: "outputs/genbank/{species}.x.genbank.gather.csv"
    output:
        conf="conf/{species}-genome-grist-conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        with open(output.conf, 'wt') as fp:
           print(f"""\
sample:
- {wildcards.species}
outdir: outputs
metagenome_trim_memory: 1e9
""", file=fp)

    
       
# specifying make_sgc_conf as the target downloads the genomes of interest,
# circumventing a bug in genome grist that prevents using the target
# download gather genomes. The sgc conf file is a dummy file -- it will be
# written to outputs/sgc, but the conf file has the wrong catlas bases.
rule download_species_assemblies:
    input: 
        gather_grist = "outputs/genbank/{species}.x.genbank.gather.csv",
        conf = "conf/{species}-genome-grist-conf.yml"
    output: "genbank_genomes/{acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
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
    input: Checkpoint_GatherResults('outputs/prokka/{species}/{acc}/{acc}.gff')
    output: 'outputs/roary/{species}/pan_genome_reference.fa' 
    conda: 'envs/roary.yml'
    resources:
        mem_mb = 64000
    threads: 8
    params: 
        outdir = 'outputs/roary/' + wildcards.species 
    shell:'''
    roary -e -n -f {params.outdir} -p {threads} -z {input}
    '''
