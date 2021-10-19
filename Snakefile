import pandas as pd
import csv
import re

metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
SPECIES = metadata['species_no_space'].tolist()
SCALED = [100]

# Snakemake will use the ALPHA_KSIZE wildcard from rule all to generate output file names
# Then, snakemake will back propagate the strings from the final file names to solve for
# the wildcards "alphabet" and "ksize" throughout the rest of the workflow. 
# The underscore for the chrs in the list ALPHA_KSIZE separates the alphabet string from 
# the ksize string, allowing snakemake to solve {alphabet}_{ksize} wildcard strings. 
# Therefore, the chrs in the ALPHA_KSIZE list also set the alphabet names as "dayhoff" and "protein".

# set constrained k sizes
protein_ksizes = [7, 8, 9, 10, 11]
dayhoff_ksizes = [13, 15, 17]
hp_ksizes      = [27, 31]
# combine
ALPHA_KSIZE  = expand('protein-k{k}', k=protein_ksizes)
ALPHA_KSIZE += expand('dayhoff-k{k}', k=dayhoff_ksizes)
ALPHA_KSIZE += expand('hp-k{k}',      k=hp_ksizes)


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
        expand("outputs/correlate_pan_units/{alpha_ksize}_scaled{scaled}/{species}_genes.tsv", alpha_ksize = ALPHA_KSIZE, scaled = SCALED, species = SPECIES)

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
        mem_mb = 10000
    threads: 8
    benchmark: "benchmarks/roary/{species}.txt"
    params: 
        outdir = lambda wildcards: 'outputs/roary/' + wildcards.species 
    shell:'''
    roary -e -n -f {params.outdir} -p {threads} -z {input}
    mv {params.outdir}_[0-9]*/* {params.outdir}/
    rmdir {params.outdir}_[0-9]*
    '''

rule sourmash_sketch_species_genomes:
    input: 'outputs/prokka/{species}/{acc}/{acc}.faa'
    output: 'outputs/sourmash_sketch/{alpha}-k{ksize}_scaled{scaled}/{species}/{acc}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_sig/{alpha}-k{ksize}_scaled{scaled}_{species}_{acc}.txt"
    shell:"""
    sourmash sketch protein -p k={wildcards.ksize},scaled={wildcards.scaled},{wildcards.alpha} -o {output} --name {wildcards.acc} {input}
    """

rule convert_signature_to_csv:
    input: 'outputs/sourmash_sketch/{alpha}-k{ksize}_scaled{scaled}/{species}/{acc}.sig'
    output: 'outputs/sourmash_sketch/{alpha}-k{ksize}_scaled{scaled}/{species}/{acc}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_csv/{alpha}-k{ksize}_scaled{scaled}_{species}_{acc}.txt"
    resources:
        mem_mb=2000
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long:
    input: 
        Checkpoint_GatherResults("outputs/sourmash_sketch/{{alpha}}-k{{ksize}}_scaled{{scaled}}/{{species}}/{acc}.csv")
    output: csv = "outputs/sourmash_sketch_tables/{alpha}-k{ksize}_scaled{scaled}/{species}_long.csv"
    conda: 'envs/r.yml'
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_tables_long/{alpha}-k{ksize}_scaled{scaled}_{species}.txt"
    resources:
        mem_mb=64000
    script: "scripts/sketch_csv_to_long.R"

rule make_hash_table_wide:
    input:  "outputs/sourmash_sketch_tables/{alpha}-k{ksize}_scaled{scaled}/{species}_long.csv"
    output: "outputs/sourmash_sketch_tables/{alpha}-k{ksize}_scaled{scaled}/{species}_wide.feather"
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_tables_wide/{alpha}-k{ksize}_scaled{scaled}_{species}.txt"
    resources:
        mem_mb=32000
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
        mers="outputs/sourmash_sketch_tables/{alpha}-k{ksize}_scaled{scaled}/{species}_wide.feather"
    output:
        genes="outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_genes.tsv",
        genes_pdf="outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_genes.pdf",
        unique="outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_unique.tsv",
        unique_pdf="outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_unique.pdf",
        specaccum_pdf = "outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_specaccum.pdf"
        mantel="outputs/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}/{species}_mantel.tsv",
    threads: 1
    resources: mem_mb=6000
    benchmark: "benchmarks/correlate_pan_units/{alpha}-k{ksize}_scaled{scaled}_{species}.txt"
    conda: "envs/r_cor_pan.yml"
    script: "scripts/correlate_pan_units.R"
