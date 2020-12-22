#!/usr/bin/env python3

from Bio import Seq
from pathlib import Path
import pandas
import snakemake


#############
# FUNCTIONS #
#############

# quick reverse complement with biopython
def rc(x):
    return(str(Seq.Seq(x).reverse_complement()))

def get_fastq_paths(wildcards):
    return(samples.loc[wildcards.sample_number].to_dict())

###########
# GLOBALS #
###########

database = 'data/its1-58024_dada2.fa'
sample_csv = 'data/samples.csv'

samples = pandas.read_csv(
    sample_csv,
    index_col='sample')
all_samples = sorted(set(samples.index))

its_f = 'ATGCGATACTTGGTGTGAAT'
its_r = 'GACGCTTCTCCAGACTACAAT'

# singularity
pipeline = '/home/tom/Projects/pollen-its-wrapper/img/pollen-its-wrapper.sif'


#########
# RULES #
#########


rule target:
    input:
        'output/030_dada2/taxa_with_counts.csv'

# run dada2 workflow
rule merge_taxa_with_counts:
    input:
        nochim = 'output/030_dada2/seqtab_nochim.Rds',
        taxa = 'output/030_dada2/taxa.Rds',
    output:
        taxa_with_counts = 'output/030_dada2/taxa_with_counts.csv'
    log:
        'output/logs/030_dada2/merge_taxa_with_counts.log'
    container:
        pipeline
    script:
        'src/merge_taxa_with_counts.R'

rule run_dada2:
    input:
        r1 = expand('output/020_preprocess/s{sample_number}_cutadapt_r1.fq.gz',
                    sample_number=all_samples),
        r2 = expand('output/020_preprocess/s{sample_number}_cutadapt_r2.fq.gz',
                    sample_number=all_samples),
        fa = database
    output:
        seqtab = 'output/030_dada2/seqtab.Rds',
        nochim = 'output/030_dada2/seqtab_nochim.Rds',
        taxa = 'output/030_dada2/taxa.Rds',
        errF_plot = 'output/030_dada2/errF.pdf',
        errR_plot = 'output/030_dada2/errR.pdf'
    log:
        'output/logs/030_dada2/run_dada2.log'
    threads:
        workflow.cores
    container:
        pipeline
    script:
        'src/run_dada2.R'


# use cutadapt as recommended
rule cutadapt:
    input:
        r1 = 'output/020_preprocess/s{sample_number}_filtn_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_filtn_r2.fq.gz'
    output:
        r1 = 'output/020_preprocess/s{sample_number}_cutadapt_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_cutadapt_r2.fq.gz'
    params:
        primer_f = its_f,
        primer_r = its_r,
        rcf = rc(its_f),
        rcr = rc(its_r)
    log:
        'output/logs/020_preprocess/cutadapt_{sample_number}.log'
    container:
        pipeline
    shell:
        'cutadapt '
        '-g {params.primer_f} '
        '-a {params.rcr} '
        '-G {params.primer_r} '
        '-A {params.rcf} '
        '-n 2 '
        '-o {output.r1} '
        '-p {output.r2} '
        '{input.r1} '
        '{input.r2} '
        '&> {log}'


# use dada2's default trimming
rule d2_filter_and_trim:
    input:
        unpack(get_fastq_paths)
    output:
        r1 = 'output/020_preprocess/s{sample_number}_filtn_r1.fq.gz',
        r2 = 'output/020_preprocess/s{sample_number}_filtn_r2.fq.gz'
    log:
        'output/logs/020_preprocess/d2-filter-and-trim_{sample_number}.log'
    container:
        pipeline
    script:
        'src/d2_filter_and_trim.R'
