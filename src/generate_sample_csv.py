#!/usr/bin/env python3

'''
Quick script to glob the files in read_base, match them to the index column in
sample_csv and write the paths to output_csv_file
'''

from pathlib import Path
import pandas
import snakemake

output_csv_file = 'data/samples.csv'
sample_csv = 'data/OG3550-P2 sample info - SAMPLE TEMPLATE .csv'
read_base = Path('data/OG3550P2-209261052',
                 'FASTQ_Generation_2020-11-09_21_26_03Z-339796457')

sample_data = pandas.read_csv(
    sample_csv,
    skiprows=48,
    header=None,
    index_col=0)

all_samples = sorted(set(sample_data.index))

# generate csv
with open(output_csv_file, 'wt') as f:
    f.write('sample,r1_path,r2_path\n')
    for x in all_samples:
        my_r1_path = Path(
            read_base,
            f'{x}_{{lane1}}-ds.{{hash}}',
            f'{x}_{{S}}_{{lane2}}_R1_{{lane3}}.fastq.gz')
        my_r2_path = Path(
            read_base,
            f'{x}_{{lane1}}-ds.{{hash}}',
            f'{x}_{{S}}_{{lane2}}_R2_{{lane3}}.fastq.gz')
        my_r1 = snakemake.io.glob_wildcards(my_r1_path)
        my_r2 = snakemake.io.glob_wildcards(my_r2_path)
        my_r1_files = snakemake.io.expand(my_r1_path, zip, **my_r1._asdict())
        my_r2_files = snakemake.io.expand(my_r2_path, zip, **my_r2._asdict())
        for i in range(0, len(my_r1_files)):
            f.write(f'{x},{my_r1_files[i]},{my_r2_files[i]}\n')
