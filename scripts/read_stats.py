

import os
import sys
import re
import glob
import pandas as pd
import subprocess


def create_sample_dict(samplesheet):

    samples = pd.read_csv(samplesheet).set_index("Feature_BC_ID", drop=False)
    
    sample_dict = {}

    for sample, sample_id in zip(samples["Name"], samples["Feature_BC_ID"]):

        sample_dict[sample] = {'id': sample_id}
    
    return sample_dict

def get_fastq_file(fastq_dir, sample_id):
    fq = glob.glob(os.path.join(fastq_dir, "**", sample_id + "_*R1_*.fastq.gz"),
        recursive=True)
    if len(fq) == 0:
        raise ValueError(sample_id + "has no matching input fastq file")
    if len(fq) > 1:
        raise ValueError(sample_id + "has more than one matching input fastq file")
    return fq[0]


def collect_dedup_count(counts_file):

    # n_deduped = 0
    
    # with open(counts_file, 'r') as c:
    #     header = c.readline()
    #     for ref in c:
    #         n_deduped += int(ref.rstrip().split('\t')[1])

    # c.close()

    counts = pd.read_csv(counts_file, sep='\t')

    n_deduped = counts['count'].sum()

    return n_deduped

def collect_alns_count(alns_file):
    
    flagstat_output = subprocess.run('samtools flagstat ' + alns_file,
        stdout=subprocess.PIPE, shell=True)

    n_mapped = int(re.search("(\d+) \+ \d+ mapped", flagstat_output.stdout.decode()).group(1))

    return n_mapped

def collect_fastq_count(fastq_file):

    n_lines = subprocess.run('gunzip -c ' + fastq_file + ' | wc -l',
        stdout=subprocess.PIPE, shell=True)

    n_reads = int(int(n_lines.stdout)/4)
    
    return n_reads


def write_counts_to_output(samples_stats, out_file):


    with open(out_file, 'w') as o:

        # print('name', 'id', 'reads', 'trimmed', 'mapped', *dedup_keys, sep=',', file=o)
        print('name', 'id', 'reads', 'trimmed', 'mapped', 'dedup', sep=',', file=o)

        for sample in samples_stats.keys():

            sample_stats = samples_stats[sample]

            print(sample, sample_stats['id'], sample_stats['reads'], sample_stats['trimmed'], sample_stats['mapped'],
                # *[sample_stats[d] if not sample_stats.get(d) is None else "" for d in dedup_keys],
                sample_stats['dedup'],
                sep=',', file=o)

    o.close()




def process_stats(args):
    if (args.samplesheet is None):
        samples_stats = {}
    else:
        samples_stats = create_sample_dict(args.samplesheet)

        # dedup_keys = []


        for sample in samples_stats.keys():

            fastq_file = get_fastq_file(args.fastq_dir, samples_stats[sample]['id'])
            trim_file = os.path.join(args.trim_dir, sample + "_R1.fastq.gz")
            alns_file = os.path.join(args.alns_dir, sample + ".bam")
            # counts_files = glob.glob(os.path.join(counts_dir, sample + "*.txt"))
            counts_file = os.path.join(args.counts_dir, sample + ".txt")
            

            samples_stats[sample]['reads'] = collect_fastq_count(fastq_file)
            samples_stats[sample]['trimmed'] = collect_fastq_count(trim_file)
            samples_stats[sample]['mapped'] = collect_alns_count(alns_file)
            samples_stats[sample]['dedup'] = collect_dedup_count(counts_file)

            # for count_file in counts_files:
            #     dedup_method = re.search(sample + "\.?(.*).txt", count_file).group(1)
            #     dedup_key = "dedup"
            #     if dedup_method != "":
            #         dedup_key += "_" + dedup_method
            #     if not dedup_key in dedup_keys: dedup_keys.append(dedup_key)
            #     samples_stats[sample][dedup_key] = collect_dedup_count(count_file)


    if (args.pdna_fastq is not None):
        samples_stats['pDNA'] = {"id": "pDNA"}
        samples_stats['pDNA']['reads'] = collect_fastq_count(args.pdna_fastq)
        samples_stats['pDNA']['trimmed'] = collect_fastq_count(args.pdna_trim)
        samples_stats['pDNA']['mapped'] = collect_alns_count(args.pdna_alns)
        samples_stats['pDNA']['dedup'] = collect_dedup_count(args.pdna_counts)

    
    write_counts_to_output(samples_stats, args.out_file)


if __name__ == "__main__":

##############################################################################
    # Initiates argument parser
    import argparse

    parser = argparse.ArgumentParser(description='Get read-level stats for feature barcode reads')


    parser.add_argument('out_file', help = 'Where to write output', type=str)

    
    # Optional arguments

    parser.add_argument('--samplesheet', help='Path to CSV with ID and Name header',
                        default = None, type=str)

    parser.add_argument('--pdna-fastq', help = 'Where to find pdna fastq',
                        default = None, type=str)

    parser.add_argument('--fastq-dir', help = 'Where to find demultiplexed fastqs',
                        default="fastqs/", type=str)


    parser.add_argument('--trim-dir', help = 'Where to find aligned bam files',
                        default="outs/trim/", type=str)

    parser.add_argument('--alns-dir', help = 'Where to find aligned bam files',
                        default="outs/alns/", type=str)

    parser.add_argument('--counts-dir', help = 'Where to find deduped count files',
                        default="outs/feature_counts/", type=str)

    
    parser.add_argument('--pdna-trim', default = "outs/pdna/trim/pDNA.fastq.gz", type=str)
    parser.add_argument('--pdna-alns', default = "outs/pdna/alns/pDNA.bam", type=str)
    parser.add_argument('--pdna-counts', default = "outs/pdna/feature_counts/pDNA.txt", type=str)

    args = parser.parse_args()

##############################################################################


    process_stats(args)






