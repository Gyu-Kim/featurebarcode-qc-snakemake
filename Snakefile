# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import re
import gzip
import glob
import csv
import pandas as pd
from snakemake.utils import validate, min_version
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


##### load config and sample sheets #####


configfile: "config.yaml"
# report: "report/workflow.rst"

validate(config, schema="schemas/config.schema.yaml")

if config['samplesheet'] != "":
    samples = pd.read_csv(config['samplesheet']).set_index("Name", drop=False)
    validate(samples, schema="schemas/samples.schema.yaml")
    
    FEATURE_BC_IDS = dict(zip(samples["Name"], samples["Feature_BC_ID"]))
    SAMPLE_BC_IDS = dict(zip(samples["Name"], samples["Sample_BC_ID"])) if "Sample_BC_ID" in samples.columns else {}
else:
    samples = None
    FEATURE_BC_IDS = {}



# LIBRARY_FASTA = bowtie
# LIBRARY_BASENAME = "bowtie/" + os.path.splitext(os.path.basename(LIBRARY_FASTA))[0]
# LIBRARY_INDEX_DONEFILE = LIBRARY_BASENAME + ".done"

THREADS = config['threads']
FASTQ_DIR = config['fastq_dir']
SAMPLE_BC_FASTQ_DIR = config['sample_bc_fastq_dir'] if 'sample_bc_fastq_dir' in config.keys() else FASTQ_DIR



# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"


# include: "rules/other.smk"

wildcard_constraints:
    directory=".+\/",
    sample="[^\/]+"


localrules: all, combine_feature_counts, merge_fastqs, merge_fastqs_sample_bc, sample_bc_flatfile, make_library_fasta


rule all:
    input:
        "outs/feature_counts.txt",
        "outs/read_stats.csv",
        "outs/featurebarcode-qc-report.html"
    run:
        print("workflow complete!")

def get_report_deps(wildcards):
    inputs = {'counts' : "outs/feature_counts.txt",
              'stats' : 'outs/read_stats.csv'}
    # if SAMPLE_BC_IDS:
        # inputs['sample_bc_counts'] = "outs/sample_bc_counts.txt"
    return inputs

rule create_report:
    input: unpack(get_report_deps)
    output:
        "outs/featurebarcode-qc-report.html"
    script:
        "report/featurebarcode-qc-report.Rmd"

rule get_read_stats:
    input:
        trim = expand("outs/trim/{sample}_R1.fastq.gz", sample=FEATURE_BC_IDS.keys()),
        alns = expand("outs/alns/{sample}.bam", sample=FEATURE_BC_IDS.keys()),
        counts = expand("outs/feature_counts/{sample}.txt", sample=FEATURE_BC_IDS.keys()),
        # pdna_trim = "outs/pdna/trim/pDNA.fastq.gz",
        # pdna_alns = "outs/pdna/alns/pDNA.bam",
        pdna_counts = "outs/pdna/feature_counts/pDNA.txt"
    output:
        "outs/read_stats.csv"
    params:
        samplesheet = lambda wildcards: "--samplesheet " + config['samplesheet'] if samples is not None else "",
        pdna_fastq = lambda wildcards: "--pdna-fastq " + config['pdna_fastq'] if config['pdna_fastq'] is not "" else "",
        fastq_dir = lambda wildcards: "--fastq-dir " + config['fastq_dir'] if config['fastq_dir'] is not "" else ""
    shell:
        "python scripts/read_stats.py {output} {params.samplesheet} "
        "{params.pdna_fastq} {params.fastq_dir}"

rule combine_feature_counts:
    input:
        expand("outs/feature_counts/{sample}.txt",sample=FEATURE_BC_IDS),
        "outs/pdna/feature_counts/pDNA.txt"
    output: "outs/feature_counts.txt"
    params:
        input_string = lambda wildcards, input: ','.join(input)
        # input_string = ','.join(expand("outs/feature_counts/{sample}.txt", sample=SAMPLES.keys()) +
        #     ["outs/pdna/feature_counts/pDNA.txt"])
    shell:
        "python scripts/combine_feature_counts.py {params.input_string} {output}"



    

rule extract_umi_pdna:
    input: config['pdna_fastq']
    output: "outs/pdna/umi/pDNA.fastq.gz"
    shell:
        "umi_tools extract --extract-method string "
        "--bc-pattern NNNNNNNNNN "
        "--stdin {input} --stdout {output}"

rule trim_reads_pdna:
    input: "outs/pdna/umi/pDNA.fastq.gz"
    output: "outs/pdna/trim/pDNA.fastq.gz"
    threads: THREADS
    params:
        u6_promoter = config['trimming']['u6_promoter'],
        sgrna_scaffold = config['trimming']['sgrna_scaffold'],
        error_rate = config['trimming']['error_rate']
        min_len = config['trimming']['min_len']
    shell:
        "cutadapt -j {threads} -m {params.min_len} --discard-untrimmed "
        "-g \"{params.u6_promoter}...{params.sgrna_scaffold};max_error_rate={params.error_rate}\" "
        "-o {output} {input}"




def check_pdna_fastq(wildcards):
    if (config['pdna_fastq'] is not ""):
        return {'bam': "outs/pdna/alns/pDNA.bam", 'bamidx': "outs/pdna/alns/pDNA.bam.bai"}
    else:
        return {'feature_ref': config['feature_ref']}


rule feature_counts_pdna:
    input: unpack(check_pdna_fastq)
    output:
        counts="outs/pdna/feature_counts/pDNA.txt"
    log:
        log="outs/pdna/feature_counts/pDNA.log"
    run:
        if (hasattr(input, 'bam')):
            os.system("umi_tools count --per-contig --method " + config['dedup_method'] +
                " --stdin=" + input.bam + " --stdout=" + output.counts + " --log=" + log.log)
        else:
            with open(output.counts, mode='w') as output_file:
                output_csv = csv.writer(output_file, delimiter='\t')
                output_csv.writerow(['gene', 'count'])
                with open(input.feature_ref,mode='r') as input_file:
                    input_csv = csv.DictReader(input_file)
                    for row in input_csv:
                        output_csv.writerow([row['id'], '1'])
                    input_file.close()
                    output_file.close()




rule bowtie_align_pdna:
    input:
        donefile = "bowtie/feature_ref.done",
        fastq = "outs/pdna/trim/pDNA.fastq.gz"
    output:
        sam = "outs/pdna/alns/pDNA.sam"
    threads: THREADS
    params:
        index = "bowtie/feature_ref"
    shell:
        "bowtie --sam -v 1 -y -a --best -t -p {threads} {params.index} {input.fastq} {output.sam}"

def get_paired_fqs_sample_bc(wildcards):
    sample_bc_id = SAMPLE_BC_IDS[wildcards.sample]
    r1 = glob.glob(os.path.join(SAMPLE_BC_FASTQ_DIR, "**", sample_bc_id + "_*R1_*.fastq.gz"),
        recursive=True)
    r2 = glob.glob(os.path.join(SAMPLE_BC_FASTQ_DIR, "**", sample_bc_id + "_*R2_*.fastq.gz"), 
        recursive=True)
    if len(r1) == 0:
        raise ValueError(sample_bc_id + "has no matching sample barcode fastq file")
    if len(r1) != len(r2):
        raise ValueError(sample_bc_id + "has more than one matching sample barcode fastq file")
    return {"read1": sorted(r1), "read2": sorted(r2)}

rule merge_fastqs_sample_bc:
    input: unpack(get_paired_fqs_sample_bc)
    output:
        read1 = temp("outs/sample_bc/umi/{sample}_merged_R1.fastq.gz"),
        read2 = temp("outs/sample_bc/umi/{sample}_merged_R2.fastq.gz")
    shell:
        "cat {input.read1} > {output.read1}; "
        "cat {input.read2} > {output.read2}; "


rule extract_umi_sample_bc:
    input:
        read1 = "outs/sample_bc/umi/{sample}_merged_R1.fastq.gz",
        read2 = "outs/sample_bc/umi/{sample}_merged_R2.fastq.gz"
    output: "outs/sample_bc/umi/{sample}.fq.gz"
    params:
        whitelist = config['cell_barcode']['whitelist']
    shell:
        "umi_tools extract --extract-method regex --read2-stdout "
        "--bc-pattern '(?P<cell_1>.{{16}})(?P<umi_1>.{{12}})' " 
        "--bc-pattern2 '.{{8}}(?P<discard_1>.*)' "
        "--stdin {input.read1} --read2-in {input.read2} --stdout {output} "
        "--filter-cell-barcode --whitelist  {params.whitelist}"


rule sample_bc_flatfile:
    input: fq="outs/sample_bc/umi/{sample}.fq.gz"
    output: tsv="outs/sample_bc/umi/{sample}.tsv"
    run:
        o = open(output.tsv, 'w')
        fq = gzip.open(input.fq, 'rt')
        for bc in SeqIO.parse(fq, 'fastq'):
            # UMI and Cell BC need to be switched for current version of UMI-tools
            # extract appends _CB_UMI to the read ID
            # whereas count_tab command expects format _UMI_CB
            read_id = re.search("^(.*)_([ACGTN]+)_([ACGTN]+)$", bc.id)
            new_id = read_id.group(1) + "_" + read_id.group(3) + "_" + read_id.group(2)
            o.write(new_id + '\t' + str(bc.seq) + '\n')
        o.close()
        os.system("sort -k2,2 -o " + output.tsv + " " + output.tsv)


rule sample_bc_counts:
    input:
        tsv="outs/sample_bc/umi/{sample}.tsv"
    output:
        counts="outs/sample_bc/{sample}.txt",
        log="outs/sample_bc/{sample}.log"
    shell:
        "umi_tools count_tab --per-cell --method {config[dedup_method]} "
        "--stdin={input.tsv} --stdout={output.counts} --log={output.log}"

def get_paired_fqs(wildcards):
    sample_id = FEATURE_BC_IDS[wildcards.sample]
    r1 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R1_*.fastq.gz"),
        recursive=True)
    r2 = glob.glob(os.path.join(FASTQ_DIR, "**", sample_id + "_*R2_*.fastq.gz"), 
        recursive=True)
    if len(r1) == 0:
        raise ValueError(sample_id + " has no matching input fastq file")
    if len(r1) != len(r2):
        raise ValueError(sample_id + " has different numbers of R1 and R2 fastq files")
    return {"read1": sorted(r1), "read2": sorted(r2)}

rule merge_fastqs:
    input: unpack(get_paired_fqs)
    output:
        read1 = temp("outs/trim/{sample}_merged_R1.fastq.gz"),
        read2 = temp("outs/trim/{sample}_merged_R2.fastq.gz")
    shell:
        "cat {input.read1} > {output.read1}; "
        "cat {input.read2} > {output.read2}; "

rule trim_reads:
    input:
        read1 = "outs/trim/{sample}_merged_R1.fastq.gz",
        read2 = "outs/trim/{sample}_merged_R2.fastq.gz"
    output:
        read1 = "outs/trim/{sample}_R1.fastq.gz",
        read2 = "outs/trim/{sample}_R2.fastq.gz"
    threads: THREADS
    params:
        tso = config['trimming']['tso'],
        sgrna_scaffold = config['trimming']['sgrna_scaffold'],
        error_rate = config['trimming']['error_rate']
    shell:
        "cutadapt -j {threads} --discard-untrimmed -m 18 --pair-filter=first "
        "-g \"{params.tso}...{params.sgrna_scaffold};max_error_rate={params.error_rate}\" "
        "-o {output.read2} -p {output.read1} {input.read2} {input.read1}"

rule extract_umi:
    input:
        read1="outs/trim/{sample}_R1.fastq.gz",
        read2="outs/trim/{sample}_R2.fastq.gz"
    output:
        "outs/umi/{sample}.fastq.gz"
    params:
        whitelist = config['cell_barcode']['whitelist']
    shell:
        "umi_tools extract --extract-method regex --read2-stdout "
        "--bc-pattern '(?P<cell_1>.{{16}})(?P<umi_1>.{{12}})' "
        "--bc-pattern2 '(?P<discard_1>.*).{{20}}' "
        "--stdin {input.read1} --read2-in {input.read2} --stdout {output} "
        "--filter-cell-barcode --whitelist  {params.whitelist}"


rule make_library_fasta:
    input:
        feature_ref = config['feature_ref']
    output:
        fasta = "bowtie/feature_ref.fa"
    run:
        feature_ref = pd.read_csv(input.feature_ref).set_index("id", drop=False)
        feature_seqs = [SeqRecord(Seq(s), id=i, description='')
            for s, i in zip(feature_ref['sequence'], feature_ref['id'])]
        with open (output.fasta, 'w') as fa:
            SeqIO.write(feature_seqs, fa, "fasta")


rule bowtie_build:
    input:
        fasta = "bowtie/feature_ref.fa"
    output:
        touch("bowtie/feature_ref.done")
    params:
        basename = "bowtie/feature_ref"
    shell:
        "bowtie-build {input.fasta} {params.basename}"

rule bowtie_align:
    input:
        donefile = "bowtie/feature_ref.done",
        fastq = "outs/umi/{sample}.fastq.gz"
    output:
        sam = "outs/alns/{sample}.sam"
    threads: THREADS
    params:
        index = "bowtie/feature_ref"
    shell:
        "bowtie --sam -v 1 -y -a --best -t -p {threads} {params.index} {input.fastq} {output.sam}"

rule feature_counts:
    input:
        bam="outs/alns/{sample}.bam",
        bamidx="outs/alns/{sample}.bam.bai",
    output:
        counts="outs/feature_counts/{sample}.txt",
        log="outs/feature_counts/{sample}.log"
    shell:
        "umi_tools count --per-contig --per-cell --method {config[dedup_method]} "
        "--stdin={input.bam} --stdout={output.counts} --log={output.log}"


rule sam_to_bam_sort:
    input:
        "{directory}{sample}.sam"
    output:
        "{directory}{sample}.bam"
    threads: THREADS
    shell:
        "samtools view -b --threads {threads} {input} | "
        "samtools sort -@ {threads}  -o {output}"

    
rule samtools_index:
    input:
        "{directory}{sample}.bam"
    output:
        "{directory}{sample}.bam.bai"
    threads: THREADS
    shell:
        "samtools index -@ {threads} {input}"
