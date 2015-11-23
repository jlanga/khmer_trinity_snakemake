shell.prefix("set -euo pipefail;")
configfile: "config.yaml.example"

# List variables
ENDS = "1 2 u".split()

# Folder variables
RAW_DIR      = "data/fastq_raw"
TRIM_DIR     = "data/fastq_trimmed"
NORM_DIR     = "data/fastq_norm"
ASSEMBLY_DIR = "data/assembly"


# Path to programs (or element on path)

trinity     = config["software"]["trinity"]
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]

rule all:
    input:
        expand(
            TRIM_DIR + "/{sample}_{end}.fastq.gz",
            sample = config["samples_pe"],
            end = ENDS
        )

rule clean:
    """
    rm -rf data/fastq_trimmed
    rm -rf data/fastq_norm
    rm -rf data/fastq_assembly
    """



rule trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and remove
    low quality regions and reads.
    """
    input:
        forward = lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward     = TRIM_DIR + "/{sample}_1.fastq.gz",
        reverse     = TRIM_DIR + "/{sample}_2.fastq.gz",
        unpaired    = TRIM_DIR + "/{sample}_u.fastq.gz"
    params:
        unpaired_1  = TRIM_DIR + "/{sample}_3.fastq.gz",
        unpaired_2  = TRIM_DIR + "/{sample}_4.fastq.gz",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"]
    benchmark:
        "benchmarks/trimmomatic/{sample}.json"
    log:
        "logs/trimmomatic/{sample}.log" 
    threads:
        4 # It doesn't perform well above this value
    shell:
        """
        {trimmomatic} PE                            \
            -threads 4                              \
            -{params.phred}                         \
            {input.forward}                         \
            {input.reverse}                         \
            {output.forward}                        \
            {params.unpaired_1}                     \
            {output.reverse}                        \
            {params.unpaired_2}                     \
            ILLUMINACLIP:{params.adaptor}:2:30:10   \
            2> {log}
            
        zcat {params.unpaired_1} {params.unpaired_2}    |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """



# rule trimmomatic_se
