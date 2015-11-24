shell.prefix("set -euo pipefail;")

configfile: "config.yaml.example"

# List variables
ENDS = "1 2 u".split()
PAIRS = "pe se".split()

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
            TRIM_DIR + "/{sample}.final.{pairs}.fq.gz",
            sample = config["samples_pe"],
            pairs = PAIRS
        )

rule clean:
    shell:
        """
        rm -rf data/fastq_trimmed
        rm -rf data/fastq_norm
        rm -rf data/fastq_assembly
        """



rule QC_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and remove
    low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from Trimmomatic. Further on they are catted and compressed with gzip/pigz (level 9).
    """
    input:
        forward = lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward     = temp(TRIM_DIR + "/{sample}_1.fq.gz"),
        reverse     = temp(TRIM_DIR + "/{sample}_2.fq.gz"),
        unpaired    = temp(TRIM_DIR + "/{sample}.se.fq.gz")
    params:
        unpaired_1  = TRIM_DIR + "/{sample}_3.fq.gz",
        unpaired_2  = TRIM_DIR + "/{sample}_4.fq.gz",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"]
    benchmark:
        "benchmarks/qc/trimmomatic_pe_{sample}.json"
    log:
        "logs/qc/trimmomatic_pe_{sample}.log" 
    threads:
        4 # It doesn't perform well above this value
    shell:
        """
        {trimmomatic} PE \
            -threads {threads} \
            -{params.phred} \
            <({gzip} -dc {input.forward} ) \
            <({gzip} -dc {input.reverse} ) \
            >({gzip} -9 > {output.forward} ) \
            {params.unpaired_1} \
            >({gzip} -9 > {output.reverse} ) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            2> {log}
            
        zcat {params.unpaired_1} {params.unpaired_2}    |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """



# rule QC_trimmomatic_se # Not dealing with SE reads from the sequencer


rule QC_interleave_and_filter_pe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Pipe the inputs, interleave, filter the stream and compress
    """
    input:
        forward= TRIM_DIR + "/{sample}_1.fq.gz",
        reverse= TRIM_DIR + "/{sample}_2.fq.gz"
    output:
        interleaved= temp(TRIM_DIR + "/{sample}.qc.pe.fq.gz")
    threads:
        4
    params:
        ""
    log:
        "logs/qc/interleave_and_filter_{sample}.log"
    benchmark:
        "benchmarks/qc/interleave_and_filter_{sample}.json"
    shell:
        """
        ( python bin/interleave-reads.py \
            <({gzip} -dc {input.forward}) \
            <({gzip} -dc {input.reverse}) |
        ./bin/fastq_quality_filter \
            -Q33 \
            -q 30 \
            -p 50 |
        {gzip} -9 > {output.interleaved} ) \
        2> {log}
        """



rule QC_extract_pe_and_se_reads:
    input:
        interleaved= TRIM_DIR + "/{sample}.qc.pe.fq.gz"
    output:
        paired= TRIM_DIR + "/{sample}.final.pe.fq.gz",
        single= temp(TRIM_DIR + "/{sample}.qc.se.fq.gz")
    params:
        ""
    threads:
        1
    log:
        "logs/qc/extract_paired_reads_{sample}.log"
    benchmark:
        "benchmarks/qc/extract_paired_reads_{sample}.json"
    shell:
        """
        python bin/extract-paired-reads.py \
            --output-paired {output.paired} \
            --output-single {output.single} \
            --gzip \
            {input.interleaved} \
        2> {log}
        """



rule QC_filter_and_merge_se:
    input:
        without_qc = TRIM_DIR + "/{sample}.se.fq.gz",
        with_qc    = TRIM_DIR + "/{sample}.qc.se.fq.gz"
    output:
        fastq      = TRIM_DIR + "/{sample}.final.se.fq.gz"
    threads:
        1
    log:
        "logs/qc/filter_and_merge_se_{sample}.log"
    benchmark:
        "benchmarks/qc/filter_and_merge_pe_{sample}.json"
    shell:
        """
        cp {input.with_qc} {output.fastq}
        
        ( {gzip} -dc {input.without_qc} |
        ./bin/fastq_quality_filter \
            fastq_quality_filter \
            -Q33 \
            -q 30 \
            -p 50 |
         {gzip} -9 >> {output.fastq} ) \
         2> {log}
        """

