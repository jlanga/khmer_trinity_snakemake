shell.prefix("set -euo pipefail;")

configfile: "config.yaml.example"

# List variables
ENDS = "1 2 u".split()
PAIRS = "pe se".split()
SAMPLES_PE = config["samples_pe"]


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
            NORM_DIR + "/{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        )

rule clean:
    shell:
        """
        rm -rf data/fastq_trimmed
        rm -rf data/fastq_norm
        rm -rf data/fastq_assembly
        rm -rf logs
        rm -rf benchmarks
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
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
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
            {params.trimmomatic_params}
            2> {log}
            
        zcat {params.unpaired_1} {params.unpaired_2}    |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """



# rule QC_trimmomatic_se # Not dealing with SE reads from the sequencer


rule QC_interleave_and_filter_pe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Pipe the inputs, interleave, filter the stream and compress.
    """
    input:
        forward= TRIM_DIR + "/{sample}_1.fq.gz",
        reverse= TRIM_DIR + "/{sample}_2.fq.gz"
    output:
        interleaved= temp(TRIM_DIR + "/{sample}.qc.pe.fq.gz")
    threads:
        4
    params:
            minimum_quality= config["fastx_params"]["minimum_quality"],
            minimum_percent_quality= config["fastx_params"]["minimum_percent_quality"],
            minimum_length= config["fastx_params"]["minimum_length"]
    log:
        "logs/qc/interleave_and_filter_{sample}.log"
    benchmark:
        "benchmarks/qc/interleave_and_filter_{sample}.json"
    shell:
        """
        ( python bin/interleave-reads.py \
            <({gzip} -dc {input.forward} | cut -f 1 -d " ") \
            <({gzip} -dc {input.reverse} | cut -f 1 -d " ") |
        ./bin/fastq_quality_filter \
            -Q33 \
            -q {params.minimum_quality} \
            -p {params.minimum_percent_quality} |
        {gzip} -9 > {output.interleaved} ) \
        2> {log}
        """



rule QC_extract_pe_and_se_reads:
    """
    Extract PE and SE reads from the QC_interleave_and_filter_pe rule.
    """
    input:
        interleaved= TRIM_DIR + "/{sample}.qc.pe.fq.gz"
    output:
        paired= protected(TRIM_DIR + "/{sample}.final.pe.fq.gz"),
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
    """
    Filter the SE set of reads from the QC_trimmomatic_pe rule and merge it with the SE reads from the QC_extract_pe_and_se_reads.
    """
    input:
        without_qc = TRIM_DIR + "/{sample}.se.fq.gz",
        with_qc    = TRIM_DIR + "/{sample}.qc.se.fq.gz"
    output:
        fastq      = protected(TRIM_DIR + "/{sample}.final.se.fq.gz")
    threads:
        1
    params:
        minimum_quality= config["fastx_params"]["minimum_quality"],
        minimum_percent_quality= config["fastx_params"]["minimum_percent_quality"],
        minimum_length= config["fastx_params"]["minimum_length"]
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
            -q {params.minimum_quality} \
            -p {params.minimum_percent_quality} |
         {gzip} -9 >> {output.fastq} ) \
         2> {log}
        """



rule diginorm_load_into_counting:
    input:
        fastqs = expand(
            TRIM_DIR + "/{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        )
    output:
        table = temp(NORM_DIR + "/diginorm_table.kh"),
        info  = temp(NORM_DIR + "/diginorm_table.kh.info")
    threads:
        24
    log:
        "logs/diginorm/load_into_counting.log"
    benchmark:
        "benchmarks/diginorm/load_into_counting.json"
    params:
        ksize=         config["diginorm_params"]["ksize"],
        max_tablesize= config["diginorm_params"]["max_tablesize"],
        n_tables=      config["diginorm_params"]["n_tables"]
    shell:
        """
        ./bin/load-into-counting.py \
            --ksize {params.ksize} \
            --n_tables {params.n_tables} \
            --max-tablesize {params.max_tablesize} \
            --threads {threads} \
            --no-bigcount \
            {output.table} \
            {input.fastqs} \
        2>  {log}
        """



rule diginorm_normalize_by_median_pe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = TRIM_DIR + "/{sample}.final.pe.fq.gz",
        table = NORM_DIR + "/diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "/{sample}.keep.pe.fq.gz")
    threads:
        24 # Block RAM
    params:
        cutoff = config["diginorm_params"]["cutoff"]
    log:
        "logs/diginorm/normalize_by_median_pe_{sample}.log"
    benchmark:
        "benchmarks/diginorm/normalize_by_median_pe_{sample}.json"
    shell:
        """
        ( ./bin/normalize-by-median.py \
            --paired \
            --loadgraph {input.table} \
            --cutoff {params.cutoff} \
            --output - \
            {input.fastq} |
        pigz -9 \
        > {output.fastq} ) \
        2> {log}
        """



rule diginorm_normalize_by_median_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = TRIM_DIR + "/{sample}.final.se.fq.gz",
        table = NORM_DIR + "/diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "/{sample}.keep.se.fq.gz")
    threads:
        24 # Block excessive RAM usage
    params:
        cutoff = config["diginorm_params"]["cutoff"]
    log:
        "logs/diginorm/normalize_by_median_se_{sample}.log"
    benchmark:
        "benchmarks/diginorm/normalize_by_median_se_{sample}.json"
    shell:
        """
        ( ./bin/normalize-by-median.py \
            --loadgraph {input.table} \
            --cutoff {params.cutoff} \
            --output - \
            {input.fastq} |
        {gzip} -9 \
        > {output.fastq} ) \
        2> {log}
        """



rule diginorm_filter_abund:
    """
    Removes erroneus k-mers.
    """
    input:
        fastq = NORM_DIR + "/{sample}.keep.{pair}.fq.gz",
        table = NORM_DIR + "/diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "/{sample}.abundfilt.{pair}.fq.gz")
    threads:
        24
    params:
        ""
    log:
        "logs/diginorm_filter_abund_{sample}_{pair}.log"
    benchmark:
        "benchmarks/diginorm_filter_abunt_{sample}_{pair}.json"
    shell:
        """
        ( ./bin/filter-abund.py \
            --variable-coverage \
            --threads {threads} \
            --output - \
            {input.table} \
            {input.fastq} |
        {gzip} -9 \
        > {output.fastq} ) \
        2> {log}
        """



rule diginorm_extract_paired_reads:
    input:
        fastq = NORM_DIR + "/{sample}.abundfilt.pe.fq.gz"
    output:
        fastq_pe = NORM_DIR + "/{sample}.final.pe.fq.gz",
        fastq_se = temp(NORM_DIR + "/{sample}.temp.pe.fq.gz")
    threads:
        1
    log:
        "logs/diginorm/extract_paired_reads_{sample}.log"
    benchmark:
        "benchmarks/diginorm/extract_paired_reads_{sample}.json"
    shell:
        """
        ./bin/extract-paired-reads.py \
            --output-paired {output.fastq_pe} \
            --output-single {output.fastq_se} \
            --gzip \
            {input.fastq} \
        2> {log}
        """



rule diginorm_merge_single_reads:
    input:
        from_norm=   NORM_DIR + "/{sample}.abundfilt.se.fq.gz",
        from_paired= NORM_DIR + "/{sample}.temp.pe.fq.gz"
    output:
        fastq = NORM_DIR + "/{sample}.final.se.fq.gz"
    threads:
        1
    log:
        "logs/diginorm/merge_single_reads_{sample}.log"
    benchmark:
        "benchmarks/diginorm/merge_single_reads_{sample}.json"
    shell:
        """
        cp {input.from_norm} {output.fastq}
        {gzip} -dc {input.from_paired} |
        {gzip} -9 >> {output.fastq} 
        """
