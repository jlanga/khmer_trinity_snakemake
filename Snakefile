shell.prefix("set -euo pipefail;")
configfile: "config.yaml"



# List variables
ENDS = "1 2 u".split()
PAIRS = "pe_pe pe_se".split()
SAMPLES_PE  = config["samples_pe"]
SAMPLES_SE  = config["samples_se"]



# Folder variables
RAW_DIR         = "data/fastq_raw"
FASTQC_RAW_DIR  = "data/fastqc_raw"
TRIM_DIR        = "data/fastq_trimmed"
FASTQC_TRIM_DIR = "data/fastqc_trimmed"
NORM_DIR        = "data/fastq_norm"
FASTQC_NORM_DIR = "data/fastqc_norm"
ASSEMBLY_DIR    = "data/assembly"
PART_DIR        = "data/partitions"



# Path to programs (or element on path)
trinity     = config["software"]["trinity"]
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]



# Special variables - To make sense on why use all
ALL_TRHEADS   = 999999 # In case you want to use all threads
BLOCK_THREADS = 999999 # In case you only want to block excessive RAM usage



rule all:
    input:
        PART_DIR + "/Trinity.fasta",
        expand( # fastqc zip and html for raw PE data
            FASTQC_RAW_DIR + "/{sample}_{end}.fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            FASTQC_RAW_DIR + "/{sample}.fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for PE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for PE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        )



rule clean:
    shell:
        """
        rm -rf data/fastqc_raw
        rm -rf data/fastq_trimmed
        rm -rf data/fastqc_trimmed
        rm -rf data/fastq_norm
        rm -rf data/fastqc_norm
        rm -rf data/assembly
        rm -rf data/partitions
        rm -rf logs
        rm -rf benchmarks
        """



rule raw_make_links_pe:
    input:
        forward= lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse= lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward= RAW_DIR + "/{sample}_1.fq.gz",
        reverse= RAW_DIR + "/{sample}_2.fq.gz"
    threads:
        1
    log:
        "logs/raw/make_links_pe_{sample}.log"
    benchmark:
        "benchmarks/raw/make_links_pe_{sample}.json"
    shell:
        """
        ln -rs {input.forward} {output.forward} 2> {log}
        ln -rs {input.reverse} {output.reverse} 2>> {log}
        """ 



rule raw_make_links_se:
    input:
        single= lambda wildcards: config["samples_se"][wildcards.sample]["single"],
    output:
        single= RAW_DIR + "/{sample}.fq.gz"
    threads:
        1
    log:
        "logs/raw/make_links_se_{sample}.log"
    benchmark:
        "benchmarks/raw/make_links_se_{sample}.json"
    shell:
        """
        ln -rs {input.single} {output.single} 2>  {log}
        """ 



rule raw_fastqc_sample:
    input:
        fastq = RAW_DIR + "/{sample}.fq.gz"
    output:
        html= FASTQC_RAW_DIR + "/{sample}.fastqc.html",
        zip=  FASTQC_RAW_DIR + "/{sample}.fastqc.zip"
    threads:
        1
    params:
        html= RAW_DIR + "/{sample}_fastqc.html",
        zip=  RAW_DIR + "/{sample}_fastqc.zip"
    log:
        "logs/raw/fastqc_{sample}.log"
    benchmark:
        "benchmarks/raw/fastqc_{sample}.json"
    shell:
        """
        fastqc \
            --nogroup \
            {input.fastq} \
        2> {log} 1>&2
        
        mv {params.html} {output.html}
        mv {params.zip}  {output.zip}
        """



rule raw_fastqc:
    input:
        expand(
            FASTQC_RAW_DIR + "/{sample}_{end}.fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = "zip html".split()
        ) + expand(
            FASTQC_RAW_DIR + "/{sample}.fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "zip html".split()
        )

####
# Quality control rules
####

rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from 
    Trimmomatic. Further on they are catted and compressed with gzip/pigz 
    (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
    header. It is done posterior to the trimming since the output comes 
    slower than the input is read.
    """
    input:
        forward = RAW_DIR + "/{sample}_1.fq.gz",
        reverse = RAW_DIR + "/{sample}_2.fq.gz"
    output:
        forward     = temp(TRIM_DIR + "/{sample}_1.fq.gz"),
        reverse     = temp(TRIM_DIR + "/{sample}_2.fq.gz"),
        unpaired    = protected(TRIM_DIR + "/{sample}.final.pe_se.fq.gz")
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
            >(cut -f 1 -d " " | {gzip} -9 > {output.forward} ) \
            {params.unpaired_1} \
            >(cut -f 1 -d " " | {gzip} -9 > {output.reverse} ) \
            {params.unpaired_2} \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}
            
        zcat {params.unpaired_1} {params.unpaired_2} |
        cut -f 1 -d " " |
        {gzip} -9 > {output.unpaired}
        
        rm {params.unpaired_1} {params.unpaired_2}
        """



rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and 
        remove low quality regions and reads.
    Input is piped through gzip/pigz.
    Output is piped to gzip.
    """
    input:
        single = RAW_DIR + "/{sample}.fq.gz",
    output:
        single     = protected(TRIM_DIR + "/{sample}.final.se.fq.gz")
    params:
        adaptor     = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/qc/trimmomatic_se_{sample}.json"
    log:
        "logs/qc/trimmomatic_se_{sample}.log" 
    threads:
        4 # It doesn't perform well above this value
    shell:
        """
        {trimmomatic} SE \
            -threads {threads} \
            -{params.phred} \
            <({gzip} -dc {input.single}) \
            >(cut -f 1 -d " " | {gzip} -9 > {output.single}) \
            ILLUMINACLIP:{params.adaptor}:2:30:10 \
            {params.trimmomatic_params} \
        2> {log}
        """
    


rule qc_interleave_pe_pe:
    """
    From the adaptor free _1 and _2 , interleave the reads.
    Read the inputs, interleave, filter the stream and compress.
    """
    input:
        forward= TRIM_DIR + "/{sample}_1.fq.gz",
        reverse= TRIM_DIR + "/{sample}_2.fq.gz"
    output:
        interleaved= protected(TRIM_DIR + "/{sample}.final.pe_pe.fq.gz")
    threads:
        2 # One for the pairer and other for gzip
    log:
        "logs/qc/interleave_pe_{sample}.log"
    benchmark:
        "benchmarks/qc/interleave_pe_{sample}.json"
    shell:
        """
        ( python bin/interleave-reads.py \
            {input.forward} \
            {input.reverse} |
        {gzip} -9 > {output.interleaved} ) \
        2> {log}
        """



rule qc_fastqc:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = TRIM_DIR + "/{sample}.final.{pair}.fq.gz"
    output:
        zip   = protected(FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.zip"),
        html  = protected(FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.html")
    threads:
        1
    params:
        zip   = TRIM_DIR + "/{sample}.final.{pair}_fastqc.zip",
        html  = TRIM_DIR + "/{sample}.final.{pair}_fastqc.html"
    log:
        "logs/qc/fastqc_trimmed_{sample}_{pair}.log"
    benchmark:
        "benchmarks/qc/fastqc_trimmed_{sample}_{pair}.json"
    shell:
        """
        fastqc \
            --nogroup \
            {input.fastq} \
        > {log} 2>&1
        
        mv {params.zip} {output.zip}
        mv {params.html} {output.html}
        """



rule qc:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    """
    input:
        expand( # fastqc zip and html for raw PE data
            FASTQC_RAW_DIR + "/{sample}_{end}.fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            FASTQC_RAW_DIR + "/{sample}.fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "zip html".split()
        ),
        expand( # fq.gz files for PE data
            TRIM_DIR + "/{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair   = PAIRS
        ),
        expand( # fq.gz files for SE data
            TRIM_DIR + "/{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_SE,
            pair = ["se"]
        ),
        expand( # fastqc zip and html files for PE data
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        )



###
# Diginorm rules
###

rule diginorm_load_into_counting:
    """
    Build the hash table data structure from all the trimmed reads.
    Caution!: The --help says that it can be multithreaded but it raises
    errors!
    """
    input:
        fastqs = expand(
            TRIM_DIR + "/{sample}.final.{pair}.fq.gz",
            sample = SAMPLES_PE,
            pair = PAIRS
        ) + expand(
            TRIM_DIR + "/{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        table = temp(NORM_DIR + "/diginorm_table.kh"),
        info  = temp(NORM_DIR + "/diginorm_table.kh.info")
    threads:
        1
    log:
        "logs/diginorm/load_into_counting.log"
    benchmark:
        "benchmarks/diginorm/load_into_counting.json"
    params:
        ksize=    config["diginorm_params"]["ksize"],
        hashsize= config["diginorm_params"]["hashsize"],
        n_hashes= config["diginorm_params"]["n_hashes"]
    shell:
        """
        ./bin/load-into-counting.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --no-bigcount \
            {output.table} \
            {input.fastqs} \
        >  {log} 2>&1
        """



rule diginorm_normalize_by_median_pe_pe:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = TRIM_DIR + "/{sample}.final.pe_pe.fq.gz",
        table = NORM_DIR + "/diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "/{sample}.keep.pe_pe.fq.gz")
    threads:
        BLOCK_THREADS # Block
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.pe_pe.fq.gz.keep"
    log:
        "logs/diginorm/normalize_by_median_pe_pe_{sample}.log"
    benchmark:
        "benchmarks/diginorm/normalize_by_median_pe_pe_{sample}.json"
    shell:
        """
        ./bin/normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --paired \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1
         
        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}
        
        rm {params.keep_fq}
        """



rule diginorm_normalize_by_median_pe_se:
    """
    Normalizes by median EACH FILE.
    Therefore one loads once per file the hash table.
    """
    input:
        fastq = TRIM_DIR + "/{sample}.final.pe_se.fq.gz",
        table = NORM_DIR + "/diginorm_table.kh"
    output:
        fastq = temp(NORM_DIR + "/{sample}.keep.pe_se.fq.gz")
    threads:
        BLOCK_THREADS # Block excessive RAM usage
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.pe_se.fq.gz.keep"
    log:
        "logs/diginorm/normalize_by_median_pe_se_{sample}.log"
    benchmark:
        "benchmarks/diginorm/normalize_by_median_pe_se_{sample}.json"
    shell:
        """
        ./bin/normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1
         
        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}
        
        rm {params.keep_fq}
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
        BLOCK_THREADS # Block excessive RAM usage
    params:
        cutoff   = config["diginorm_params"]["cutoff"],
        ksize    = config["diginorm_params"]["ksize"],
        n_hashes = config["diginorm_params"]["n_hashes"],
        hashsize = config["diginorm_params"]["hashsize"],
        keep_fq  = "{sample}.final.se.fq.gz.keep"
    log:
        "logs/diginorm/normalize_by_median_se_{sample}.log"
    benchmark:
        "benchmarks/diginorm/normalize_by_median_se_{sample}.json"
    shell:
        """
        ./bin/normalize-by-median.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --cutoff {params.cutoff} \
            --loadhash {input.table} \
            {input.fastq} \
        > {log} 2>&1
         
        pigz -9c \
            {params.keep_fq} \
        > {output.fastq} \
        2>> {log}
        
        rm {params.keep_fq}
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
        BLOCK_THREADS # BLOCK
    params:
        abundfilt_fq = "{sample}.keep.{pair}.fq.gz.abundfilt"
    log:
        "logs/diginorm/filter_abund_{sample}_{pair}.log"
    benchmark:
        "benchmarks/diginorm_filter_abunt_{sample}_{pair}.json"
    shell:
        """
        ./bin/filter-abund.py \
            --variable-coverage \
            {input.table} \
            {input.fastq} \
        > {log} 2>&1
            
        {gzip} -9c \
            {params.abundfilt_fq} \
        > {output.fastq}  \
        2>> {log}
        
        rm {params.abundfilt_fq}
        """



rule diginorm_extract_paired_reads:
    """
    Split the filtered reads into PE and SE.
    """
    input:
        fastq = NORM_DIR + "/{sample}.abundfilt.pe_pe.fq.gz"
    output:
        fastq_pe = protected(NORM_DIR + "/{sample}.final.pe_pe.fq.gz"),
        fastq_se = temp(NORM_DIR + "/{sample}.temp.pe_se.fq.gz")
    threads:
        1
    params:
        fastq_pe = "{sample}.abundfilt.pe_pe.fq.gz.pe",
        fastq_se = "{sample}.abundfilt.pe_pe.fq.gz.se"
    log:
        "logs/diginorm/extract_paired_reads_{sample}.log"
    benchmark:
        "benchmarks/diginorm/extract_paired_reads_{sample}.json"
    shell:
        """
        ./bin/extract-paired-reads.py \
            {input.fastq} \
        > {log} 2>&1
        
        {gzip} -9c {params.fastq_pe} > {output.fastq_pe}
        {gzip} -9c {params.fastq_se} > {output.fastq_se}
        
        rm {params.fastq_pe} {params.fastq_se}
        """



rule diginorm_merge_pe_single_reads:
    """
    Put together the SE reads from the same sample
    """
    input:
        from_norm=   NORM_DIR + "/{sample}.abundfilt.pe_se.fq.gz",
        from_paired= NORM_DIR + "/{sample}.temp.pe_se.fq.gz"
    output:
        fastq = protected(NORM_DIR + "/{sample}.final.pe_se.fq.gz")
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



rule dignorm_get_former_se_reads:
    """
    Move the result of diginorm_extract_paired_reads for true SE reads 
    to their final position.
    """
    input:
        single= NORM_DIR + "/{sample}.abundfilt.se.fq.gz"
    output:
        single= NORM_DIR + "/{sample}.final.se.fq.gz"
    threads:
        1
    log:
        "logs/diginorm/get_former_se_reads_{sample}.log"
    benchmark:
        "benchmarks/diginorm/get_former_se_reads_{sample}.json"
    shell:
        """
        mv {input.single} {output.single}
        """



rule diginorm_fastqc:
    """
    Do FASTQC reports
    Uses --nogroup!
    One thread per fastq.gz file
    """
    input:
        fastq = NORM_DIR + "/{sample}.final.{pair}.fq.gz"
    output:
        zip   = protected(FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.zip"),
        html  = protected(FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.html") 
    threads:
        1
    params:
        zip   = NORM_DIR + "/{sample}.final.{pair}_fastqc.zip",
        html  = NORM_DIR + "/{sample}.final.{pair}_fastqc.html"
    log:
        "logs/diginorm/fastqc_trimmed_{sample}_{pair}.log"
    benchmark:
        "benchmarks/diginorm/fastqc_trimmed_{sample}_{pair}.json"
    shell:
        """
        fastqc \
            --nogroup \
            {input.fastq} \
        > {log} 2>&1
        
        mv {params.zip} {output.zip}
        mv {params.html} {output.html}
        """



rule diginorm:
    """
    Rule to do the Quality Control and the Digital Normalization:
    raw:
        - raw_fastqc
    qc:
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    Diginorm:
        - diginorm_load_into_counting
        - diginorm_normalize_by_median_pe_pe
        - diginorm_normalize_by_median_pe_se
        - diginorm_filter_abund
        - diginorm_extract_paired_reads
        - diginorm_merge_pe_single_reads
        - diginorm_get_former_se_reads
        - diginorm_fastqc
    """
    input:
        expand( # fastqc zip and html for raw PE data
            FASTQC_RAW_DIR + "/{sample}_{end}.fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            FASTQC_RAW_DIR + "/{sample}.fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "zip html".split()
        ),
        expand( # pe_pe
            NORM_DIR + "/{sample}.final.pe_pe.fq.gz",
            sample = SAMPLES_PE
        ),
        expand( # pe_se
            NORM_DIR + "/{sample}.final.pe_se.fq.gz",
            sample = SAMPLES_PE
        ) + expand( # se
            NORM_DIR + "/{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        ),
        expand( # fastqc zip and html files for PE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for PE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        )
        
        

###
# Assembly rules
###

rule assembly_split_pe_files:
    """
    Split pe_pe files into _1 and _2.
    """
    input:
        fastq_pe = NORM_DIR + "/{sample}.final.pe_pe.fq.gz"
    output:
        left  = temp(ASSEMBLY_DIR + "/{sample}_1.fq.gz"),
        right = temp(ASSEMBLY_DIR + "/{sample}_2.fq.gz")
    threads:
        1
    params:
        left  = "{sample}.final.pe_pe.fq.gz.1",
        right = "{sample}.final.pe_pe.fq.gz.2"
    log:
        "logs/assembly/split_pe_files_{sample}.log"
    benchmark:
        "benchmarks/assembly/split_pe_files_{sample}.json"
    shell:
        """
        ./bin/split-paired-reads.py \
            {input.fastq_pe} \
        > {log} 2>&1
        
        {gzip} -9c {params.left}  > {output.left}
        {gzip} -9c {params.right} > {output.right}
        
        rm {params.left} {params.right}
        """



rule assembly_merge_right_and_left:
    """
    Generate the left.fq and right.fq
    left.fq = /1 reads from PE + all SE reads
    right.fq = /2 reads form PE
    """
    input:
        forward=  expand(ASSEMBLY_DIR + "/{sample}_1.fq.gz",
            sample = SAMPLES_PE),
        reverse=  expand(ASSEMBLY_DIR + "/{sample}_2.fq.gz",
            sample = SAMPLES_PE),
        unpaired= expand( # pe_se
            NORM_DIR + "/{sample}.final.pe_se.fq.gz",
            sample = SAMPLES_PE
        ) + expand( # se
            NORM_DIR + "/{sample}.final.se.fq.gz",
            sample = SAMPLES_SE
        )
    output:
        left=  temp(ASSEMBLY_DIR + "/left.fq"),
        right= temp(ASSEMBLY_DIR + "/right.fq")
    threads:
        1
    log:
        "logs/assembly/merge_right_and_left.log"
    benchmark:
        "benchmarks/assembly/merge_right_and_left.json"
    shell:
        """
        {gzip} -dc {input.forward} {input.unpaired} > {output.left} 2> {log}
        {gzip} -dc {input.reverse} > {output.right} 2>> {log}
        """



rule assembly_run_trinity:
    """
    Assembly reads with Trinity.
    Notes on hardcoded settings:
        - Runs on paired end mode
        - Expect fastq files as inputs (left and right)
        - Does the full cleanup so it only remains a fasta file.
    """
    input:
        left  = ASSEMBLY_DIR + "/left.fq",
        right = ASSEMBLY_DIR + "/right.fq"
    output:
        fasta = protected(ASSEMBLY_DIR + "/Trinity.fasta")
    threads:
        BLOCK_THREADS
    params:
        memory= config["trinity_params"]["memory"],
        outdir= ASSEMBLY_DIR + "/trinity_out_dir"
    log:
        "logs/assembly/run_trinity.log"
    benchmark:
        "benchmarks/assembly/run_trinity.log"
    shell:
        """
        ./bin/Trinity \
            --seqType fq \
            --max_memory {params.memory} \
            --left {input.left} \
            --right {input.right} \
            --CPU {threads} \
            --full_cleanup \
            --output {params.outdir} \
        > {log}
        
        mv {params.outdir}.Trinity.fasta {output.fasta}
        
        rm {input.left}.readcount {input.right}.readcount
        """



rule assembly:
    """
    Runs from raw data to assembly:
    QC:
        - qc_trimmomatic_pe
        - QC_interleave_pe
        
    Diginorm:
        - diginorm_load_into_counting
        - diginorm_normalize_by_median_pe
        - diginorm_normalize_by_median_se
        - diginorm_filter_abund
        - diginorm_extract_paired_reads
        - diginorm_merge_single_reads
    Assembly:
        - assembly_prepare_reads
        - assembly_run_trinity
    """
    input:
        ASSEMBLY_DIR + "/Trinity.fasta",
        expand( # fastqc zip and html for raw PE data
            FASTQC_RAW_DIR + "/{sample}_{end}.fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "1 2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            FASTQC_RAW_DIR + "/{sample}.fastqc.{extension}",
            sample = SAMPLES_SE,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for PE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data DIGINORM
            FASTQC_NORM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for PE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_PE,
            pair   = PAIRS,
            extension = "zip html".split()
        ),
        expand( # fastqc zip and html files for SE data QC
            FASTQC_TRIM_DIR + "/{sample}.final.{pair}_fastqc.{extension}",
            sample = SAMPLES_SE,
            pair = ["se"],
            extension = "zip html".split()
        )



rule part_do_partition:
    """
    Construct a DBG from the reads and mark in the header of the fasta the
    component.
    """
    input:
        assembly = ASSEMBLY_DIR + "/Trinity.fasta"
    output:
        assembly = temp(PART_DIR + "/Trinity.fasta.part")
    params:
        assembly_tmp=  "Trinity.fasta.part",
        prefix=        config["prefix"],
        ksize=         config["diginorm_params"]["ksize"],
        hashsize=      config["diginorm_params"]["hashsize"],
        n_hashes=      config["diginorm_params"]["n_hashes"]
    threads:
        BLOCK_THREADS
    log:
        "logs/part/do_partition.log"
    benchmark:
        "benchmarks/part/do_partition.json"
    shell:
        """
        ./bin/do-partition.py \
            --ksize {params.ksize} \
            --n_hashes {params.n_hashes} \
            --hashsize {params.hashsize} \
            --threads {threads} \
            {params.prefix} \
            {input.assembly} \
        > {log} 2>&1
        
        mv {params.assembly_tmp} {output.assembly}
        mv {params.prefix}.info {PART_DIR}
        """



rule part_rename_with_partitions:
    """
    Rename de fasta assemblty according to the parition name and component.
    """
    input:
        assembly = PART_DIR + "/Trinity.fasta.part"
    output:
        assembly = protected(PART_DIR + "/Trinity.fasta")
    params:
        prefix       = config["prefix"],
        assembly_tmp = "Trinity.fasta.part.renamed.fasta.gz"
    threads:
        1
    log:
        "logs/part/rename_with_partitions.log"
    benchmark:
        "benchmarks/part/rename_with_partitions.json"
    shell:
        """
        python src/eel-pond/rename-with-partitions.py \
            {params.prefix} \
            {input.assembly} \
        > {log}
        
        {gzip} -dc \
            {params.assembly_tmp} \
        >  {output.assembly} \
        2>> {log}
        
        rm {params.assembly_tmp}
        """
