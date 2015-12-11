#!/bin/bash
set -euo pipefail

# Environment


## Folders
raw_dir=data/fastq_raw
adapt_dir=data/adapters
adapt_log=logs/adapters

out_file=${adapt_dir}/adaptors.tsv
log_file=${adapt_log}/adaptors.log


# Create folders
mkdir -p $adapt_dir $adapt_log

# fastq_to_sequence: from a piped fastq, return only the sequence field (those in the 4n + 2 lines)
# Input: a piped fastq file. Pipe it from gzip or whatever compressor: `pigz -dc file.fastq.gz | fastq_to_sequence`
# Output: only the nucleotide sequences, through stdout
# Note: For better performance, use parallel, with --keep-order optional
fastq_to_sequence(){
    awk ' NR % 4 == 2 { print } '
}


# count_adapters: for each sequence, find if it contains any adapter
# Input: raw sequences, one per line, through stdin
# Output: the number of adapters detected
# Note: combine with parallel for faster results
count_adapters(){
    awk '                                   \
    BEGIN { adapters[ "TruSeq2-SE"   ] = 0;   \
            adapters[ "TruSeq2-PE"   ] = 0;   \
            adapters[ "TruSeq3-SE"   ] = 0;   \
            adapters[ "TruSeq3-PE"   ] = 0;   \
            adapters[ "TruSeq3-PE-2" ] = 0;   \
            adapters[ "NexteraPE"    ] = 0; } \
    /AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT/       { adapters[ "TruSeq2-PE"   ] += 1 };  \
    /AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC/                               { adapters[ "TruSeq3-PE-2" ] += 1 };  \
    /AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG/                                { adapters[ "TruSeq2-SE"   ] += 1 };  \
    /AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG/    { adapters[ "TruSeq2-PE"   ] += 1 };  \
    /AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT/                                { adapters[ "TruSeq2-SE"   ] += 1 };  \
    /AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA/                               { adapters[ "TruSeq3-PE-2" ] += 1 };  \
    /AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT/       { adapters[ "TruSeq2-PE"   ] += 1 };  \
    /AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG/                               { adapters[ "TruSeq2-SE"   ] += 1 };  \
    /AGATGTGTATAAGAGACAG/                                              { adapters[ "NexteraPE"    ] += 1 };  \
    /CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT/    { adapters[ "TruSeq2-PE"   ] += 1 };  \
    /CTGTCTCTTATACACATCTCCGAGCCCACGAGAC/                               { adapters[ "NexteraPE"    ] += 1 };  \
    /CTGTCTCTTATACACATCTGACGCTGCCGACGA/                                { adapters[ "NexteraPE"    ] += 1 };  \
    /GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG/                               { adapters[ "NexteraPE"    ] += 1 };  \
    /GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT/                               { adapters[ "TruSeq3-PE"   ] += 1 ; adapters[ "TruSeq3-PE-2" ] += 1 };  \
    /TACACTCTTTCCCTACACGACGCTCTTCCGATCT/                               { adapters[ "TruSeq3-PE"   ] += 1 ; adapters[ "TruSeq3-PE-2" ] += 1 };  \
    /TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG/                                { adapters[ "NexteraPE"    ] += 1 };  \
    /TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC/                          { adapters[ "TruSeq2-PE"   ] += 1 };  \
    /TTTTTTTTTTCAAGCAGAAGACGGCATACGA/                                  { adapters[ "TruSeq2-PE"   ] += 1 };  \
    END {   print   adapters["TruSeq2-SE"] "\t" adapters["TruSeq2-PE"] "\t" adapters["TruSeq3-SE"] "\t" \
                    adapters["TruSeq3-PE"] "\t" adapters["TruSeq3-PE-2"] "\t" adapters["NexteraPE"] } '
}

# collect_results: when done in parallel count_adapters, sums the partial results
# Input: the execution of count_adapters done in parallel. If not, this is redundant
# Output: a table with the

collect_results(){
    awk \
        'BEGIN {    adapters[ "TruSeq2-SE"   ] = 0;   \
                    adapters[ "TruSeq2-PE"   ] = 0;   \
                    adapters[ "TruSeq3-SE"   ] = 0;   \
                    adapters[ "TruSeq3-PE"   ] = 0;   \
                    adapters[ "TruSeq3-PE-2" ] = 0;   \
                    adapters[ "NexteraPE"    ] = 0; } \
        {   adapters[ "TruSeq2-SE"   ] += $1 ;      \
            adapters[ "TruSeq2-PE"   ] += $2 ;      \
            adapters[ "TruSeq3-SE"   ] += $3 ;      \
            adapters[ "TruSeq3-PE"   ] += $4 ;      \
            adapters[ "TruSeq3-PE-2" ] += $5 ;      \
            adapters[ "NexteraPE"    ] += $6 ;  }   \
        END {   print   adapters["TruSeq2-SE"] "\t" \
                        adapters["TruSeq2-PE"] "\t" \
                        adapters["TruSeq3-SE"] "\t" \
                        adapters["TruSeq3-PE"] "\t" \
                        adapters["TruSeq3-PE-2"] "\t" \
                        adapters["NexteraPE"]     } '
}

find_adapters(){

    fastq_to_sequence |
    count_adapters    |
    collect_results

}

export -f fastq_to_sequence
export -f count_adapters
export -f collect_results
export -f find_adapters

echo -e \
    File"\t"TruSeq2-SE"\t"TruSeq2-PE"\t"TruSeq3-SE"\t"TruSeq3-PE"\t"TruSeq3-PE-2"\t"NexteraPE \
> $out_file

parallel --keep-order --tag \
    gzip -dc {} "|" find_adapters ::: ${raw_dir}/*.fq.gz \
>>  $out_file \
2>  $log_file

