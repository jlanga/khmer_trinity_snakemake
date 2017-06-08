# khmer_trinity_snakemake

## NOTE

This is an old repo. For a more efficient and automatic approach, go to [smsk_khmer_trinity](https://github.com/jlanga/smsk_khmer_trinity/)


## What's this?

This is a repo containint a Snakemake workflow to assemble RNA-Seq reads.

The general workflow is as follows:

- Read trimming with `Trimmomatic`

- Read normalization with the `khmer`/`diginorm` Python package.

- Transcriptome assembly with Trinity.

## How to run it?

1. Clone this repo

2. Make a virtualenv with python>=3.4:

```
virtualenv --python=python3.5 .
```

3. Activate it

```
source bin/activate
```

4. Install locally the required software:

```
bash scripts/install_software.sh
```

5. Test the pipeline with the sample data to see if it works (Trimmomatic, Trinity, khmer and additional pip packages):

```
ln -s config.yaml.example config.yaml
snakemake -j 24
```

## How to analyse your data?

Copy the config.yaml.example into config.yaml and fill the file with your data:

- Paths to your FASTQ files

- Adaptors used

- Phred Scores

- Additional parameters for Trimmomatics, Khmer and Trinity


Have fun!
