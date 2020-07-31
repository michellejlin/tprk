# tprK pipeline
This pipeline was designed to take Illumina and PacBio files straight off the sequencer to a final comparison table of all the different variable regions with their relative frequencies, as well as various pretty plots along the way.

## Table of Contents
* [Setup](#Setup)
* [Input Files](#Input-Files)
* [Usage](#Usage)
* [Common Errors](#Common-Errors)

## Setup

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/). The first time running this program will take a while, as the docker image will take some time to build, but this is a one time thing!


## Input Files
Put the following things in one folder:
- **All the sequence files to run analysis on**
    - PacBio Q20 reads, gzipped
    - Single-end Illumina reads trimmed and run through Trimmomatic, gzipped
    - By default the pipeline expects both PacBio and Illumina files for every sample. Running the pipeline with just PacBio or just Illumina files is possible with the `--pacbio` and `--illumina` flags respectively. However, some plots require both files to be generated and these plots will not be output.
- **Metadata file.** This should be a .csv with three columns: SampleName, PacBio, Illumina, shown in the table below. Make sure to include absolute paths to PacBio and Illumina files!
    - This file should be placed in the same folder as your files to be analyzed.
    - There MUST be a newline character at the end of this file to be read as a valid csv. Simply hit enter in the last row to ensure there is a valid new line.
    - Ensure that there are no special characters, including hyphens! Underscores are okay.
    - If running just Illumina or just PacBio, simply leave those columns blank (but make sure to have commas as appropriate).
    - Example metadata files (for both, just Illumina, and just PacBio) are provided in the `example/` folder. The general format of the metadata file should be three columns, separated by commas, as shown:

| SampleName  | Illumina  | PacBio |
| ------------- | ------------- | ------------- |
| This will largely be the name used for generating tables and plots. | Should be in format Ill_[sample name].fastq.gz. The Illumina file specified for the sample name. This must match exactly the name of the matching file in the folder. This should be a trimmed file run through Trimmomatic. | Should be in format PB_[sample name].fasta.gz. The PacBio file specified for the sample name. This must match exactly the name of the matching file in the folder. This should be a Q20 file.  | 

## Usage
- Example command for just Illumina files in current directory on a laptop without many CPUs: ```nextflow run michellejlin/tprk -r nextflow --INPUT ./ --OUTDIR output/ --ILLUMINA --METADATA metadata.csv -resume -with-docker ubuntu:18.04 -with-trace -profile laptop```
- Example command for comparing PacBio and Illumina files with specified cutoffs on the cloud with a large dataset: ```AWS_PROFILE=covid nextflow run michellejlin/tprk -r nextflow --INPUT example/ --OUTDIR example/output/ --METADATA metadata.csv --LARGE -resume -with-docker ubuntu:18.04 -with-trace -c ~/nextflow.covid.config -profile Cloud```

For a list of arguments, you can also run ```nextflow run michellejlin/tprk -r nextflow --help``` .

| Command  | Description |
| ---      | ---         | 
| --INPUT  | Input folder where gzipped fastqs are located. For current  directory, `./` can be used.
| --OUTDIR | Output folder where .bams and consensus fastas will be piped into.
| --METADATA | Path to metadata file with specific format. 
| --PACBIO | Specify that there are only PacBio files to be read.
| --ILLUMINA | Specify that there are only Illumina files to be read.
| --LARGE | Specify that this is a large dataset. Splitting of visualizations will be done.
|--RF_FILTER | Specify relative frequency filter. Default is 0.2.
|--COUNT_FILTER | Specify count filter. Default is 5.
|--ILLUMINA_FILTER | Specify whether PacBio reads should be filtered to only include files supported by Illumina reads that reach the cutoff.
| -resume  | nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-docker ubuntu:18.04 | Runs command with Ubuntu docker.
| -with-trace | Outputs a trace.txt that shows which processes end up in which work/ folders. 

## Common Errors
1. `incomplete final line found by readTableHeader on '/Users/uwvirongs/Documents/tprk/metadata.csv'`
   Make sure your metadata file has a new line at the end. You can do this by simply pressing enter on the last line of your file and saving. 
