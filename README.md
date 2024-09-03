# vmp_sqb: a viral metagenomics pipeline for squidbase

This repository contains a Nextflow pipeline for processing FASTQ files, running Kraken2 for taxonomic classification, extracting relevant reads, and performing mapping and coverage analysis. The pipeline is designed for multiplexed nanopore runs, specifically tailored for viral detection.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Processes](#about)


## Installation

### Prerequisites
- [Nextflow](https://www.nextflow.io/)
- [Conda](https://docs.conda.io/en/latest/)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)
- [seqtk](https://github.com/lh3/seqtk)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://www.htslib.org/)
- [pandas](https://pandas.pydata.org/)

Note that this software does not need to be installed when running the pipeline with the default conda option; an `vmp_env.yaml` has been provided to run the pipeline.

### Clone the Repository
```bash
git clone https://github.com/Cuypers-Wim/vmp_sqb.git
cd vmp_sqb
```

### Configure the Pipeline
Edit the nextflow.config file to point to your data directories and set any necessary parameters.

## Usage 
```bash
nextflow run main.nf 
```

Running this pipeline generates a folder named `results`, which contains the following subdirectories:

- `mapping`: Contains mapping statistics and BAM files related to the read alignment and filtering processes.
- `reads`: Contains the extracted reads in FASTQ format that correspond to the target species.
- `references`: Stores the reference genomes used for mapping.
- `POD5`: Includes POD5 files containing raw nanopore signal data ("squiggles") specific to the identified species.

## About

This pipeline is designed to process data generated from a nanopore sequencing run. Typically, such a run produces folders containing FASTQ and POD5 files, where each barcode folder may contain multiple subfolders with these files. The pipeline is particularly useful when you need to identify and analyze a specific species within a sample, such as the Chikungunya virus with NCBI taxon identifier 37124, expected to be found in barcode one.

To run the pipeline, you need to provide a CSV file that specifies the expected species and their corresponding NCBI taxon identifiers for each barcode. The pipeline uses this information to extract reads corresponding to the target species from each barcode folder. It then performs mapping to confirm the presence of the species and eliminate any false positives.

The pipeline outputs mapping statistics for the reference genome or the best matching reference genome. Additionally, if required, the pipeline can generate POD5 files containing the raw nanopore signal data (the so-called "squiggles") specific to the identified species. These POD5 files can be further uploaded to platforms like Squidbase for detailed analysis.

###  Example input CSV file

| Barcode   | Virus Name | Taxon ID |
|-----------|------------|----------|
| barcode01 | CHIKV-1    | 37124    |
| barcode02 | Denv1-5    | 11053    |
| barcode03 | Denv2-7    | 11060    |
| barcode04 | Denv3-10   | 11069    |
| barcode05 | EEEV-16    | 11021    |


### Descriptions of Key Processes:

- `CONCATENATE_FASTQ`: Merges all FASTQ files in a barcode folder into a single FASTQ file for downstream analysis.
- `RUNKRAKEN2`: Executes Kraken2 on the concatenated FASTQ files for taxonomic classification, generating reports and extracting read IDs.
- `SUBSET_FASTQ`: Extracts reads of interest from the concatenated FASTQ files based on the Kraken2 classification results.
- `RETRIEVE_GENOMEREFS`: Downloads and formats reference genomes based on specified Taxon IDs.
- `MAP_READS`: Aligns the subsetted FASTQ reads to the reference genomes using minimap2 (to remove false-positives), followed by sorting the output BAM file.
- `BAM_PROCESSING`: Marks duplicates in the BAM file and creates an index for efficient access.
- `MAPPING_STATS_ALL`: Generates mapping statistics for the BAM files, including flagstat and idxstats reports.
- `EXTRACT_READIDS`: Filters the BAM file to retain reads that map to the target species and extracts the corresponding read IDs.
- `SAMTOOLS_DEPTH`: Computes the depth of coverage across the reference genomes using the filtered BAM files.
- `DEPTH_OF_COVERAGE`: Analyzes the depth of coverage data to calculate coverage statistics, outputting the results in CSV format.
- `SUBSET_POD5`: Filters POD5 files to retain only the raw signals/squiggles corresponding to the extracted read IDs.

Each process is designed to handle specific aspects of the analysis workflow, making the pipeline flexible and easy to extend.

