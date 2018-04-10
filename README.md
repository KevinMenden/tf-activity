# TF-activity workflow
Pipeline to infer active transcription factors from CAGE-seq data using transcription factor binding motifs.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker Repository on Dockerhub](https://img.shields.io/badge/docker-available-green.svg "Docker Repository on Dockerhub")](https://hub.docker.com/r/kevinmenden/tf-activity/)
![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)


## Description
Given a set of CAGE-peaks of interest (e.g. up-regulated in some condition), this pipeline will extract the genomic sequence
in a given range around these peaks and look for enrichment of TF motifs in these sequences. As background either the
shuffled input sequences will be used or a background extracted from user supplied CAGE-peaks that are not of interest.

The pipeline will annotate the peaks and create TF-gene mappings. Additionally, the CAGE-peaks are overlapped
with ChIP-seq peaks extracted from the [ENCODE](https://www.encodeproject.org/) project.

## Usage
As minimum input for this pipeline you will need:
* CAGE peaks of interest in bed format
* fasta of reference genome
* gtf of reference genome

#### `--fasta`
Used to specify the path to the reference genome.

#### `--gtf`
Used to specify path to the GTF file of the reference genome.

#### `--peaks`
Your cage peaks of interest in BED file format.

#### `--background`
Peaks that are not differentially expressed can be used here as background peaks. They have to be in BED format
just like your peaks of interest.

#### `--pfms`
A file containing all the TF motifs to use. This file must be in homer format. If you do not have a file in homer format,
you can instead specifiy a TF motif file in Jaspar format using the `--pfms_jaspar` flag. If none of these two flags
is set, the pipeline will use all motifs from the Jaspar core collection.

#### `-profile`
Which profile to use. Use `docker` to use the docker container provided.




Additional input options:


Example pipeline call:
```bash
nextflow run kevinmenden/tf-activity -profile docker --fasta path/to/genome.fa --gtf path/to/gtf/genome.gtf \
--peaks peaks.bed --background background.bed
```

In the above example, the pipeline will use the [docker image](https://hub.docker.com/r/kevinmenden/tf-activity/) from
dockerhub to run. Because no motifs are specified, the default Jaspar core motifs will be used.
