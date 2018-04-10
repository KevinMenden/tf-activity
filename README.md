# TF-activity workflow
Pipeline to infer active transcription factors from CAGE-seq data using transcription factor binding motifs.


## Description
Given a set of CAGE-peaks of interest (e.g. up-regulated in some condition), this pipeline will extract the genomic sequence
in a given range around these peaks and look for enrichment of TF motifs in these sequences. As background either the
shuffled input sequences will be used or a background extracted from user supplied CAGE-peaks that are not of interest.

The pipeline will annotate the peaks and create TF-gene mappings. Additionally, the CAGE-peaks are overlapped
with ChIP-seq peaks extracted from the [ENCODE] project.

## Usage
As minimum input for this pipeline you will need:
* CAGE peaks of interest in bed format
* fasta of reference genome
* gtf of reference genome


[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](https://www.nextflow.io/)
![Docker Repository on Dockerhub](https://img.shields.io/badge/docker-available-green.svg "Docker Repository on Dockerhub")
![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)

