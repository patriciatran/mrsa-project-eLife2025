![Static Badge](https://img.shields.io/badge/Manuscript-Accepted-159e2b)
![Static Badge](https://img.shields.io/badge/Code-Currently_Verifying-fcd253)

# MRSA Project

This is the code for the project https://elifesciences.org/reviewed-preprints/102743v1.

This code repository contains the container recipes and scripts using in a few of the major bioinformatics analyses in the paper:
1. Reference MRSA genomes variant identification ([./variants](./variants))
2. Characterization of phages ([./phages](./phages))
3. Differential gene expression analysis ([./rnaseq](./rnaseq))

This main README contains a summary of the repository, individual readmes for each are also available ([./variants/README.md](./variants/README.md), [./phages/README.me](./phages/README.md), and [./rnaseq/README.md](./rnaseq/README.md)) for more detailed information such as software versions, description of code, etc.

# Description

## Genomic variant annotation

The goal of this pipeline is to perform variant calling on long-read (e.g. Oxford Nanopore ONT technologies) sequencing data against a reference bacterial genome of interests. The program generated multiple alignment files in the .SAM and .BAM format. Additionally, it includes two steps to functionally annotate the bacterial genomes using multiple reference genome databases. This can help match in which genes the variants tends to occur.

## RNA-seq

The purpose of the RNA-seq pipeline is to map short-reads (Illumina) transcripts onto reference *Staph aureus* genomes, count transcripts across the different treatments. The data is processed using DESeq2 to identify differentially expressed genes between treatments.

## Phage annotations

The goal of this pipeline is to assemble phage genomes, annotate genes encoded by the genomes, and assess genome quality and taxonomic classification of the phages. 

# Cyberinfrastructure

Each folder contains `HTCondor` submit and executable files under `folder/scripts`. Additionally definition files to create `Apptainer` container images are available under `folder/recipes`.

# Citation

Bacteriophage infection drives loss of β-lactam resistance in methicillin-resistant Staphylococcus aureus
My Tran, Angel J. Hernandez Viera, Patricia Q. Tran, Erick Nilsen, Lily Tran, Charlie Y. Mo
[https://doi.org/10.7554/eLife.102743.1](https://doi.org/10.7554/eLife.102743.1)

The code in this repository was created and implemented by Patricia Q. Tran

