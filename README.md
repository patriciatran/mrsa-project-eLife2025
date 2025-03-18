# MRSA Project

This is the code for the project https://elifesciences.org/reviewed-preprints/102743v1.
[include project summary here]

This code repository contains the container recipes and scripts using in a few of the major bioinformatics analyses in the paper:
1. Reference MRSA genomes variant identification (`./variants`)
2. Characterization of phages (`./phages`)
3. Differential gene expression analysis (`./rnaseq`)

This main README contains a summary of the repository, but you go to each folder (`./variants`, `./phages`, and `./rnaseq`) for more detailed information such as software versions, description of code, etc.

# Description

## Variant annotation

The goal of this pipeline is to perform variant calling on long-read (e.g. Oxford Nanopore ONT technologies) sequencing data against a reference bacterial genome of interests. The program generated multiple alignment files in the .SAM and .BAM format. Additionally, it includes two steps to functionally annotate the bacterial genomes using multiple reference genome databases. This can help match in which genes the variants tends to occur.

## RNA-seq

The purpose of the RNA-seq pipeline is to map short-reads transcript onto reference genomes, count transcripts across the different treatments. The data is processed using DESeq2 to identify differentially expressed genes between treatments.

## Phage annotations
To add?

# Cyberinfrastructure

# Citation

Bacteriophage infection drives loss of β-lactam resistance in methicillin-resistant Staphylococcus aureus
My Tran, Angel J Hernandez Viera, Patricia Q Tran, Charlie Y Mo
Department of Bacteriology, University of Wisconsin-Madison, USA
https://doi.org/10.7554/eLife.102743.1
