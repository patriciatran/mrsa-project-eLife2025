# RNAseq

# Purpose

The purpose of the RNA-seq pipeline is to map short-reads (Illumina) transcripts onto reference *Staph aureus* genomes, count transcripts across the different treatments. The data is processed using DESeq2 to identify differentially expressed genes between treatments.

# Cyberinfrastructure & Implementation

The RNA sequence processing was done using `HTCondor` on CHTC using submit and executable files available under `./scripts`. Apptainer definition files are available under `./recipes`. After the data was processed in the form of a gene count matrix, the data was exported and analysed using R and Rstudio using the code under `./Rcode`. 

# Workflow steps

## CHTC
1. Mapping reads to reference genomes 
2. Count transcripts mapping to gene features
3. Use mmseqs2 to find protein similarity between the two reference genomes.

## R
1. Perform DEseq2 analysis on both genomes
2. Combine and merge results to compare up-regulated and down-regulated genomes between the two genomes.


# Repository files

This repository contains 2 folders: `recipes` and `scripts`.
The `recipes` folder contains the Apptainer definition files needed to create the Apptainer sif files. 
The `scripts` folder contains the HTcondor `.sh` and `.sub` files.

# Input folder set-up

# Expected output directory structure

##  Building containers

To build the software containers, you will need to start an interactive job, build the container, test it, and move it to a location accessible by the working nodes (e.g. staging, not home).
For detailed instructions, visit https://github.com/UW-Madison-Bacteriology-Bioinformatics/chtc-containers. 

brief instructions:
```
cd recipes
nano build.sub
# change the file listed in the transfer_input_files line
condor_submit -i build.sub
# replace "container" with the name of your choice
apptainer build container.sif container.def
apptainer shell -e container.sif
# test container by typing the -h --help command.
exit
mv container.sif /staging/netid/apptainer/.
exit
cd ..
```

## Run code

## Next steps
This workflow will create transcripts count files. I recommend using Globus.org to transfer files to your ResearchDrive or to your personal endpoint.
For instructions, please visit: https://chtc.cs.wisc.edu/uw-research-computing/globus

## Importing to R / Rstudio for analysis with DESeq2

An Rscript is available at [./Rcode](./RCode) showing how the data was processed for the paper.

## References
This pipeline uses the following tools:

- Bowtie2 v.2.5.4: https://github.com/BenLangmead/bowtie2
- FeatureCount v.2.0.8 : https://subread.sourceforge.net/featureCounts.html
- R v.4.4.0 : https://cran.r-project.org/bin/macosx/base/
- DESeq v.1.44.0: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
- Tidyverse 2.0.0: https://www.tidyverse.org/
- EnhancedVolcanoPlots (version 1.22.0): https://github.com/kevinblighe/EnhancedVolcano






