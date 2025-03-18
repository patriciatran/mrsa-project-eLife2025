# Phages annotation pipeline

# Purpose

# Cyberinfrastructure & Implementation

These scripts are meant to be run by HTCondor, a workflow manager that takes in an executable file and a submit file. It is meant to be run on the UW-Madison shared campus-computing infrastructure CHTC, but could work on other systems with few modifications. The pipeline takes advantage of the high-throughput scaling abilities of HTCondor to submit multiple jobs at the same time. For example, if we had an experimental design of 3 treatment, 3 replicate and 3 reference genomes and 10 time points, we would have to perform the steps 3 * 3 * 3 * 10 = 270 times. Instead, we write the metadata (information about the sample design) in a comma separated file containing 270 rows, and HTCondor will submit 270 jobs at the same time for us.

# Workflow steps


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

An Rscript is available at XXXXXX showing how the data was processed for the paper.

## References
This pipeline uses the following tools:

- Software: link
- Software: link
- Software: link
- Software: link





