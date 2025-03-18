# Phages annotation pipeline

# Purpose

# Cyberinfrastructure & Implementation


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





