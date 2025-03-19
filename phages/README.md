# Phages annotation pipeline

# Purpose

The goal of this pipeline is to assemble phage genomes, annotate genes encoded by the genomes, and assess genome quality and taxonomic classification of the phages. 

# Workflow steps

1. Trim sequences using Trimmomatic
2. Genome assembly: remove host contaminants from reads, and used leftover (unmapped) reads for assembly using SPADES.
3. Perform genome annotation using Bakta.
4. Find closely related genomes on NCBI using BLASTn. 
5. Assess phage genome quality (completeness) using CheckV.
6. Assign taxonomy to phages using TaxMyPhage.

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

## References
This pipeline uses the following tools:

- Trimmomatic v0.39 : https://github.com/timflutre/trimmomatic
- SPADES v4.0.0: https://github.com/ablab/spades
- Blastn v.2.16.0: https://blast.ncbi.nlm.nih.gov/Blast.cgi
- CheckV v.1.0.3: https://bitbucket.org/berkeleylab/checkv/src/master/
- TaxMyPhage v.0.3.4: https://github.com/amillard/tax_myPHAGE/blob/main/README.md





