#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("msa")

library(msa)

# load phage genome:
phage <- system.file("~/Downloads/Alignment_MutationinEvo2.msa")
phage_seq<- readDNAStringSet("~/Downloads/Alignment_MutationinEvo2.msa",
                                   format="fasta")
phage_seq

#if (!requireNamespace("devtools", quietly=TRUE))
#  install.packages("devtools")
#devtools::install_github("YuLab-SMU/ggmsa")

# Plot it:
library(ggmsa)
msa_plot <- ggmsa(phage_seq, 
      start = 200, 
      end = 240, 
      char_width = 0.5,
      seq_name = TRUE,
      color = "Shapely_NT") + geom_seqlogo() + geom_msaBar()

msa_plot


class(msa_plot)
library(ggplotify)
as.ggplot(msa_plot)

library(tidyverse)

ggsave(plot = as.ggplot(msa_plot),
       filename = "13Aug2024_PhageMutation/mutation.pdf", 
       width=7,
       height=1.2)
