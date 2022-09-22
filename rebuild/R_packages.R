#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DNAcopy", update=FALSE, ask=FALSE)