# Singularity Image for SMURFSeq - Snakemake pipeline

## Why do I want this:

The SMURFsnake pipeline can just be used as is, given that the 'config.yaml' was adjusted properly. However, the requirements comprise tools (bwa, samtools), frameworks (python, R), and software libraries (for R and python) that need to be available. A singularity image provides all these requirements.

## Why isn't the image a part of the repository:

The finished image-file is roughly 800MB in size and it is more convenient to provide the recipe for construction than to put it in the repository.

## What are these files:

### SMURFsnake.def

The singularity recipe that is required for the build-process.

### environment.yaml

Specifies all the libraries and tools that are part of the created conda environment.

### R_packages.R

Specifies the packages that are installed for the R-framework.

### post.bash

A workaround that allows to set-up and enable conda environments during image building.

### SMURFSeq.zip

This archive contains the SMURFSeq computational scripts and bin files for cnv-analyses. These are stored within the singularity image.