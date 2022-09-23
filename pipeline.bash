## Execution script for use on computing cluster

# Adjust to desired output and available cores prior to execution. 
# Syntax is {sample}_{reference}_{bin_size}
snakemake --cores 8 results/cnvs/sampleID_hg38_20k.pdf
