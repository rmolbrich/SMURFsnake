### SMURFSeq protocol in ingularity using snakemake ###
1. 'config.yaml' must be updated with filenames for samples 
2. for alignment with 'bwa', the path to the reference must be updated
2.1. alignment works best with hard masked reference
2.2. at the moment bed-files and binning scheme is build for hg19 only!
3. adjust 'SLURM' parameters in cluster-script
4. adjust '$WORKDIR' variable in cluster-script
5. 'pipeline.bash' contains the Snakemake execution command
5.1. adjust execution command in this file
5.2. or use the commented code-block in cluster-script to dynamically create exec-cmd

# Execute with
singularity exec --cleanenv SMURFsnake.sif bash pipeline.bash
# Code in pipeline.bash
snakemake --cores <INT> results/cnvs/A.pdf


# Also works as shell
singularity shell --cleanenv SMURFsnake.sif
