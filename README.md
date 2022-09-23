# SMURFsnake

`SMURFSeq` pipeline for CNV-analysis of Nanopore as `Snakemake` pipeline within a `Singularity` container


## Set-up

* Update `config.yaml` with correct paths for reference genomes. Use expressive tags, i.e., `hg19` and `hg38`, to designate the choices. Tags may be added or removed, keep in mind that they will become part of the output filenames.

```
ref: 
    hg19: path/to/reference/hg19.fa
    hg38: path/to/reference/hg38.fa.masked
```

* The remaining paths in the `config.yaml` point to directories within the Singularity container. If you choose to run the pipeline without the image, the required scripts and binning-files are located in the `SMURFSeq.zip` within the `rebuild` directory. Alternatively, the files may be obtained from the original source at [smithlabcode](https://github.com/smithlabcode/smurfseq_scripts.git).

* Construction of the singularity image requires a Unix machine with sudo rights and the singularity library in version >= 1.6. Additional information is provided by the ReadMe located in the rebuild directory.

* The file `pipeline.bash` contains the Snakemake execution command


## Run as job on a cluster

Put the following line in the cluster-script to run pipeline as intended.

```
singularity exec --cleanenv SMURFsnake.sif bash pipeline.bash

# Code in pipeline.bash - Syntax is {sample}_{reference}_{bin_size}
snakemake --cores <INT> results/cnvs/{sample}_{reference}_{bin_size}.pdf
```

## Run in shell

```
singularity shell --cleanenv SMURFsnake.sif
```


