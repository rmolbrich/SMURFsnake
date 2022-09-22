# Fail fast
set -o errexit

# Fill this file with your environment creation
mamba env create -q --name snakemake --file /opt/environment.yaml

# Not sure about this.. works so keep it
conda activate snakemake