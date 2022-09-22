### Implementation of SMURFSeq protocol from Smith et. al 2019 to work with Snakemake
## by Michael Olbrich 2022


## Project configuration ##

# Project config to set-up variables and input files.
configfile: "config.yaml"

# Define final output files to facilitate execution of pipeline with single command.
# rule all:
# 	input:
# 		"results/cnvs/{sample}.pdf"

## ToDo: understand why exactly this approach is required
# Retrieve input files from config file.
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]


## Project rules ##

# Use 'bwa' for alignment of SMURFed reads
## ToDo: use minimap2
rule bwa_map:
	input:
		get_bwa_map_input_fastqs
	output:
		"results/mapped_reads/{sample}.sam"
	params:
		ref=config["ref"]
	log:
		"log/bwa_mem/{sample}.log"
	threads: 8
	shell:
		"bwa mem -x ont2d -k 12 -W 12 -A 4 -B 10 -O 6 -E 3 -T 120 -t {threads} {params.ref} {input} > {output}"

# Use samtools to filter ambiguously mapped reads
rule filt_ambig:
	input: 
		"results/mapped_reads/{sample}.sam"
	output:
		"results/unambig_reads/{sample}.sam"
	shell:
		"samtools view -h -q 1 -e '[AS]>=120' {input} > {output}" 

# Calculate bin-counts
rule count_bins:
	input:
		"results/unambig_reads/{sample}.sam"
	output:
		bcbed="results/bin_counts/{sample}_bin_counts.bed",
		bstxt="results/bin_counts/{sample}_bin_stats.txt"
	params:
		script=config["scripts"]["count"],
		chrs=config["chrs"],
		bins=config["bins"]["bin_size"]
	shell:
		"{params.script} -i {input} -c {params.chrs} -b {params.bins} -o {output.bcbed} -s {output.bstxt}"

# Analyse CNVs
rule analyse_cnv:
	input:
		"results/bin_counts/{sample}_bin_counts.bed"
	output:
		pdf="results/cnvs/{sample}.pdf",
		dat="results/cnvs/{sample}.data.txt",
		sho="results/cnvs/{sample}.short.txt"
	params:
		script=config["scripts"]["analysis"],
		bingc=config["bins"]["bin_gc"],
		binex=config["bins"]["bin_ex"]
	shell:
		"{params.script} {input} results/cnvs/{wildcards.sample} {params.bingc} {params.binex}"























