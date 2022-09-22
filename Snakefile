### Implementation of SMURFSeq protocol from Smith et. al 2019 to work with Snakemake
## by Michael Olbrich 2022


## Project configuration ##
configfile: "config.yaml"

## Use tags 'hg19' or 'hg38' to select for reference genome as listed in 'config.yaml'
# 
# This is used in 'rule bwa_map'
def get_ref_genome(wildcards):
	return config["ref"][wildcards.reference]

## Use tags 'hg19' or 'hg38' to select chrom sizes for reference
# 
# This is used in 'rule count_bins'
def get_ref_chrom(wildcards):
	return config["chrs"][wildcards.reference]

## Use tags for genome-reference and bin-size ('5k','20k','50k') to find correct files
#
# This is used in 'rule count_bins'
def get_bin_size(wildcards):
	return config["bin"][wildcards.reference][wildcards.bin_size]

## Use tags for genome-reference and bin-size ('5k','20k','50k') to find correct files
#
# This is used in 'rule analyse_cnv'
def get_bin_ex(wildcards):
	return config["bex"][wildcards.reference][wildcards.bin_size]

## Use tags for genome-reference and bin-size ('5k','20k','50k') to find correct files
#
# This is used in 'rule analyse_cnv'
def get_bin_gc(wildcards):
	return config["bgc"][wildcards.reference][wildcards.bin_size]

## Project rules ##

# Use 'bwa' for alignment of SMURFed reads
## ToDo: use minimap2
rule bwa_map:
	input:
		"data/raw_reads/{sample}.fastq"
	output:
		"data/mapped_reads/{sample}_{reference}.sam"
	params:
		ref=get_ref_genome
	log:
		"log/bwa_mem/{sample}_{reference}.log"
	threads: 8
	shell:
		"bwa mem -x ont2d -k 12 -W 12 -A 4 -B 10 -O 6 -E 3 -T 120 -t {threads} {params.ref} {input} > {output}"

# Use samtools to filter ambiguously mapped reads
rule filt_ambig:
	input: 
		"data/mapped_reads/{sample}_{reference}.sam"
	output:
		"data/unambig_reads/{sample}_{reference}.sam"
	shell:
		"samtools view -h -q 1 -e '[AS]>=120' {input} > {output}" 

# Calculate bin-counts
rule count_bins:
	input:
		"data/unambig_reads/{sample}_{reference}.sam"
	output:
		bcbed="data/bin_counts/{sample}_{reference}_{bin_size}_bin_counts.bed",
		bstxt="data/bin_counts/{sample}_{reference}_{bin_size}_bin_stats.txt"
	params:
		script=config["scripts"]["count"],
		chrs=get_ref_chrom,
		bins=get_bin_size
	shell:
		"{params.script} -i {input} -c {params.chrs} -b {params.bins} -o {output.bcbed} -s {output.bstxt}"

# Analyse CNVs
rule analyse_cnv:
	input:
		"data/bin_counts/{sample}_{reference}_{bin_size}_bin_counts.bed"
	output:
		pdf="data/cnvs/{sample}_{reference}_{bin_size}.pdf",
		dat="data/cnvs/{sample}_{reference}_{bin_size}.data.txt",
		sho="data/cnvs/{sample}_{reference}_{bin_size}.short.txt"
	params:
		script=config["scripts"]["analysis"],
		bingc=get_bin_gc,
		binex=get_bin_ex
	shell:
		"{params.script} {input} data/cnvs/{wildcards.sample}_{wildcards.reference}_{wildcards.bin_size} {params.bingc} {params.binex}"























