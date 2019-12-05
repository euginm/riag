#!/bin/bash

# general parameters

## temp directory, optional, use if dont want to use temp dir from PATH
#tmp_dir=/path/to/tmp/dir

## output directory
out_dir=/path/to/out/dir

## output prefix, use "rdna_pipeline" if not set
#prefix='my_prefix'

## path to the folder with long reads fasta/q files, files can be compressed
lr_dir=/path/to/long/reads

## Long reads approximate coverage (need for rDNA ropy number estimation), will be computed if not set
## To avoid computation, set to approximate value
#lr_coverage=40 

## Long reads aligned to genome assembly (BAM) path, will be produced if not set
#lr_to_assembly_bam_path=/path/to/long/read/alignments

## path to the genome assembly fasta file for 02_predict_rdna_in_genome.sh (also required for 01_predict_rdna.sh if lr_coverage or lr_to_assembly_bam_path variables are not provided)
#assembly_path=/path/to/assembly

## number of threads, use all avaliable CPUs if not specified
#threads=20

## similarity e-value cut-off (default '1e-06'); all HMM predictions with higher e-value will be discarded
#evalue='1e-06'

## rdna gene estimated size (with intergenic spacer) for assembly (in x0.5-x2 range, defalut value 30k)
#rdna_repeat_size='30k'

## type of long reads, possible values: 'pacbio-raw', 'nano-raw' 
#lr_type=pacbio-raw

## number of polishing iterations after assembly, default is 3
#n_polish_iterations=3

## rDNA operon fasta file, required for 02_predict_rdna_in_genome.sh; this file is normally a result of 01_predict_rdna.sh script
#rdna_operon_path=/path/to/rdna/operon/fasta


# tools excecutable variables: use only if the tools are not installed (i.e., not in the PATH),
# or want to use a custom build

## barrnap ( https://github.com/tseemann/barrnap )
#barrnap=/path/to/barrnap

## samtools ( https://github.com/samtools/samtools )
#samtools=/path/to/samtools

## bedtools >= 2.27.0 ( https://github.com/arq5x/bedtools2/releases )
#bedtools=/path/to/bedtools

## minimap2 ( https://github.com/lh3/minimap2 ) - needed if long read alignments to assembly not provided and long read coverage variable not set
#minimap2=/path/to/minimap2

## nhmmer (part of HMMER 3.x: http://hmmer.org/ )
#nhmmer=/path/to/nhmmer

## flye ( https://github.com/fenderglass/Flye )
#flye=/path/to/flye

## BLAST+: blastn and makeblastdb (provide only if it is not in the same directory as blastn) ( ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ) - for 02_predict_rdna_in_genome.sh
#blastn=/path/to/blastn
#makeblastdb//path/to/makeblastdb
