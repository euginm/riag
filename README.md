# rDNA in a Genome

The rDNA-in-a-Genome pipeline performs a qualitative and quantitative assessment of rDNA from a eukaryotic genome in a semi-automatic way. The pipeline is designed to use long sequencing reads and draft assembly of a newly sequenced genome.

## Pipeline features
- Finds long reads that contain rDNA subunits (rDNA reads)
- Filters out the rDNA reads with mitochondrial/bacterial or non-canonical subunits composition
- Calls an assembler to produce the rRNA gene consensus
- Gives the rDNA copy number lower bound
- Predict rDNA in the genome assembly and maps the rRNA gene consensus to it

## Quick start
Clone the git:
```sh
$ git clone https://github.com/euginm/riab
```
Fill `config.sh` and run `01_predict_rdna.sh`:
```sh
$ ./01_predict_rdna.sh config.sh
```
Review the rRNA gene assembly results, update the `config.sh` and start `02_find_rdna_in_genome.sh`:
```sh
$ ./02_find_rdna_in_genome.sh config.sh
```
## Pipeline description

1. Predict rDNA subunits in every long read.
2. Quantify and evaluate found rDNA reads, pick best candidates for consensus assembly.
3. Assemble rDNA consensus sequence from selected rDNA reads.
4. Estimate the number of rDNA copies based on found rDNA reads the amount and average long read coverage depth.
5. Predict rDNA subunits in genome assembly, find rDNA consensus occurrences.

## Dependencies
* [Barrnap][barrnap]: nhmmer wrapper that provides eukaryotic rRNA HMM-profile and conviniently transform nhmmer output into a GFF annotation file. It requires:
    * [nhmmer][nhmmer]: part of HMMER suite, used for HMM-based prediction of rDNA in sequences.
    * [Perl 5.xx][perl 5.xx]: barrnap wrapper is written in Perl.
    * [bedtools][bedtools]: used to convert nhmmer output to GFF.
* [Flye][flye]: de novo assembler, designed to work with error-prone long reads. 
* [Samtools][samtools]: reads and alignments manipulation.
* [Biopython][biopython]: FASTA/Q parcing and writing.
* [Minimap2][minimap2]: used to map long reads to genome assembly.
* [BLAST+][BLAST+]: `blastn` and `makeblastdb` commands, used to locate rDNA consensus in the genome assembly.


The tools are looked up in `PATH`, or paths to their executables are taken from `config.sh`.
## Input
The input arguments are provided via `config.sh`.
`01_predict_rdna.sh` arguments:

| Argument | Description |
| ------: | ------ |
| **`out_dir`** | Directory for output files, will be created if not found |
| **`lr_dir`** | Directory containing raw long reads (multiple) files in FASTA/Q format, can be gzipped |
| **`tmp_dir`** | Directory for temp files, will be created if not found. **Default**: `TMPDIR` environment variable |
| **`prefix`** | Prefix for all output files. **Default**: 'rdna_pipeline' |
| **`lr_coverage`** | Long reads sequencing depth, used for rDNA copy number estimation |
| **`lr_to_assembly_bam_path`** | BAM file with all long reads mapped to (draft) genome assembly, used to compute PacBio sequencing depth, if `pacbio_coverage` not provided |
| **`assembly_path`** | Path to (draft) genome assembly, used as a reference for mapping long reads and sequencing depth estimation, if neither `pacbio_coverage` nor `pacbio_to_assembly_bam_path` provided |
| **`threads`** | Number of CPUs to use (40 or more recommended). **Default**: number of available CPUs |
| **`evalue`** | Maximum e-value for predicted rDNA subunits to keep them in the analysis. **Default**: '1e-06' |
| **`rdna_repeat_size`** | Approximate length of rDNA repeat (rRNA gene + intergenic spacer). Used by Flye assembler,  **Default**: '30k' |
| **`n_polish_iterations`** | Number of assembly polishing iterations, Flye parameter. **Default**: 3 |
| **`lr_type`**| Parameter for Flye assembler, possible values: 'pacbio-raw', 'nano-raw'. **Default**: 'pacbio-raw' |

`02_find_rdna_in_genome.sh` ***additional*** arguments:

| Argument | Description |
| ------: | ------ |
| **`rdna_operon_path`** | Path to fasta file with rDNA operon. Normally, it is produced in previous pipeline stage. |

## Output
All output files are placed in `out_dir` and start with `prefix`. `01_predict_rdna.sh` output:

| Output file | Description
| ------: | ------
| rdna_reads.gff | Annotation file with rDNA subunits predictions, barrnap output.
| rdna_reads.stats | Summary of found rDNA reads.
| rdna_reads.json | JSON file with rDNA read IDs, grouped by rDNA subunits composition.
| rdna_reads_for_assembly.fasta | Long reads, selected for rRNA gene consensus assembly.
| rdna_flye_assembly.fasta | Flye assembler output.
| rdna_flye_assembly.gff | Predicted rDNA subunits in Flye assembly

`02_find_rdna_in_genome.sh` output:

| Output file | Description
| ------: | ------
| (*genome assembly basename*)_rdna_prediction.gff | Annotation file with rDNA subunit predictions in genome assembly, barrnap output.
| rdna_operon_to_(*genome assembly basename*).tab | BLAST output with all occurrences of rDNA operon in genome assembly.
## Notes
* This pipeline is designed to provide an overview of rDNA in an organism. However, it can't substitute a proper rDNA analysis.
* It worth to review `rdna_reads.stats` file. A big amount of discarded rDNA reads may indicate that the organism has unconventional rDNA subunits composition, or HMM-based rDNA predictions do not work properly.
* If the pipeline exits due to error, it is possible to continue it, by using the same config file.

## Example run
In the example run, we use genomic data and genome assembly of a [leatherback sea turtle][vgp_turtle], provided by [Vertebrate Genomes Project][vgp]) 

Download assembly and PacBio reads:
```sh
$ mkdir Dermochelys_coriacea && cd Dermochelys_coriacea
$ wget https://s3.amazonaws.com/genomeark/species/Dermochelys_coriacea/rDerCor1/assembly_curated/rDerCor1.pri.cur.20190930.f asta.gz && gunzip rDerCor1.pri.cur.20190930.fasta.gz
$ mkdir pacbio && cd pacbio
$ aws s3 --no-sign-request sync s3://genomeark/species/Dermochelys_coriacea/rDerCor1/genomic_data/pacbio/ . --exclude "*ccs.bam*"
```
Reads are in .bam format, convert them to fasta:
```sh
$ parallel "bam2fasta -o {} {}.subreads.bam" ::: `ls *.bam | cut -f1 -d.`
$ rm *.bam*
```
Copy `config.sh` to `turtle_config.sh`
```sh
$ cd ../
$ cp ../config.sh turtle_config.sh
```
Fill the `turtle_config.sh` and start rDNA pipeline:
```sh
$ ./01_predict_rdna.sh Dermochelys_coriacea/turtle_config.sh
```
Review the assembly and rDNA predictions in it, cut the rDNA operon:
```sh
$ samtools faidx Dermochelys_coriacea/rDNA_pipeline/sea_turtle_rdna_flye_assembly.fasta contig_1:69100-77700 > Dermochelys_coriacea/rDNA_pipeline/sea_turtle_rdna_operon_from_flye_assembly.fasta
```
Update `turtle_config.sh` with rDNA operon path and continue the pipeline:
```sh
$ ./02_find_rdna_in_genome.sh Dermochelys_coriacea/turtle_config.sh
```

[python]: <https://www.python.org/downloads>
[perl 5.XX]: <https://www.perl.org/get.html>
[barrnap]: <https://github.com/tseemann/barrnap>
[samtools]: <https://github.com/samtools/samtools>
[minimap2]: <https://github.com/lh3/minimap2>
[bedtools]: <https://github.com/arq5x/bedtools2/releases>
[nhmmer]: <http://hmmer.org/>
[flye]: <https://github.com/fenderglass/Flye>
[BLAST+]: <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>
[biopython]: <https://biopython.org/wiki/Download>
[vgp_turtle]: <https://vgp.github.io/genomeark/Dermochelys_coriacea>
[vgp]: <https://vertebrategenomesproject.org>
