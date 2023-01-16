# NF-PrepareGenome

This pipeline allows the user to generate various genome resource datafiles required by the different pipelines used at the UMCU/UBEC such as NF-IAP, RNASeq-NF and Sarek.

The absolute minimum requirement is the fasta file for the genome for which these files are generated, but depending on the steps to execute additional files might be required (although usually at most an annotation file in the form of a GTF or GFF file). When there is no corresponding index for the fasta file present, the pipeline will generate it.

Various configurations are available specicying the required steps for each pipeline, but also for simply running all steps or none of them with the latter requiring the users to specify the steps themselves.

## Installing and setup

The following tools and steps are required for this pipeline:

1. [Nextflow v21+](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. [Singularity](https://sylabs.io/guides/3.5/admin-guide/)
3. [Pull/Clone NF-PrepareGenome](#pull-or-clone)

## Pull or Clone (clone NextflowModules seperately)
```
git clone git@github.com:UMCUGenetics/NF-PrepareGenome.git
git -C ./NF-PrepareGenome clone git@github.com:UMCUGenetics/NextflowModules.git
```

## Running the pipeline
```
# Minimal example
nextflow run prepare_genome.nf --genome_fasta /path/Sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --out_dir . -profile slurm --work ./work -c config/nf_iap.config

# Example including GTF file
nextflow run prepare_genome.nf --genome_fasta /path/Sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --out_dir . -profile slurm --work ./work -c config/nf_manual.config --star_path /path/to/star_genome/folder/ --annotation /path/to/Annotation/GRCh38_latest_genomic.gtf 
```

## Options

The following options are vailable when running the pipeline:

**--chr_files**  
These are used by Control-FREEC

**--len**  
Length file

**--dict**  
Dictionary file, used by GATK

**--bed**  
BED file

**--intervalList**  
intervalList file

**--star**  
STAR genome generation

**--bwa**  
BWA indices generation

## TODO

Generation of mappability tracks for Control-FREEC
