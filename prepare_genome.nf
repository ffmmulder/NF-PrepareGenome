#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/* Example
pipeline="/hpc/ubec/tools/pipelines/dev/NF-PrepareGenome/"
genome_path="/hpc/ubec/resources/genomes/GRCh38_GCA_000001405.15_no_alt/"
/hpc/ubec/tools/nextflow/21.04.3.5560/nextflow run ${pipeline}/prepare_genome.nf --genome_fasta ${genome_path}/Sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --out_dir . -profile slurm --work ./work -c  ${pipeline}/config/nf_manual.config --star_path /hpc/ubec/resources/tools/star/2.7.3a/GRCh38_GCA_000001405.15_no_alt_ncbigtf --annotation ${genome_path}/Annotation/GRCh38_latest_genomic.gtf 

/hpc/ubec/tools/nextflow/21.04.3.5560/nextflow run /hpc/ubec/tools/pipelines/dev/NF-PrepareGenome/prepare_genome.nf --genome_fasta genome.fa --out_dir . -profile slurm --work ./work -c /hpc/ubec/tools/pipelines/dev/NF-PrepareGenome/config/nf_manual.config --bwa true --intervalList true --len true --dict true
*/

/*===========================================================
                        Prepare Genome
===========================================================
#### Homepage / Documentation

----------------------------------------------------------------------------------------
*/
def helpMessage() {
    // Log colors ANSI codes
    c_reset = "\033[0m";
    c_green = "\033[0;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";

    log.info"""
    Usage:
      The typical command for running the pipeline is as follows:
      nextflow run prepare_genome.nf -c </path/to/run.config> --fasta <fasta_dir> --out_dir <output_dir> -profile slurm

      ${c_blue}Mandatory arguments:${c_reset}
       ${c_yellow}--genome_fasta [str]${c_reset}       						Path to a the genome fasta file to be used for genome preparation
       ${c_yellow}--out_dir [str]${c_reset}                                                     The output directory where results will be saved.
       ${c_yellow}-c [str]${c_reset}                                                            The path to the run.config file. ${c_green}Technically all mandatory and optional arguments can be in the run.config file.${c_reset}
       ${c_yellow}-profile [str]${c_reset}                                                      The preconfigured environment to run on, choose from sge (legacy) or slurm (recommended).

      ${c_blue}Optional workflow arguments:${c_reset}
       ${c_yellow}--subset [string]${c_reset}                                                   Generate subset based on given chr list (Space seperated string with the chromosomes to be included, for example: '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X Y MT') for the fasta. (Default: empty)
       ${c_yellow}--chrFiles [bool]${c_reset}                                            	Split fasta in chr subfiles (will be placed in chr_files subdir) (Default: true)
       ${c_yellow}--bwa [bool]${c_reset}                                          		Run BWA indexing. (Default: true)
       ${c_yellow}--bed [bool]${c_reset}                                                  	Create BED file. (Default: true)
       ${c_yellow}--dict [bool]${c_reset}                                                 	Create sequence dictionary .dict file. (Default: true)
       ${c_yellow}--intervalList [bool]${c_reset}                                               Creates interval_list file. (Default: true)
       ${c_yellow}--len [bool]${c_reset}                                                 	Create .len file. (Default: true)
       ${c_yellow}--star [bool]${c_reset}                                               	Generate STAR genome.
       ${c_yellow}--annotation [string]${c_reset}                                               Path to annotation GFF or GTF file, needed for STAR.
       ${c_yellow}--overwrite [bool]${c_reset}                                                  Overwrite existing files when they exist. (Default: false)

      ${c_blue}Optional nextflow arguments:${c_reset}
       ${c_yellow}-resume [str]${c_reset}                                                       Resume a previous nf-iap run.
       ${c_yellow}-N [str]${c_reset}                                                            Send notification to this email when the workflow execution completes.
    """.stripIndent()
//${c_yellow}--star_path [string]${c_reset}                                               	Path to star folder where STAR genome must be generated. When empty this step is skipped.
//       ${c_yellow}--star [bool]${c_reset}                                               	Generate STAR genome, only when --star_path is also supplied (and --star_path being supplied always generates STAR genome).
       
}

/*=================================
          Input validation
=================================*/


// Show help message and exit.
if(params.help){
  helpMessage()
  exit 0
}

if ( !params.genome_fasta){
  exit 1, "Please provide a 'genome_fasta'. You can provide these parameters either in the <analysis_name>.config file or on the commandline (add -- in front of the parameter)."
}

if (!params.out_dir){
  exit 1, "No 'out_dir' parameter found in <analysis_name>.config file or on the commandline (add -- in front of the parameter)."
}

if( params.star ){
    if (!params.annotation){
        exit 1, "No 'annotation' parameter found in <analysis_name>.config file or on the commandline (add -- in front of the parameter). Required when running STAR genome generation"
    }
    //if (!params.star_path){
    //    exit 1, "No 'star_path' parameter passed. Required when running STAR genome generation and needs to be passed on the command line (--star) because of the order in which the flow is processed..."
    //}
    //if (params.star_path.charAt(0) != '/'){
    //    print "Please make sure a full starpath is specified! (/path/to/star)... Otherwise it will be placed in the working directory (NOT the specified output directory if any)."
    //}
} else {
  //  if (params.star_path){
  //      print "star_path specified, overriding star being set to false..."
  //  }
}


/*=================================
          Run workflow
=================================*/
workflow {
  main :
    genome_fasta = Channel
        .fromPath( params.genome_fasta, checkIfExists: true )
        .ifEmpty { exit 1, "FASTA file not found: ${params.genome_fasta}"}    

   genome_index = Channel
        .fromPath( params.genome_fasta+".fai", checkIfExists: true )
        .ifEmpty { 
            println "No "+params.genome_fasta+".fai exists, generating one..."
            include { Faidx } from './NextflowModules/Samtools/1.10/Faidx.nf' params( params )
            Faidx( genome_fasta, "", "")
            genome_index = Faidx.out.genome_faidx
        }    

    if ( params.chr_files ){
        include { CreateChrFiles } from './NextflowModules/Utils/CreateChrFiles.nf' params( params )
        CreateChrFiles( genome_fasta )
    }
    if ( params.len ){
        include { CreateLen } from './NextflowModules/Utils/CreateLen.nf' params( params )
        CreateLen( genome_index )
	genome_len = CreateLen.out.genome_len
    }

    if ( params.dict ){
        include { CreateSequenceDictionary } from './NextflowModules/Picard/2.22.0/CreateSequenceDictionary.nf' params( params )
        CreateSequenceDictionary( genome_fasta )
        genome_dict = CreateSequenceDictionary.out.genome_dict
    }

    if ( params.bed ){
        include { CreateBed } from './NextflowModules/Utils/CreateBed.nf' params( params )
        CreateBed( genome_index)
	genome_bed = CreateBed.out.genome_bed
    }

    if ( params.intervalList ){
        include { CreateIntervalList } from './NextflowModules/Utils/CreateIntervaList.nf' params( params )
        CreateIntervalList( genome_index, genome_dict)
	genome_interval_list = CreateIntervalList.out.genome_interval_list
    }

    //if ( params.star_path ){
    if ( params.star ){
        include { GenomeGenerate } from './NextflowModules/STAR/2.7.3a/GenomeGenerate.nf' params( params : "$params")    

        gtf_file = Channel
            .fromPath( params.annotation, checkIfExists: true )
            .ifEmpty { exit 1, "Annotation file not found: ${params.genome_fasta}"}  

        if ( "gff".equals(params.annotation.split("\\.")[-1]) ){
            include { ConvertGffToGtf } from './NextflowModules/AGAT/0.8.1/ConvertGffToGtf.nf' params( params )
            ConvertGffToGtf( gtf_file )
            GenomeGenerate ( genome_fasta, ConvertGffToGtf.out.gtf )
        } else {         
            GenomeGenerate ( genome_fasta, Channel.fromPath(params.annotation) )
        }
    }

    if ( params.bwa ){
        include { Index } from './NextflowModules/BWA/0.7.17/Index.nf' params(optional: "${params.bwaindex.optional}", genome_fasta : "${params.genome_fasta}")
        Index(genome_fasta)
    }    

}
