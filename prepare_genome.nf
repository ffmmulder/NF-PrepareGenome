#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*===========================================================
                        Prepare Genome
===========================================================
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

if ( !params.nextflow_modules_path ){
  println "No NextflowModules path specified, assumging './NextflowModules'. Otherwise please specify it using --nextflow_modules_path or by setting it in the nextflow.config file"
  params.nextflow_modules_path = './NextflowModules'
  Channel
        .fromPath( params.nextflow_modules_path, checkIfExists: true )
        .ifEmpty { exit 1, "NextflowModules path not found: ${params.nextflow_modules_path}"}    
  
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
}


/*=================================
          Run workflow
=================================*/
workflow {
  main :
    genome_fasta = Channel
        .fromPath( params.genome_fasta, checkIfExists: true )
        .ifEmpty { exit 1, "FASTA file not found: ${params.genome_fasta}"}    
   
    if( !file(params.genome_fasta+".fai").exists() || params.fai )
    {
       if (! params.fai ){
           println params.genome_fasta+".fai not found, creating!"
       }
       include { Faidx } from params.nextflowmodules_path+'/Samtools/1.10/Faidx.nf' params( params )
       Faidx( genome_fasta )
       genome_index = Faidx.out.genome_faidx
    }   

    if ( params.chr_files ){
        include { CreateChrFiles } from params.nextflowmodules_path+'//Utils/CreateChrFiles.nf' params( params )
        CreateChrFiles( genome_fasta )
    }
    if ( params.len ){
//        include { CreateLen } from './NextflowModules/Utils/CreateLen.nf' params( params )
        include { CreateLen } from params.nextflowmodules_path+'/Utils/CreateLen.nf' params( params )
        CreateLen( genome_index )
	genome_len = CreateLen.out.genome_len
    }

    if ( params.dict ){
        include { CreateSequenceDictionary } from params.nextflowmodules_path+'/Picard/2.22.0/CreateSequenceDictionary.nf' params( params )
        CreateSequenceDictionary( genome_fasta )
        genome_dict = CreateSequenceDictionary.out.genome_dict
    }

    if ( params.bed ){
        include { CreateBed } from params.nextflowmodules_path+'/Utils/CreateBed.nf' params( params )
        CreateBed( genome_index)
	genome_bed = CreateBed.out.genome_bed
    }

    if ( params.intervalList ){
        include { CreateIntervalList } from params.nextflowmodules_path+'/Utils/CreateIntervaList.nf' params( params )
        CreateIntervalList( genome_index, genome_dict)
	genome_interval_list = CreateIntervalList.out.genome_interval_list
    }

    //if ( params.star_path ){
    if ( params.star ){
        include { GenomeGenerate } from params.nextflowmodules_path+'//STAR/2.7.3a/GenomeGenerate.nf' params( params : "$params")    

        gtf_file = Channel
            .fromPath( params.annotation, checkIfExists: true )
            .ifEmpty { exit 1, "Annotation file not found: ${params.genome_fasta}"}  

        if ( "gff".equals(params.annotation.split("\\.")[-1]) ){
            include { ConvertGffToGtf } from params.nextflowmodules_path+'//AGAT/0.8.1/ConvertGffToGtf.nf' params( params )
            ConvertGffToGtf( gtf_file )
            GenomeGenerate ( genome_fasta, ConvertGffToGtf.out.gtf )
        } else {         
            GenomeGenerate ( genome_fasta, Channel.fromPath(params.annotation) )
        }
    }

    if ( params.bwa ){
        include { Index } from params.nextflowmodules_path+'//BWA/0.7.17/Index.nf' params(optional: "${params.bwaindex.optional}", genome_fasta : "${params.genome_fasta}")
        Index(genome_fasta)
    }    

    // Create log files: Repository versions and Workflow params
    VersionLog()
    Workflow_ExportParams()
}

// Workflow completion notification
workflow.onComplete {
    // HTML Template
    def template = new File("$baseDir/assets/workflow_complete.html")
    def binding = [
        runName: run_name,
        workflow: workflow
    ]
    def engine = new groovy.text.GStringTemplateEngine()
    def email_html = engine.createTemplate(template).make(binding).toString()

    // Send email
    if (workflow.success) {
        def subject = "Deeplexicon Workflow Successful: ${run_name}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    } else {
        def subject = "Deeplexicon Workflow Failed: ${run_name}"
        sendMail(to: params.email.trim(), subject: subject, body: email_html)
    }
}

//from other pipeline, modify for this
process VersionLog {
    // Custom process to log repository versions
    tag {"VersionLog ${analysis_id}"}
    label 'VersionLog'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a new version log when restarting the workflow.

    output:
        path('repository_version.log', emit: log_file)

    script:
        """
        echo 'PrepareGenome' > repository_version.log
        git --git-dir=${workflow.projectDir}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log

        echo 'NextflowModules' >> repository_version.log
        git --git-dir=${params.nextflowmodules_path}/.git log --pretty=oneline --decorate -n 2 >> repository_version.log
        """
}