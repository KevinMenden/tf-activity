#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         NF-CAGEseq
========================================================================================
 NF-CAGEseq Analysis Pipeline. Started 2018-03-09.
 #### Homepage / Documentation
 https://github.com/KevinMenden/NF-CAGEseq
 #### Authors
 Kevin Menden KevinMenden <kevin.menden@dzne.de> - https://github.com/KevinMenden>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     TF-activity v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run kevinmenden/tf-activity --peaks 'cage_peaks.bed' -profile docker


    Mandatory arguments:
      --peaks                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. docker / aws

    Options:
      --pairedEnd                   Specifies that the input is paired end reads (Currently not supported)
      --aligner                     [ star ] Specify which aligner should be used. (Currently only 'star' supported)

    Trimming:
      --cutEcop                     [true|false] Whether the 5' EcoP15I regognition site should be removed. Default is true.
      --cutLinker                   [true|false] Whether the 3' linker should be removed. Default is true.
      --trimming                    [true | false] Whether trimming should be performed. Default is true.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF reference
      --genome                      Name of genome

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = params.version

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.outdir = './results'
params.gtf = false
params.saveReference = false


output_docs = file("$baseDir/docs/output.md")





// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}


"""


// Header log info
log.info "========================================="
log.info " NF-CAGEseq v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Fasta Ref']    = params.fasta
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
            "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}

/**
 * Load and validate inputs
 */
// Load FASTA file
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
// Load GTF file
if( params.gtf ){
    Channel
            .fromPath(params.gtf)
            .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
            .into { gtf_makeSTARindex; gtf_star}
}







/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    tag "$name"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion notification
 */
workflow.onComplete {
    log.info "[TF Activity] Pipeline Complete"
}
