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
      --background                  Path to non-significant CAGE peaks
      -profile                      Hardware config to use. docker / aws

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
 * SETP CONFIGURATION VARIABLES
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
params.range = 300
//params.pfms = "/JASPAR2018_CORE_vertebrates_nr_pfms.homer"
params.pfms = false
params.pfms_jaspar = false
params.background = false
params.encode = "/merged_concat_tfs/"


//output_docs = file("$baseDir/docs/output.md")

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}



// Header log info
log.info "========================================="
log.info " TF-activity v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Fasta Ref']    = params.fasta
summary['Motif File']   = params.pfms
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
// Load BED file
if ( params.peaks ){
    Channel.fromPath(params.peaks).ifEmpty{ exit 1, "Cannot find peak file"}.into{ peak_file_extraction; peak_file_annotation}
} else {
    exit 1, "Specify the peak file!"
}
// Load background
if ( params.background ){
    Channel.fromPath(params.background).ifEmpty{exit 1, "Cannot find background file"}.set{background}
}

// Load FASTA file
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
// Load GTF file
if( params.gtf ){
    gtf = file(params.gtf)
    if ( !gtf.exists() ) exit 1, "GTF file not found:" ${params.gtf}
}
// Load pfms
if (!params.pfms_jaspar && !params.pfms){
    // Use default PFMs
    process copy_pfms {
        publishDir "${params.outdir}/motifs", mode: 'copy'

        output:
        file "*homer" into pfms

        script:
        """
        cp "/JASPAR2018_CORE_vertebrates_nr_pfms.homer" .
        """
    }
}
if( params.pfms ){
    pfms = file(params.pfms)
    if ( !pfms.exists() ) exit 1, "PFMs not found:" ${params.pfms}
}
// If Jaspar motifs given, transform them to Homer format
else if (params.pfms_jaspar ) {
    pfms_jaspar = file(params.pfms_jaspar)
    if ( !pfms_jaspar.exists()) exit 1, "Jaspar PFMs not found:" ${params.pfms_jaspar}

    process jaspar_to_homer {
        publishDir "${params.outdir}/motifs", mode: 'copy'

        input:
        file pfms_jaspar from pfms_jaspar

        output:
        file "*.homer" into pfms

        script:
        """
        jaspar_to_homer_motif.py $pfms_jaspar
        """
    }
}

/**
 * STEP 1 Extract the extended regions
 */
process extract_regions {
    tag "$peaks.baseName"
    publishDir "${params.outdir}/regions", mode: 'copy'

    input:
    file peaks from peak_file_extraction
    file fasta from fasta

    output:
    file "*.fasta" into extended_peaks

    script:
    """
    extract_regions.py $peaks -g $fasta -r $params.range
    """
}
// Create different output channels
extended_peaks.into{ ext_peaks_background; ext_peaks_enrichment; ext_peaks_bed; ext_peaks_targets }

/**
 * STEP 2 Make background sequences
 * If non-DE peaks are availble, extend their regions
 * Else shuffle the input peaks
 */
if (params.background){
    process make_background {
        publishDir "${params.outdir}/background", mode: 'copy'

        input:
        file peaks from background
        file fasta from fasta

        output:
        file "*.fasta" into background_seq

        script:
        """
        extract_regions.py $peaks -g $fasta -r $params.range -f background_seq.fasta
        """
    }
}
else {
    process make_shuffled_background {
        tag "$ext_peaks.baseName"
        publishDir "${params.outdir}/background", mode: 'copy'

        input:
        file ext_peaks from ext_peaks_background

        output:
        file "*.fasta" into background_seq

        script:
        """
        fasta-shuffle-letters -dna $ext_peaks shuffled_background_seq.fasta
        """
    }
}


/**
 * STEP 3 Calculate the enrichment
 */
process enrichment {
    tag "$peaks.baseName"
    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    file peaks from ext_peaks_enrichment
    file background from background_seq
    file pfms from pfms

    output:
    file "*.txt" into enrichment_result

    script:
    """
    cat $pfms > motifs_used.pfm
    homer2 known -i $peaks -b $background -m $pfms -opt -stat hypergeo -p ${task.cpus} > homer2_enrichment_result.txt
    """
}

/**
 * STEP 4 Create BED files from extended region
 */
process region_bed {
    tag "${peaks.baseName}"
    publishDir "${params.outdir}/regions", mode: 'copy'

    input:
    file peaks from ext_peaks_bed

    output:
    file "*.bed" into region_bed
    file "*.merged.bed" into merged_region_bed
    file "*sorted.bed" into sorted_peaks_bed

    script:
    """
    region_to_bed.py $peaks ${peaks.baseName}.bed -r $params.range
    bedtools sort -i ${peaks.baseName}.bed > ${peaks.baseName}.sorted.bed
    bedtools merge -i ${peaks.baseName}.sorted.bed > ${peaks.baseName}.merged.bed
    """
}

/**
 * STEP 5 Calculate intersection with TF ChIP-seq peaks from ENCODE
 */
process encode_intersect {
    tag "${region_bed.baseName}"
    publishDir "${params.outdir}/encode", mode: 'copy'

    input:
    file region_bed from merged_region_bed

    output:
    file "*.txt" into encode_intersection

    script:
    """
    intersect_chipseq_data.py $region_bed $params.encode -out ${region_bed.baseName}.encodeIntersect.txt
    """
}

/**
 * STEP 6 Find TF motif instances in sequences
 */
process find_motifs {
    publishDir "${params.outdir}/tf_targets", mode: 'copy'
    cpus 4

    input:
    file peaks from ext_peaks_targets
    file pfms from pfms

    output:
    file "*.txt" into motif_instances

    script:
    """
    homer2 find -i $peaks -m $pfms -p ${task.cpus} > motif_instances_homer2.txt
    """
}

/**
 * STEP 7 Annotate peaks
 */
process annotate {
    tag "${peaks.baseName}"
    publishDir "${params.outdir}/annotate", mode: 'copy'

    input:
    file peaks from peak_file_annotation
    file gtf from gtf

    output:
    file "*annotated_peaks.txt" into annotated_peaks

    script:
    """
    bed_to_peak.py $peaks id_peak_file.txt
    annotatePeaks.pl id_peak_file.txt $params.fasta -gtf $gtf > annotated_peaks.txt
    """
}


/**
 * STEP 8 TF-target filtering
 */
process tf_targets {
    tag "${enriched.baseName}"
    publishDir "${params.outdir}/tf_targets", mode: 'copy'

    input:
    file instances from motif_instances
    file enriched from enrichment_result
    file anno_peaks from annotated_peaks

    output:
    file "*.txt" into tf_target_results

    script:
    """
    tf_targets.py $enriched $instances $anno_peaks
    """
}

/*
 * Completion notification
 */
workflow.onComplete {
    log.info "[TF Activity] Pipeline Complete"
}
