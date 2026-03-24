/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { MULTIQC                } from '../modules/nf-core/multiqc/main'

// Local modules
include { BAM_TO_GFF             } from '../modules/local/bam_to_gff'

// Subworkflows
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome'
include { EXTRACT_FEATURES       } from '../subworkflows/local/extract_features'

// nf-core utilities
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_crossmap_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CROSSMAP {

    take:
    ch_samplesheet // channel: [ meta, assembly, annotation ]
    ch_sources     // channel: [ meta, assembly, annotation ] — species providing gene models
    ch_targets     // channel: [ meta, assembly, annotation ] — species receiving gene models

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Prepare genomes — download (optional) + alias normalization
    //
    // TODO: Enable PREPARE_GENOME once annocli is available
    // PREPARE_GENOME ( ch_sources )
    // ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // For now, use source files directly from samplesheet
    // (assumes annotation is already alias-matched or alias step will be added)
    ch_prepared_sources = ch_sources

    //
    // SUBWORKFLOW: Extract features from source annotations
    //
    // Parse feature_types parameter into a channel list
    ch_feature_types = Channel.fromList(params.feature_types.tokenize(','))

    EXTRACT_FEATURES (
        ch_prepared_sources,
        ch_feature_types
    )
    ch_versions = ch_versions.mix(EXTRACT_FEATURES.out.versions)

    //
    // Step 8: MINIMAP2 — all-on-all spliced alignment
    // Create cartesian product: source features × target genomes
    //
    // source features: [ meta_source, feature_type, exon_fasta ]
    // target genomes:  [ meta_target, assembly, annotation ]
    //
    ch_source_features = EXTRACT_FEATURES.out.feature_sequences

    ch_target_genomes = ch_targets
        .map { meta, assembly, annotation -> [ meta, assembly ] }

    // Cartesian product of source features × target genomes
    ch_mapping_input = ch_source_features
        .combine(ch_target_genomes)
        .map { meta_source, feature_type, exon_fasta, meta_target, target_assembly ->
            def meta_mapping = [
                id:             "${meta_target.id}_${feature_type}_${meta_source.id}",
                source_id:      meta_source.id,
                source_taxid:   meta_source.taxid,
                target_id:      meta_target.id,
                target_taxid:   meta_target.taxid,
                feature_type:   feature_type,
                target_role:    meta_target.role
            ]
            return [ meta_mapping, target_assembly, exon_fasta ]
        }

    // TODO: Step 8 — Install and call MINIMAP2_ALIGN nf-core module
    //   MINIMAP2_ALIGN ( ch_mapping_input.map { meta, ref, reads -> [ meta, reads, ref, true, 'bam', false, false ] } )
    //   ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    // TODO: Step 9 — BAM_TO_GFF
    //   BAM_TO_GFF ( MINIMAP2_ALIGN.out.bam )
    //   ch_versions = ch_versions.mix(BAM_TO_GFF.out.versions.first())

    // TODO: Step 10 — GFFCOMPARE (conditional: when target is also a source)
    //   ch_for_gffcompare = BAM_TO_GFF.out.gff
    //       .filter { meta, gff -> meta.target_role == 'both' }
    //   GFFCOMPARE ( ch_for_gffcompare joined with reference GFF )

    // TODO: Step 11 — GFFREAD_EXTRACT (extract mapped transcript sequences)
    //   GFFREAD ( BAM_TO_GFF.out.gff joined with target assembly )

    // TODO: Step 12 — BUSCO (source + mapped transcripts)
    //   BUSCO ( source_transcripts.mix(mapped_transcripts) )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'crossmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
