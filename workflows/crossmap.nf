/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules
include { MINIMAP2_ALIGN         } from '../modules/nf-core/minimap2/align/main'
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

    //
    // Step 8: MINIMAP2_ALIGN — spliced alignment
    //
    MINIMAP2_ALIGN (
        ch_mapping_input.map { meta, ref, reads -> [ meta, reads ] },   // reads (source exon FASTA)
        ch_mapping_input.map { meta, ref, reads -> [ meta, ref ] },     // reference (target genome)
        true,       // bam_format
        'csi',      // bam_index_extension
        false,      // cigar_paf_format
        false       // cigar_bam
    )

    //
    // Step 9: BAM_TO_GFF — convert spliced BAM to GFF3 gene models
    //
    BAM_TO_GFF ( MINIMAP2_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_TO_GFF.out.versions.first())

    //
    // Step 10: GFFCOMPARE — compare mapped vs reference (when target is also source)
    // Only run when target species has role 'both' (has its own annotation)
    //
    // TODO: Wire GFFCOMPARE once nf-core module is installed
    // ch_for_gffcompare = BAM_TO_GFF.out.gff
    //     .filter { meta, gff -> meta.target_role == 'both' && !params.skip_gffcompare }
    // Need to join with the correct reference GFF (matching feature_type for target species)
    // GFFCOMPARE ( ch_for_gffcompare, ch_reference_gff )

    //
    // Step 11: GFFREAD — extract mapped transcript sequences
    // Needs target genome FASTA joined with predicted GFF
    //
    // TODO: Wire GFFREAD once nf-core module is installed
    // ch_for_gffread = BAM_TO_GFF.out.gff
    //     .map { meta, gff -> [ meta.target_id, meta, gff ] }
    //     .combine(ch_target_genomes.map { meta, assembly -> [ meta.id, assembly ] }, by: 0)
    //     .map { target_id, meta, gff, assembly -> [ meta, gff, assembly ] }
    // GFFREAD ( ch_for_gffread )

    //
    // Step 12: BUSCO — completeness assessment (source + mapped transcripts)
    //
    // TODO: Wire BUSCO once nf-core module is installed
    // ch_for_busco = EXTRACT_FEATURES.out.feature_sequences
    //     .map { meta, ft, fasta -> [ meta + [busco_tag: "source_${ft}"], fasta ] }
    //     .mix(GFFREAD.out.fasta.map { meta, fasta -> [ meta + [busco_tag: "mapped"], fasta ] })
    // BUSCO ( ch_for_busco, params.busco_lineage, [], [], [] )

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
