/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_GENOME: Download genomes and normalize aliases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOCLI_DOWNLOAD } from '../../../modules/local/annocli_download'
include { ANNOCLI_ALIAS    } from '../../../modules/local/annocli_alias'

workflow PREPARE_GENOME {

    take:
    ch_species   // channel: [ meta, assembly, annotation ]

    main:

    ch_versions = Channel.empty()

    //
    // Separate species that need download from those with pre-provided files
    //
    ch_species
        .branch {
            meta, assembly, annotation ->
                needs_download: !assembly
                ready:          true
        }
        .set { ch_branched }

    //
    // MODULE: Download genome + annotation via annocli (when files not provided)
    //
    if (!params.skip_download) {
        ANNOCLI_DOWNLOAD (
            ch_branched.needs_download.map { meta, assembly, annotation -> [ meta ] }
        )
        ch_versions = ch_versions.mix(ANNOCLI_DOWNLOAD.out.versions.first())

        // Combine downloaded files with metadata
        ch_downloaded = ANNOCLI_DOWNLOAD.out.fasta
            .join(ANNOCLI_DOWNLOAD.out.gff)
            .map { meta, fasta, gff -> [ meta, fasta, gff ] }

        // Merge downloaded and pre-provided species
        ch_with_files = ch_branched.ready.mix(ch_downloaded)
    } else {
        ch_with_files = ch_branched.ready
    }

    //
    // MODULE: Normalize sequence aliases between GFF and FASTA
    //
    ch_for_alias = ch_with_files
        .filter { meta, assembly, annotation -> annotation }
        .map { meta, assembly, annotation -> [ meta, annotation, assembly ] }

    ANNOCLI_ALIAS ( ch_for_alias )
    ch_versions = ch_versions.mix(ANNOCLI_ALIAS.out.versions.first())

    // Combine alias-matched GFF with assembly
    ch_prepared = ANNOCLI_ALIAS.out.gff
        .join(ch_with_files.map { meta, assembly, annotation -> [ meta, assembly ] })
        .map { meta, gff, assembly -> [ meta, assembly, gff ] }

    // For target-only species (no annotation), pass through assembly
    ch_targets_only = ch_with_files
        .filter { meta, assembly, annotation -> !annotation }

    emit:
    prepared      = ch_prepared       // channel: [ meta, assembly, alias_matched_gff ]
    targets_only  = ch_targets_only   // channel: [ meta, assembly, [] ]
    versions      = ch_versions
}
