/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXTRACT_FEATURES: Filter, deduplicate, and extract sequences from source annotations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Steps:
      3. FILTER_FEATURES — extract feature types (lnc_RNA, mRNA) from GFF3
      4. AGAT_LONGEST_ISOFORM — keep longest isoform per gene
      5. GET_NON_OVERLAPPING — lncRNA not overlapping mRNA (when both types present)
      6. GENERATE_DECOYS — intergenic intervals + relocate loci (optional)
      7. AGAT_EXTRACT_SEQUENCES — GFF3 + genome → exon FASTA
*/

include { FILTER_FEATURES          } from '../../../modules/local/filter_features'
include { AGAT_LONGEST_ISOFORM     } from '../../../modules/local/agat_longest_isoform'
include { AGAT_EXTRACT_SEQUENCES   } from '../../../modules/local/agat_extract_sequences'
include { GET_NON_OVERLAPPING      } from '../../../modules/local/get_non_overlapping'
include { GET_INTERGENIC_INTERVALS } from '../../../modules/local/get_intergenic_intervals'
include { RELOCATE_LOCI            } from '../../../modules/local/relocate_loci'

workflow EXTRACT_FEATURES {

    take:
    ch_sources       // channel: [ meta, assembly, annotation_gff ]
    feature_types    // list: [ 'lnc_RNA', 'mRNA' ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 3: FILTER_FEATURES — extract each feature type from source annotations
    //
    ch_source_gff = ch_sources.map { meta, assembly, gff -> [ meta, gff ] }

    FILTER_FEATURES (
        ch_source_gff,
        feature_types
    )
    ch_versions = ch_versions.mix(FILTER_FEATURES.out.versions.first())

    //
    // Step 4: AGAT_LONGEST_ISOFORM — keep longest isoform per gene
    //
    AGAT_LONGEST_ISOFORM ( FILTER_FEATURES.out.gff )
    ch_versions = ch_versions.mix(AGAT_LONGEST_ISOFORM.out.versions.first())

    // Collect longest isoform GFFs: [ meta, feature_type, longest_gff ]
    ch_longest = AGAT_LONGEST_ISOFORM.out.gff

    //
    // Step 5: GET_NON_OVERLAPPING — lncRNA genes not overlapping mRNA genes
    // Only when both lnc_RNA and mRNA feature types are extracted
    //
    // Split by feature type to pair lncRNA with mRNA per species
    ch_lncrna = ch_longest
        .filter { meta, feature_type, gff -> feature_type == 'lnc_RNA' }
        .map { meta, feature_type, gff -> [ meta.id, meta, gff ] }

    ch_mrna = ch_longest
        .filter { meta, feature_type, gff -> feature_type == 'mRNA' }
        .map { meta, feature_type, gff -> [ meta.id, gff ] }

    // Join lncRNA and mRNA by sample ID, then compute non-overlapping
    ch_for_nonoverlap = ch_lncrna
        .join(ch_mrna)
        .map { id, meta, lncrna_gff, mrna_gff -> [ meta, lncrna_gff, mrna_gff ] }

    GET_NON_OVERLAPPING ( ch_for_nonoverlap )
    ch_versions = ch_versions.mix(GET_NON_OVERLAPPING.out.versions.first())

    // Add non-overlapping as another feature variant
    ch_nonoverlapping = GET_NON_OVERLAPPING.out.gff
        .map { meta, gff -> [ meta, 'lnc_RNA_nonoverlapping', gff ] }

    //
    // Step 6: GENERATE_DECOYS (optional — controlled by params.skip_decoy)
    //
    ch_decoy_gff = Channel.empty()
    if (!params.skip_decoy) {
        // 6a: Get intergenic intervals from full source annotation
        ch_source_full_gff = ch_sources.map { meta, assembly, gff -> [ meta, gff ] }
        GET_INTERGENIC_INTERVALS ( ch_source_full_gff )
        ch_versions = ch_versions.mix(GET_INTERGENIC_INTERVALS.out.versions.first())

        // 6b: Relocate non-overlapping lncRNAs into intergenic regions
        ch_for_relocate = GET_INTERGENIC_INTERVALS.out.bed
            .join(GET_NON_OVERLAPPING.out.gff)

        RELOCATE_LOCI ( ch_for_relocate )
        ch_versions = ch_versions.mix(RELOCATE_LOCI.out.versions.first())

        ch_decoy_gff = RELOCATE_LOCI.out.gff
            .map { meta, gff -> [ meta, 'decoy', gff ] }
    }

    //
    // Step 7: EXTRACT_SEQUENCES — GFF3 + genome → exon FASTA for each feature set
    //
    // Combine all GFF variants: original feature types + non-overlapping + decoy
    ch_all_gffs = ch_longest
        .mix(ch_nonoverlapping)
        .mix(ch_decoy_gff)

    // Join with assembly to get [ meta, feature_type, gff, assembly ]
    ch_assemblies = ch_sources.map { meta, assembly, gff -> [ meta.id, assembly ] }

    ch_for_extract = ch_all_gffs
        .map { meta, feature_type, gff -> [ meta.id, meta, feature_type, gff ] }
        .combine(ch_assemblies, by: 0)
        .map { id, meta, feature_type, gff, assembly -> [ meta, feature_type, gff, assembly ] }

    AGAT_EXTRACT_SEQUENCES ( ch_for_extract )
    ch_versions = ch_versions.mix(AGAT_EXTRACT_SEQUENCES.out.versions.first())

    emit:
    feature_sequences = AGAT_EXTRACT_SEQUENCES.out.fasta  // channel: [ meta, feature_type, exon_fasta ]
    all_gffs          = ch_all_gffs                        // channel: [ meta, feature_type, gff ] — all feature variants
    longest_gffs      = ch_longest                         // channel: [ meta, feature_type, longest_gff ]
    nonoverlapping    = GET_NON_OVERLAPPING.out.gff        // channel: [ meta, nonoverlapping_gff ]
    versions          = ch_versions
}
