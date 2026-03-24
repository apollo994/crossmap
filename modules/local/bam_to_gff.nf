process BAM_TO_GFF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.gff3"), emit: gff
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: Rewrite convert_bam_to_gff logic here
    // The script converts spliced BAM alignments to GFF3 gene model annotations.
    // Needs to reconstruct exon/intron structure from CIGAR strings.
    """
    echo "ERROR: bam_to_gff logic not yet implemented — awaiting user-provided script" >&2
    echo "##gff-version 3" > ${prefix}.gff3
    exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.21
    END_VERSIONS
    """
}
