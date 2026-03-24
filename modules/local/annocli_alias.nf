process ANNOCLI_ALIAS {
    tag "$meta.id"
    label 'process_single'

    // TODO: Update container once annocli is containerized
    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(gff), path(fasta)

    output:
    tuple val(meta), path("*.aliasMatch.gff3"), emit: gff
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: Replace with actual annocli alias command once user provides interface
    """
    echo "ERROR: annocli alias not yet implemented — awaiting user-provided tool details" >&2
    cp ${gff} ${prefix}.aliasMatch.gff3
    exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annocli: unknown
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.aliasMatch.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annocli: 0.0.0
    END_VERSIONS
    """
}
