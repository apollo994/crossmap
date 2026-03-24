process RELOCATE_LOCI {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(intergenic_bed), path(nonoverlapping_gff)

    output:
    tuple val(meta), path("*.decoy.gff3"), emit: gff
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO: Rewrite relocate_loci2.py logic here
    // The script takes intergenic BED regions and non-overlapping lncRNA GFF3,
    // then relocates lncRNA gene models into intergenic regions to create decoy annotations.
    """
    echo "ERROR: relocate_loci logic not yet implemented — awaiting user-provided script" >&2
    echo "##gff-version 3" > ${prefix}.decoy.gff3
    exit 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.decoy.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12
    END_VERSIONS
    """
}
