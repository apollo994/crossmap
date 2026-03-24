process AGAT_LONGEST_ISOFORM {
    tag "$meta.id - $feature_type"
    label 'process_single'

    conda "bioconda::agat=1.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), val(feature_type), path(gff)

    output:
    tuple val(meta), val(feature_type), path("*.longest.gff3"), emit: gff
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_sp_keep_longest_isoform.pl \\
        -g ${gff} \\
        -o ${prefix}.${feature_type}.longest.gff3

    # Clean up AGAT report files
    rm -f *_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_keep_longest_isoform.pl --version 2>&1 | grep -oP 'v\\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.${feature_type}.longest.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: 1.4.2
    END_VERSIONS
    """
}
