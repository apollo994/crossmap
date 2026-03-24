process AGAT_EXTRACT_SEQUENCES {
    tag "$meta.id - $feature_type"
    label 'process_low'

    conda "bioconda::agat=1.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), val(feature_type), path(gff), path(fasta)

    output:
    tuple val(meta), val(feature_type), path("*.exons.fa"), emit: fasta
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-t exon --merge'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    agat_sp_extract_sequences.pl \\
        -g ${gff} \\
        -f ${fasta} \\
        ${args} \\
        -o ${prefix}.${feature_type}.exons.fa

    # Clean up AGAT report files
    rm -f *_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_extract_sequences.pl --version 2>&1 | grep -oP 'v\\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${feature_type}.exons.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: 1.4.2
    END_VERSIONS
    """
}
