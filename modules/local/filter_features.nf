process FILTER_FEATURES {
    tag "$meta.id - $feature_type"
    label 'process_single'

    conda "bioconda::agat=1.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)
    each feature_type

    output:
    tuple val(meta), val(feature_type), path("*.${feature_type}.gff3"), emit: gff
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract IDs of features matching the target type
    awk -F'\\t' -v type="${feature_type}" '
        /^#/ { next }
        \$3 == type {
            if (match(\$9, /ID=([^;]+)/, m)) print m[1]
        }
    ' ${gff} | sort -u > target_ids.txt

    # Check if any features were found
    if [ ! -s target_ids.txt ]; then
        echo "##gff-version 3" > ${prefix}.${feature_type}.gff3
        echo "WARNING: No features of type '${feature_type}' found in ${gff}" >&2
    else
        # Use AGAT to filter — keeps matching features and their hierarchy
        agat_sp_filter_feature_from_keep_list.pl \\
            --gff ${gff} \\
            --keep_list target_ids.txt \\
            -o ${prefix}.${feature_type}.gff3

        # Clean up AGAT report files
        rm -f *_report.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_filter_feature_from_keep_list.pl --version 2>&1 | grep -oP 'v\\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.${feature_type}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: 1.4.2
    END_VERSIONS
    """
}
