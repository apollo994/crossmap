process GET_NON_OVERLAPPING {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::agat=1.4.2 bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-agat-bedtools:latest' :
        'biocontainers/mulled-v2-agat-bedtools:latest' }"

    input:
    tuple val(meta), path(query_gff), path(reference_gff)

    output:
    tuple val(meta), path("*.nonoverlapping.gff3"), emit: gff
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 1) Extract gene-level features from query GFF to BED
    awk -F'\\t' 'BEGIN{OFS="\\t"}
        \$0 !~ /^#/ && (\$3=="gene" || \$3=="pseudogene" || \$3=="ncRNA_gene") {
            id=""
            if (match(\$9, /ID=([^;]+)/, m)) id=m[1]
            if (id=="") next
            print \$1, \$4-1, \$5, id, 0, \$7
        }' ${query_gff} | sort -k1,1 -k2,2n > query_genes.bed

    # 2) Extract gene-level features from reference GFF to BED
    awk -F'\\t' 'BEGIN{OFS="\\t"}
        \$0 !~ /^#/ && (\$3=="gene" || \$3=="pseudogene" || \$3=="ncRNA_gene") {
            id=""
            if (match(\$9, /ID=([^;]+)/, m)) id=m[1]
            if (id=="") next
            print \$1, \$4-1, \$5, id, 0, \$7
        }' ${reference_gff} | sort -k1,1 -k2,2n > reference_genes.bed

    # 3) Keep query genes with ZERO overlap against reference
    bedtools intersect -v \\
        -a query_genes.bed \\
        -b reference_genes.bed \\
        > nonoverlap_genes.bed || true

    # 4) If none remain, write empty GFF
    if [ ! -s nonoverlap_genes.bed ]; then
        printf "##gff-version 3\\n" > ${prefix}.nonoverlapping.gff3
        echo "INFO: No non-overlapping genes found." >&2
    else
        # 5) Build keep list of gene IDs
        cut -f4 nonoverlap_genes.bed | sort -u > keep_ids.txt

        # 6) Use AGAT to filter — keeps matching genes and all children
        agat_sp_filter_feature_from_keep_list.pl \\
            --gff ${query_gff} \\
            --keep_list keep_ids.txt \\
            -o ${prefix}.nonoverlapping.gff3

        # Clean up AGAT report files
        rm -f *_report.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        agat: \$(agat_sp_filter_feature_from_keep_list.pl --version 2>&1 | grep -oP 'v\\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "##gff-version 3" > ${prefix}.nonoverlapping.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
        agat: 1.4.2
    END_VERSIONS
    """
}
