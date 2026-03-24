process GET_INTERGENIC_INTERVALS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::agat=1.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/agat:1.4.2--pl5321hdfd78af_0' :
        'biocontainers/agat:1.4.2--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.intergenic.bed"), emit: bed
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 1) Convert ncRNA_gene and pseudogene -> gene for AGAT compatibility
    awk -F'\\t' -v OFS='\\t' '
        /^#/ { print; next }
        NF < 9 { print; next }
        {
            if (\$3 == "ncRNA_gene" || \$3 == "pseudogene") \$3 = "gene"
            print
        }
    ' ${gff} > fixed_types.gff3

    # 2) Add intergenic regions with AGAT
    agat_sp_add_intergenic_regions.pl \\
        --gff fixed_types.gff3 \\
        --out with_intergenic.gff3

    # 3) Extract intergenic_region features and convert to BED6
    awk -F'\\t' -v OFS='\\t' '
        function attr_val(attrs, key,   n,i,a,kv,k,v) {
            n = split(attrs, a, /;/)
            for (i=1; i<=n; i++) {
                split(a[i], kv, /=/)
                k = kv[1]; v = kv[2]
                if (k == key) return v
            }
            return ""
        }
        BEGIN { c=0 }
        /^#/ { next }
        NF < 9 { next }
        \$3 != "intergenic_region" { next }
        {
            c++
            chrom = \$1
            start0 = \$4 - 1
            end = \$5
            id = attr_val(\$9, "ID")
            nm = attr_val(\$9, "Name")
            name = (id != "" ? id : (nm != "" ? nm : ("intergenic_region_" c)))
            score = 0
            strand = "."
            print chrom, start0, end, name, score, strand
        }
    ' with_intergenic.gff3 > ${prefix}.intergenic.bed

    # Clean up
    rm -f fixed_types.gff3 with_intergenic.gff3 *_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: \$(agat_sp_add_intergenic_regions.pl --version 2>&1 | grep -oP 'v\\K[\\d.]+' || echo 'unknown')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.intergenic.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        agat: 1.4.2
    END_VERSIONS
    """
}
