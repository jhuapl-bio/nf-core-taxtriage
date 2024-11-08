process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21 bioconda::samtools=1.12' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads), path(reference)
    val bam_format
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // define mapx if paired end reads or single end illumina or single end oxford, match with regex ignore case
    def mapx = ''
    if (meta.platform =~ /(?i)illumina/) {
        mapx = '-ax sr'
    } else if (meta.platform =~ /(?i)pacbio/) {
        mapx = '-ax map-hifi'
    } else {
        mapx = '-ax map-ont'
    }
    def input_reads = reads.findAll { it != null }.join(' ')
    def I_value = "${(task.memory.toMega() * 0.85).longValue()}M" // 80% of allocated memory, append GB to end as a string
    def S_value = "${( ( task.memory.toMega() / task.cpus ) * 0.1).longValue()}M" // 20% of allocated memory, append GB to end as a string
    // set cpu_view variable to half of the cpus, and round down, minimum of 1
    def cpu_view = Math.max(1, (task.cpus / 2).round().intValue())
        
    def bam_output = bam_format ? "-a | samtools view -@ ${cpu_view} -b -h -o ${prefix}.bam" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    // if input is illumina then use -ax sr else use -ax map-ont
    
    """
    
    minimap2 \\
        -N 1 $args $mapx \\
        -t $cpu_view -I $I_value \\
        $reference \\
        $input_reads \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}