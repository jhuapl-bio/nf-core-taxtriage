// ##############################################################################################
// # Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
// # All rights reserved.
// # Permission is hereby granted, free of charge, to any person obtaining a copy of this
// # software and associated documentation files (the "Software"), to deal in the Software
// # without restriction, including without limitation the rights to use, copy, modify,
// # merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// # permit persons to whom the Software is furnished to do so.
// #
// # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// # INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// # PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// # LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
// # OR OTHER DEALINGS IN THE SOFTWARE.
// #
process PATHOGENS_FIND_SAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py39hcada746_1' :
        'biocontainers/pysam:0.21.0--py39hcada746_1' }"

    input:
    tuple val(meta), path(bamfiles), path(bai), path(mapping), path(pathogens_list)

    output:
        path "versions.yml"           , emit: versions
        tuple val(meta), path("*.txt")    , optional:false, emit: txt


    when:
    task.ext.when == null || task.ext.when




    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/

    def output = "${meta.id}.paths.txt"
    def id = meta.id
    def type = meta.type ? " -t ${meta.type} " : " -t Unknown "
    def min_reads_align = params.min_reads_align  ? " -r ${params.min_reads_align} " : " -r 3 "

    """

    match_paths.py \\
        -i $bamfiles \\
        -o $output \\
        -s $id \\
        $type \\
        $min_reads_align \\
        -p $pathogens_list \\
        -m $mapping

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version)
    END_VERSIONS

    """
}


