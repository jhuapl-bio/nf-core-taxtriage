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
process DOWNLOAD_DB {

    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:4.2.0' :
        'gcc:bullseye' }"

    input:
    val(db)
    val(url)
    val(db_dir)
    val(checksum)

    when:
    task.ext.when == null || task.ext.when


    script: // This script is bundled with the pipeline, in nf-core/taxtriage/bin/

    """
    sh $PWD/bin/download_db.sh $db $PWD/$db_dir $url $checksum
    """
}
