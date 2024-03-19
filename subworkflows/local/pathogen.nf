//
// Cross-check alignment abundances with pathogen list
//
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


include { PATHOGENS_FIND_SAMPLE } from '../../modules/local/pathogens_find'
include { ORGANISM_MERGE_REPORT } from '../../modules/local/report_merge'

workflow PATHOGENS {
    take:
        alignments
        pathogens_list
    main:
        ch_pathogens_report = Channel.empty()
        if (!pathogens_list){
            println ("No pathogens list provided, skipping pathogen detection")
        } else{
            PATHOGENS_FIND_SAMPLE(
                alignments.combine(pathogens_list)
            )
            // collect all outputs FIND_PATHOGENS.out.txt into a single channel
            // and assign it to the variable pathogens_list
            ch_empty_file = file("$projectDir/assets/NO_FILE")

            if (!params.distributions){
                distributions = ch_empty_file
            } else{
                distributions = file(params.distributions)
            }
            full_list_pathogen_files = PATHOGENS_FIND_SAMPLE.out.txt.map{m, txt -> txt}.collect()
            ORGANISM_MERGE_REPORT(
                full_list_pathogen_files,
                distributions
            )
            ch_pathogens_report = ORGANISM_MERGE_REPORT.out.report
        }
    emit:
        pathogens_list = ch_pathogens_report
}
