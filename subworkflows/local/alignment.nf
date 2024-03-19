//
// Check input samplesheet and get read channels
//

include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
// include { BAM_TO_SAM } from "../../modules/local/bam_to_sam"
include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COVERAGE } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_HIST_COVERAGE  }  from '../../modules/local/samtools_hist_coverage'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_ILLUMINA } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_INDEX  } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_OXFORD } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'
include { MERGE_FASTA } from '../../modules//local/merge_fasta'
include { SPLIT_VCF } from '../../modules/local/split_vcf'
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { RSEQC_BAMSTAT } from '../../modules/nf-core/rseqc/bamstat/main'
include { BEDTOOLS_DEPTHCOVERAGE } from '../../modules/local/bedtools_coverage'

workflow ALIGNMENT {
    take:
    fastq_reads

    main:
    ch_bams = Channel.empty()
    ch_fasta = Channel.empty()
    ch_pileups = Channel.empty()
    ch_stats = Channel.empty()
    ch_depths = Channel.empty()
    ch_pafs = Channel.empty()
    ch_sams = Channel.empty()
    ch_aligners = Channel.empty()
    collected_bams  = Channel.empty()
    ch_bamstats = Channel.empty()
    ch_merged_fasta = Channel.empty()
    ch_merged_mpileup = Channel.empty()


    fastq_reads
        .branch{
            longreads: it[0].platform =~ 'OXFORD'
            shortreads: it[0].platform =~ 'ILLUMINA'
    }.set { ch_aligners }


    ch_versions = 1


    if (params.bt2_indices){
        ch_indices = ch_aligners.shortreads.map{ m, fastq, fasta, bed, map -> return [ m, fastq, file(params.bt2_indices) ] }


        println("Using pre-built indices provided $params.bt2_indices")
    } else {
        println("Building indices from FASTA references")
        BOWTIE2_BUILD(
            ch_aligners.shortreads.map{  m, fastq, fasta, bed, map  ->
                return [ m, fasta  ]
            }

        )
        ch_aligners.shortreads.map{  m, fastq, fasta, bed, map ->
            return [ m, fastq ]
        }.join(BOWTIE2_BUILD.out.index).set { ch_indices  }
    }
    BOWTIE2_ALIGN(
        ch_indices,
        false,
        true
    )
    illumina_alignments = ch_aligners.shortreads.join(
        BOWTIE2_ALIGN.out.aligned
    )

    MINIMAP2_ALIGN(
        ch_aligners.longreads.map{ m, fastq, fasta, bed, map -> [ m, fastq, fasta ] },
        true,
        true,
        true
    )
    nanopore_alignments = ch_aligners.longreads.join(
        MINIMAP2_ALIGN.out.bam
    )

    nanopore_alignments.mix(illumina_alignments).map{
        m, fastq, fasta, bed, map, bam -> [ m, bam ]
    }.set { collected_bams }




    SAMTOOLS_INDEX(
        collected_bams
    )

    sorted_bams = collected_bams

    SAMTOOLS_DEPTH(
        sorted_bams
    )

    collected_bams.join(SAMTOOLS_INDEX.out.bai).set{ sorted_bams_with_index }

    SAMTOOLS_COVERAGE(
        sorted_bams_with_index
    )
    gcf_with_bam = sorted_bams.join(fastq_reads.map{ m, fastq, fasta, bed, map -> return [m, map] })
    SAMTOOLS_HIST_COVERAGE (
        gcf_with_bam
    )

    sorted_bams_with_index.join(fastq_reads.map{ m, fastq, fasta, bed, map -> return [m, fasta] }).set{ sorted_bams_with_index_fasta }



    if (params.get_features){

        ch_merged_bed = sorted_bams.join(
            fastq_reads.map{ m, fastq, fasta, bed, map -> return [m, bed] }
        ).map{
            meta, bam, bed -> [meta, bed, bam]
        }.filter{
            it[1] != null
        }
        BEDTOOLS_DEPTHCOVERAGE(
            ch_merged_bed
        )
    }



    // branch out the samtools_sort output to nanopore and illumina
    ch_sorted_bam_split = sorted_bams.join(fastq_reads.map{ m, fastq, fasta, bed, map -> return [m, fasta] })

    ch_sorted_bam_split.branch{
            longreads: it[0].platform =~ 'OXFORD'
            shortreads: it[0].platform =~ 'ILLUMINA'
    }.set { ch_sorted_bam_split }


    if (params.get_variants || params.reference_assembly){
        BCFTOOLS_MPILEUP_OXFORD(
            ch_sorted_bam_split.longreads.map{ m, bam, fasta -> [m, bam] },
            ch_sorted_bam_split.longreads.map{ m, bam, fasta -> fasta },
            false
        )
        ch_merged_mpileup_oxford = BCFTOOLS_MPILEUP_OXFORD.out.vcf.join(BCFTOOLS_MPILEUP_OXFORD.out.tbi)
        ch_merged_mpileup_oxford = ch_merged_mpileup_oxford.join(ch_sorted_bam_split.longreads.map{ m, bam, fasta -> [m, fasta] })

        BCFTOOLS_MPILEUP_ILLUMINA(
            ch_sorted_bam_split.shortreads.map{ m, bam, fasta -> [m, bam] },
            ch_sorted_bam_split.shortreads.map{ m, bam, fasta -> fasta },
            false
        )
        .set{ transposed_fastas_oxford }

        ch_merged_mpileup_illumina = BCFTOOLS_MPILEUP_ILLUMINA.out.vcf.join(BCFTOOLS_MPILEUP_ILLUMINA.out.tbi)
        ch_merged_mpileup_illumina = ch_merged_mpileup_illumina.join(ch_sorted_bam_split.shortreads.map{ m, bam, fasta -> [m, fasta] })

        ch_merged_mpileup = ch_merged_mpileup_illumina.mix(ch_merged_mpileup_oxford)

        if (params.reference_assembly){
            BCFTOOLS_CONSENSUS(
                ch_merged_mpileup
            )
            ch_fasta = BCFTOOLS_CONSENSUS.out.fasta
        }
    }
    RSEQC_BAMSTAT(
        sorted_bams
    )
    ch_bamstats = RSEQC_BAMSTAT.out.txt
    ch_bams =  sorted_bams_with_index

    ch_depths = SAMTOOLS_DEPTH.out.tsv

    emit:
        // sam  = ch_sams // channel: [ val(meta), [ paffile ] ] ]
        depth = ch_depths
        mpileup = ch_merged_mpileup
        bams = ch_bams
        fasta  = ch_merged_fasta
        stats = ch_stats
        bamstats = ch_bamstats
        versions = ch_versions
    }
