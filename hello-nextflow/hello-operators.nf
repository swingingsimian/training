#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input (file of input files, one per line)
// params.reads_bam = "${projectDir}/data/sample_bams.txt"
params.reads_bam = "${projectDir}/data/sample_bams.local.txt"
// Base name for final output file
params.cohort_name = "family_trio"

// Output directory
params.outdir = "results_genomics"

// Accessory files
params.reference        = "${projectDir}/data/ref/ref.fasta"
params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
params.reference_dict   = "${projectDir}/data/ref/ref.dict"
params.intervals        = "${projectDir}/data/ref/intervals.bed"

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        tuple path(input_bam), path("${input_bam}.bai")

    script:
    """
    samtools index '$input_bam'
    """
}

/*
 * Call variants with GATK HaplotypeCaller
 */
process GATK_HAPLOTYPECALLER {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir params.outdir, mode: 'symlink'

    input:
        tuple path(input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list

    output:
        path "${input_bam}.g.vcf"     , emit: vcf
        path "${input_bam}.g.vcf.idx" , emit: idx

    // why are we having to replicate here, can we not use a variable for suffix to avoid mismatch? Overkill?

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.g.vcf \
        -L ${interval_list} \
        -ERC GVCF
    """
}


/*
 * Combine GVCFs into GenomicsDB datastore
 */
process GATK_JOINTGENOTYPING {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"
    publishDir params.outdir, mode: 'copy'

    input:
        path all_gvcfs
        path all_idxs
        path interval_list
        val cohort_name
        path ref_fasta
        // Following are used implicitly, but need explitly specifying so nextflow will stage this in the container
        path ref_index
        path ref_dict

    output:
        // path "${cohort_name}_gdb" // Old output for GenomicsDBImport only process
        path "${cohort_name}.joint.vcf"     , emit: vcf
        path "${cohort_name}.joint.vcf.idx" , emit: idx

    script:
    // """
    // gatk GenomicsDBImport \
    //     -V ${all_gvcfs} \ We can't specify multiple files here with just one -V flag
    //     -L ${interval_list} \
    //     --genomicsdb-workspace-path ${cohort_name}_gdb
    // """
    // This is how to define dynamic vars!
    // def gvcfs_line = all_gvcfs.collect { gvcf -> "-V ${gvcf}" }.join(' ')
    // but done inline for now

    """
    gatk GenomicsDBImport \
        ${all_gvcfs.collect { gcvf -> "-V ${gcvf}"}.join(' ')} \
        -L ${interval_list} \
        --genomicsdb-workspace-path ${cohort_name}_gdb

    gatk GenotypeGVCFs \
        -R ${ref_fasta} \
        -V gendb://${cohort_name}_gdb \
        -L ${interval_list} \
        -O ${cohort_name}.joint.vcf
    """
}



workflow {

    // Create input channel from a text file listing input file paths
    reads_ch = Channel.fromPath(params.reads_bam).splitText()

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file,
    )

    // Collect variant calling outputs across samples
    // out is actually a list of vcf idx pairs, this collects vcfs idxs across the list separately
    // This is a java/groovy thing, can we make it more explicit with out*? No :( ERROR ~ No such variable: vcf
    all_gvcfs_ch = GATK_HAPLOTYPECALLER.out.vcf.collect()
    all_idxs_ch = GATK_HAPLOTYPECALLER.out.idx.collect()

    all_gvcfs_ch.view()

    // // Combine GVCFs into a GenomicsDB datastore
    GATK_JOINTGENOTYPING(
        all_gvcfs_ch,
        all_idxs_ch,
        intervals_file,
        params.cohort_name,
        ref_file,
        ref_index_file,
        ref_dict_file,
    )

}
