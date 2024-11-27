#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
// params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
// Primary input (array of three samples)
// params.reads_bam = [
//     "${projectDir}/data/bam/reads_mother.bam",
//     "${projectDir}/data/bam/reads_father.bam",
//     "${projectDir}/data/bam/reads_son.bam"
// ]

// Primary input (file of input files, one per line)
params.reads_bam = "${projectDir}/data/sample_bams.local.txt"

params.outdir    = "results_genomics"

//${projectDir} is a built-in Nextflow variable that points to the directory where the current Nextflow workflow script (hello-genomics.nf) is located.

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

    publishDir params.outdir, mode: 'symlink' // symlink is default and used for convenience here
    // You shouldn't do this in your final workflows, since you'll lose results when you clean up your work directory.

     input:
        path input_bam

    output:
        // path "${input_bam}.bai"
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
        // path input_bam
        // path input_bam_index
        tuple path(input_bam), path(input_bam_index)
        path ref_fasta
        path ref_index
        path ref_dict
        path interval_list
        // Note: ref_index/dict and input_bam_index are not used but still specified here?
        // We need to tell Nextflow explicitly that it has to stage those files in the working directory at 
        // runtime; otherwise it won't do it, and GATK will (correctly) throw an error about the index files 
        // being missing.

    output:
        path "${input_bam}.vcf"     , emit: vcf
        path "${input_bam}.vcf.idx" , emit: idx

    script:
    """
    gatk HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${input_bam}.vcf \
        -L ${interval_list}
    """
}



workflow {

    // // Create input channel (single file via CLI parameter)
    // reads_ch = Channel.fromPath(params.reads_bam)

    // Create input channel from a text file listing input file paths
    reads_ch = Channel.fromPath(params.reads_bam).splitText()

    // Create index file for input BAM file
    SAMTOOLS_INDEX(reads_ch)

    // temporary diagnostics
    reads_ch.view()
    SAMTOOLS_INDEX.out.view()

    // While main data inputs are streamed dynamically through channels, there are two approaches for handling 
    // accessory files. The recommended approach is to create explicit channels, which makes data flow clearer 
    // and more consistent. Alternatively, the file() function to create variables can be used for simpler cases, 
    // particularly when you need to reference the same file in multiple processes - though be aware this still 
    // creates channels implicitly.

    // Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

    // Call variants from the indexed BAM file
    GATK_HAPLOTYPECALLER(
        // These first two may not be in the same order and end up with mismatched bam and bai files
        // solution is to couple the files through the channels rather than using separate channels here
        // reads_ch,
        SAMTOOLS_INDEX.out,
        ref_file,
        ref_index_file,
        ref_dict_file,
        intervals_file
    )
    // nextflow run hello-genomics.nf -resume

}


// Why has samtools suddely started failing?
// It must have been run successfully on this machine previously
// The platform WARNING is a red herring here?
// The real issue is the gitpod paths? 

// ERROR ~ Error executing process > 'SAMTOOLS_INDEX (1)'

// Caused by:
//   Process `SAMTOOLS_INDEX (1)` terminated with an error exit status (125)


// Command executed:

//   samtools index 'reads_mother.bam'

// Command exit status:
//   125

// Command output:
//   (empty)

// Command error:
//   WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested
//   docker: Error response from daemon: Mounts denied: 
//   The path /workspace/gitpod/hello-nextflow/data/bam is not shared from the host and is not known to Docker.
//   You can configure shared paths from Docker -> Preferences... -> Resources -> File Sharing.
//   See https://docs.docker.com/desktop/settings/mac/#file-sharing for more info.

// Work dir:
//   /Users/nathanjohnson/src/nextflow-io/training/hello-nextflow/work/42/c2cc2a2b039eb4ceac452808861ebd

// Container:
//   community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464

// Tip: when you have fixed the problem you can continue the execution adding the option `-resume` to the run command line

//  -- Check '.nextflow.log' file for details