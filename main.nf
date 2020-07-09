#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    tprK pipeline, with nextflow!

    Usage: 

    An example command for running the pipeline is as follows:
    nextflow run vpeddu/lava \\

    Mandatory Arguments
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq.
        --OUTDIR        Output directory.
        --METADATA      Metadata file formatted in a .csv with columns: SampleName, Illumina, PacBio.
                        If running with --PACBIO or --ILLUMINA simply leave those columns blank (but make sure to
                        have commas as appropriate).
    
    Input Specifications
        --PACBIO        Write this flag to specify that there are only PacBio files here.
                        Comparison figures to Illumina will not be generated.
        --ILLUMINA      Write this flag to specify that there are only Illumina files here.
                        Comparison figures to PacBio will not be generated.
    
    Filtering Options
        --RF_FILTER         Optional flag for specifying what relative frequency
                            an additional filtered final merged table and visualizations should be sorted at.
                            By default this is set to 0.2.
        --COUNT_FILTER      Optional flag for specifying what count
                            an additional filtered final merged table and visualizations should be sorted at.
                            By default this is set to 5.
        --ILLUMINA_FILTER   Optional flag for specifying if PacBio reads should only include
                            Illumina-supported reads that pass the filters given. Default relative freq is
                            set to 0.2 and count is set to 5.
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.INPUT = false
params.OUTDIR= false
params.RF_FILTER = 0.2
params.COUNT_FILTER = 5
params.ILLUMINA_FILTER = false
params.PACBIO = false
params.ILLUMINA = false
params.METADATA = false

INPUT_TYPE = "both"
PACBIO_FLAG = ""
ILLUMINA_FLAG = ""

/////////////////////////////
/*    VALIDATE INPUTS      */
/////////////////////////////

// if METADATA not set
if (params.METADATA == false) {
    println("Must provide metadata file input as .csv format with three columns: \
    SampleName, Illumina, PacBio. Use --METADATA flag.")
    exit(1)
} else{
    METADATA_FILE = file(params.METADATA)
}
// if INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
   params.INPUT = "${params.INPUT}/"
}
// if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}
// Figure out input of PACBIO or ILLUMINA
if(((params.PACBIO) && (params.ILLUMINA))){ 
    println("--PACBIO and --ILLUMINA cannot be used together. Please specify only one, \
    or do not use these flags if you want to run the default way of comparing both PacBio and Illumina files.")
    exit(1)
} else if (params.PACBIO) {
    INPUT_TYPE = "pacbio"
    PACBIO_FLAG = "--pacbio"
    println("--PACBIO indicated. Will not compare to Illumina files or generate figures.")
} else if (params.ILLUMINA) {
    INPUT_TYPE = "illumina"
    ILLUMINA_FLAG = "--illumina"
    println("--ILLUMINA indicated. Will not compare to PacBio files or generate figures.")
}

// Helpful messages
println("Will filter final products for >${params.RF_FILTER} relative frequency \
and >${params.COUNT_FILTER} count.")
if (params.ILLUMINA_FILTER){
    println("Will only include PacBio reads supported by Illumina reads that pass the filter, \
    as specified by --ILLUMINA_FILTER.")
}

/////////////////////////////
/*    METADATA PARSING     */
/////////////////////////////

// Reads in Illumina/PacBio pairs from metadata
if (INPUT_TYPE == "both") {
    input_read_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina), file(row.PacBio), INPUT_TYPE) }
} else if (INPUT_TYPE == "illumina") {
    input_read_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina), file(METADATA_FILE), INPUT_TYPE) }
} else if (INPUT_TYPE == "pacbio") {
    input_read_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(METADATA_FILE), file(row.PacBio), INPUT_TYPE) }
}



////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                CREATE ALLREADS FILE                */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


// Goes from the original Illumina and PacBio reads to the all-important allreads.csv.
// Along the way, also makes frequency tables for each sample.
process createAllReads { 
    container "quay.io/greninger-lab/tprk"

	// Retry on fail at most three times 
    // errorStrategy 'retry'
    // maxRetries 3

    input:
      file METADATA_FILE
    output: 

    script:
    if (INPUT_TYPE == "both") {
        """
        Rscript ${baseDir}/og_files_to_all_reads.R -s ${baseDir} -d ${params.INPUT}
        """
    } else if (INPUT_TYPE == "illumina") {
        """
        Rscript ${baseDir}/og_files_to_all_reads.R -s ${baseDir} -d ${params.INPUT} --illumina
        """
    } else if (INPUT_TYPE == "pacbio") { 
        """
        echo ${METADATA_FILE}
        Rscript ${baseDir}/og_files_to_all_reads.R -s ${baseDir} -d ${params.INPUT} -m ${METADATA_FILE} --pacbio
        """
    }
}

// process Aligning {
//      container "quay.io/biocontainers/bbmap:38.86--h1296035_0"
//     //container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7	"

//     // Retry on fail at most three times 
//     errorStrategy 'retry'
//     maxRetries 3

//     input: 
//       tuple val(base), file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz") from Trim_out_ch
//       file REFERENCE_FASTA
//     output:
//       tuple val (base), file("${base}.bam") into Aligned_bam_ch

//     cpus 4 
//     memory '6 GB'

//     script:
//     """
//     #!/bin/bash

//     /usr/local/bin/bbmap.sh in1=${base}.R1.paired.fastq.gz in2=${base}.R2.paired.fastq.gz outm=${base}.bam ref=${REFERENCE_FASTA} -Xmx6g

//     """

//     // bwa?
//     // script:
//     // """
//     // #!/bin/bash
//     // bwa mem $workflow.projectDir/NC_045512.2.fasta ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz > aln.sam
//     // """
// }

// process NameSorting { 
//     container "quay.io/biocontainers/samtools:1.3--h0592bc0_3"

// 	// Retry on fail at most three times 
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//       tuple val (base), file("${base}.bam") from Aligned_bam_ch
//     output:
//       tuple val (base), file("${base}.sorted.sam") into Sorted_sam_ch

//     script:
//     """
//     #!/bin/bash
//     samtools sort -n -O sam ${base}.bam > ${base}.sorted.sam

//     """
// }

// process Clipping { 
//     container "quay.io/greninger-lab/swift-pipeline"

// 	// Retry on fail at most three times 
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//       tuple val (base), file("${base}.sorted.sam") from Sorted_sam_ch
//       file MASTERFILE
//     output:
//       tuple val (base), file("${base}.clipped.bam") into Clipped_bam_ch

//     script:
//     """
//     #!/bin/bash
//     /./root/.local/bin/primerclip ${MASTERFILE} ${base}.sorted.sam ${base}.clipped.sam
//     #/usr/local/miniconda/bin/samtools sort -n -O sam ${base}.clipped.sam > ${base}.clipped.sorted.sam
//     #/usr/local/miniconda/bin/samtools view -Sb ${base}.clipped.sorted.sam > ${base}.clipped.unsorted.bam
//     #/usr/local/miniconda/bin/samtools sort -o ${base}.clipped.unsorted.bam ${base}.clipped.bam
//      /usr/local/miniconda/bin/samtools sort ${base}.clipped.sam -o ${base}.clipped.bam

//     """
// }

// process generateConsensus {
//     container "quay.io/greninger-lab/swift-pipeline"

// 	// Retry on fail at most three times 
//     errorStrategy 'retry'
//     maxRetries 3

//     input:
//         tuple val (base), file(BAMFILE) from Clipped_bam_ch
//         file REFERENCE_FASTA
//     output:
//         file("${base}.fasta")
//         file("${base}.clipped.bam")

//     publishDir params.OUTDIR, mode: 'copy'

//     shell:
//     '''
//     #!/bin/bash
//     /usr/local/miniconda/bin/bcftools mpileup \\
//         --count-orphans \\
//         --no-BAQ \\
//         --max-depth 500000 \\
//         --fasta-ref !{REFERENCE_FASTA} \\
//         --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
//         !{BAMFILE} \\
//         | /usr/local/miniconda/bin/bcftools call --output-type v --ploidy 1 --keep-alts --keep-masked-ref --multiallelic-caller --variants-only -P 0 \\
//         | /usr/local/miniconda/bin/bcftools reheader --samples sample_name.list \\
//         | /usr/local/miniconda/bin/bcftools view --output-file !{base}.vcf.gz --output-type z


//     /usr/local/miniconda/bin/tabix -p vcf -f !{base}.vcf.gz


//     cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus !{base}.vcf.gz > !{base}.consensus.fa
//     cat !{REFERENCE_FASTA} !{base}.consensus.fa > align_input.fasta
//     /usr/local/miniconda/bin/mafft --auto align_input.fasta > repositioned.fasta
//     awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
    
//     python3 !{baseDir}/trim_ends.py !{base}



//     '''
// }
//     // # This approach gives a fasta with all N's (work/84/7e5c92/), vcf with lines starting from 20000s
//     // #/usr/local/miniconda/bin/samtools mpileup --max-depth 500000 -uf !{REFERENCE_FASTA} !{base}.clipped.bam | \
//     // #/usr/local/miniconda/bin/bcftools call -c -o !{base}.consensus.vcf
//     // #/usr/local/miniconda/bin/vcfutils.pl vcf2fq !{base}.consensus.vcf > !{base}.consensus.fastq 
//     // #/usr/local/miniconda/bin/seqtk seq -aQ64 -q20 -n N !{base}.consensus.fastq > !{base}.consensus.fasta

//     // # This approach gives a fasta identical to ref, blank vcf
//     // #/usr/local/miniconda/bin/bcftools mpileup -Ou --max-depth 500000 -f !{REFERENCE_FASTA} !{BAMFILE} | \
//     // #/usr/local/miniconda/bin/bcftools call -c -o !{base}.vcf
//     // ##/usr/local/miniconda/bin/bcftools call -mv -Oz -o !{base}.vcf.gz
//     // #/usr/local/miniconda/bin/bgzip !{base}.vcf
//     // #/usr/local/miniconda/bin/bcftools index !{base}.vcf.gz
//     // #cat !{REFERENCE_FASTA} | /usr/local/miniconda/bin/bcftools consensus !{base}.vcf.gz > !{base}.consensus.fasta

//     // This approach gave a bunch of ns in between
//     // /usr/local/miniconda/bin/samtools mpileup -uf !{REFERENCE_FASTA} !{BAMFILE} | /usr/local/miniconda/bin/bcftools call -c | /usr/local/miniconda/bin/vcfutils.pl vcf2fq > out.fastq
//     // /usr/local/miniconda/bin/samtools mpileup -uf !{REFERENCE_FASTA} !{BAMFILE} | /usr/local/miniconda/bin/bcftools call -c | /usr/local/miniconda/bin/bcftools view 

//     // /usr/local/miniconda/bin/seqtk seq -aQ64 -q20 -n N out.fastq > !{base}.consensus.fasta
    
