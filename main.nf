#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    tprK pipeline, with nextflow!

    Usage: 

    An example command for running the pipeline is as follows:
    nextflow run michellejlin/tprk_nextflow \\

    Mandatory Arguments
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq.
        --OUTDIR        Output directßory.
        --METADATA      Metadata file formatted in a .csv with columns: SampleName, Illumina, PacBio.
                        If running with --PACBIO or --ILLUMINA simply leave those columns blank (but make sure to
                        have commas as appropriate). Names should be in PB_<sample_name>.fastq or Ill_<sample_name>.fastq format.
                        See the example metadatas in the example/ folder for a sample.
    
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

SYPH_R = file("${baseDir}/syph_r.py")
COMPARE_DF = file("${baseDir}/compare_df.R")
FILTER_ALL_READS = file("${baseDir}/filterAllReads.py")
RECALCULATE_FREQUENCY = file("${baseDir}/recalculate_frequency.R")
SYPH_VISUALIZER = file("${baseDir}/syph_visualizer.py")
RAD_FREQUENCY = file("${baseDir}/RAD_Frequency.R")

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
    input_pacbio_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.PacBio)) }
    illumina_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina)) }
    metadata_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina), file(row.PacBio), INPUT_TYPE) }
} else if (INPUT_TYPE == "illumina") {
    illumina_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina)) }
    metadata_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(row.Illumina), file(METADATA_FILE), INPUT_TYPE) }
} else if (INPUT_TYPE == "pacbio") {
    input_pacbio_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .ifEmpty()
        .map{ row-> tuple(row.SampleName, file(row.PacBio)) }
    metadata_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(METADATA_FILE), file(row.PacBio), INPUT_TYPE) }
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*        CREATE FREQUENCY PLOTS PER SAMPLE           */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


// 
// PacBio Section
// 
if(INPUT_TYPE != "illumina") {
    // Denoises PacBio files with RAD.
    process denoisePacBioFiles {
        container "quay.io/greninger-lab/tprk"

        // Retry on fail at most three times 
        // errorStrategy 'retry'
        // maxRetries 3

        input: 
            tuple val(sample_name), file(PACBIO_FILE) from input_pacbio_ch
            file(RAD_FREQUENCY)
        output:
            tuple val(sample_name), file("PB_${sample_name}.noprimers.filtered.RAD.nolines.fix.fasta") into pacbio_ch
        
        script:
        """
        Rscript ${RAD_FREQUENCY} -s ${baseDir} -d ${params.INPUT} -m ${METADATA_FILE} -a ${PACBIO_FILE}
        """
    }

    // Create frequency tables for each PacBio sample.
    process createFrequencyTables_PacBio {
    container "quay.io/greninger-lab/tprk"

        // Retry on fail at most three times 
        // errorStrategy 'retry'
        // maxRetries 3

        input: 
            file(PACBIO_FILE) from pacbio_ch.collect()
            file(METADATA_FILE)
            file(COMPARE_DF)
        output:
            file("*final_data.csv") into pacbio_final_data_ch
            file("*final_data.csv") into final_data_ch_pb

            file "all_assignments.csv" into all_assignments_ch1
            file "compare_pacbio_df.csv" into compare_pacbio_ch
        
        script:
        """
        Rscript ${COMPARE_DF} -s ${baseDir} -m ${METADATA_FILE} -d ./ --pacbio
        """
    }
}



//
// Illumina Section
// 

// If --ILLUMINA, no all_assignments.csv will have been created in
// createFrequencyPlots_PacBio. Creates it here.
if (INPUT_TYPE == "illumina") {
    process createAllAssignments{
        input:
        output:
            file("all_assignments.csv") into all_assignments_ch1
            file("compare_pacbio_df.csv") into compare_pacbio_ch
        
        script:
        """
        touch all_assignments.csv
        touch compare_pacbio_df.csv
        """
    }
}

if (INPUT_TYPE != "pacbio") {
    // Create frequency tables for each Illumina sample.
    // Also grabs frequency tables from PacBio samples and merged the two,
    // creating allreads.csv file.
    process createFrequencyTables_Illumina {
    container "quay.io/greninger-lab/tprk"

        // Retry on fail at most three times 
        // errorStrategy 'retry'
        // maxRetries 3

        input: 
            file(ILLUMINA_FILE) from illumina_ch.collect()
            file "all_assignments.csv" from all_assignments_ch1
            file("compare_pacbio_df.csv") from compare_pacbio_ch
            file(METADATA_FILE)
            file(COMPARE_DF)
        output:
            file("*final_data.csv") into illumina_final_data_ch
            file("*final_data.csv") into final_data_ch_ill
            file("all_assignments.csv") into all_assignments_ch2
            file("allreads.csv") into allreads_ch
        
        script:
        """
        Rscript ${COMPARE_DF} -s ${baseDir} -m ${METADATA_FILE} -d ./ --illumina
        """
    }
}

// Filters allreads.csv based on set parameters and 
// recalculates relative frequencies after filter.
process filterReads {
    container "quay.io/greninger-lab/tprk"

	// Retry on fail at most three times 
    // errorStrategy 'retry'
    // maxRetries 3

    input:
      file "allreads.csv" from allreads_ch
      file METADATA_FILE
      file FILTER_ALL_READS
      file RECALCULATE_FREQUENCY
    output:
      file "allreads_filtered.csv" into allreads_filt_ch
      file "allreads_filtered_heatmap.csv" into allreads_filt_heatmap_ch

    script:
    """
    # Creates allreads_filtered.csv and recalculates the relative frequencies.
    python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a "allreads.csv"
    Rscript ${RECALCULATE_FREQUENCY} -f "allreads.csv" -m ${METADATA_FILE}
    
    # Creates allreads_filtered_heatmap.csv and recalculates the relative frequencies. This csv includes samples under the count/relative freq filters
	# if another sample shares the same read. 
    python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a "allreads.csv" -is_heatmap
    Rscript ${RECALCULATE_FREQUENCY} -f "allreads_filtered_heatmap.csv" -m ${METADATA_FILE}

    """
}


// ////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////
// /*                                                    */
// /*                CREATE VISUALIZATIONS               */
// /*                                                    */
// ////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////

// Filters allreads.csv based on set parameters and 
// recalculates relative frequencies after filter.
process createFrequencyPlots_Illumina {
    container "quay.io/greninger-lab/tprk"

	// Retry on fail at most three times 
    // errorStrategy 'retry'
    // maxRetries 3

    input:
      file(FINAL_DATA) from final_data_ch_ill
      file(SYPH_VISUALIZER)
      file(FILTER_ALL_READS)

    output:
      file("*_filtered.csv") into final_data_filtered_ch_ill
      file("*_RelativeFreqPlot*") into relative_freq_plot_ch_ill

    script:
    """
    base=`basename ${FINAL_DATA} "_final_data.csv"`
    python3 ${SYPH_VISUALIZER} ${FINAL_DATA} -t \${base} -o ./
    
    python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a ${FINAL_DATA}
    """
    }

process createFrequencyPlots_PacBio {
    container "quay.io/greninger-lab/tprk"

	// Retry on fail at most three times 
    // errorStrategy 'retry'
    // maxRetries 3

    input:
      file(FINAL_DATA) from final_data_ch_pb
      file(SYPH_VISUALIZER)
      file(FILTER_ALL_READS)

    output:
      file("*_filtered.csv") into final_data_filtered_ch_pb
      file("*_RelativeFreqPlot*") into relative_freq_plot_pb

    script:
    """
    base=`basename ${FINAL_DATA} "_final_data.csv"`
    python3 ${SYPH_VISUALIZER} ${FINAL_DATA} -t \${base} -o ./
    
    python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a ${FINAL_DATA}
    """
    }