#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    tprK pipeline, with nextflow!

    An example command for running the pipeline is as follows:
    nextflow run michellejlin/tprk -r nextflow --INPUT ./ --OUTDIR output/ -resume -with-docker ubuntu:18.04 -with-trace \\

    Mandatory Arguments
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
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
        --LARGE         For very large datasets. Will not generate comparisons between files and will output htmls for variable regions separately.
        --REFERENCE     Specify Illumina sample name (not file), to compare others to for dot-line plots. Can be used in tandem with --LARGE.
    
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
params.LARGE = false
params.REFERENCE = false

INPUT_TYPE = "both"
PACBIO_FLAG = ""
ILLUMINA_FLAG = ""

SYPH_R = file("${baseDir}/syph_r.py")
COMPARE_DF = file("${baseDir}/compare_df.R")
FILTER_ALL_READS = file("${baseDir}/filterAllReads.py")
RECALCULATE_FREQUENCY = file("${baseDir}/recalculate_frequency.R")
SYPH_VISUALIZER = file("${baseDir}/syph_visualizer.py")
RAD_FREQUENCY = file("${baseDir}/RAD_Frequency.R")
PACBIO_VS_ILLUMINA = file("${baseDir}/PacBio_v_Illumina_plots.R")
VARIABLE_REGION_COMPARE = file("${baseDir}/Variable_region_compare.R")
ALLDATA_VISUALIZER = file("${baseDir}/alldata_visualizer_alex.py")
PACBIOTREE = file("${baseDir}/PacBio2tree.R")

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
        .map{ row-> tuple(row.SampleName, file(row.PacBio)) }
    metadata_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> tuple(row.SampleName, file(METADATA_FILE), file(row.PacBio), INPUT_TYPE) }
}

sample_name_ch = Channel
        .fromPath(METADATA_FILE)
        .splitCsv(header:true)
        .map{ row-> row.SampleName }


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*        CREATE FREQUENCY TABLES PER SAMPLE          */
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
        errorStrategy 'retry'
        maxRetries 3

        input: 
            tuple val(sample_name), file(PACBIO_FILE) from input_pacbio_ch
            file(RAD_FREQUENCY)
        output:
            tuple val(sample_name), file("PB_${sample_name}.noprimers.filtered.RAD.nolines.fix.fasta") into pacbio_ch
            tuple val(sample_name), file("PB_${sample_name}.noprimers.filtered.RAD.nolines.fix.fasta") into pacbio_ch2
            val(sample_name) into pacbio_sample_name_ch
        
        script:
        """
        gunzip -d --force *.gz
        Rscript ${RAD_FREQUENCY} -s ${baseDir} -d ${params.INPUT} -m ${METADATA_FILE} -a ${PACBIO_FILE} -c ${task.cpus}
        """
    }

    // Create frequency tables for each PacBio sample.
    process createFrequencyTables_PacBio {
    container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        input: 
            file(PACBIO_FILE) from pacbio_ch.collect()
            file(METADATA_FILE)
            file(COMPARE_DF)
            file(SYPH_R)
        output:
            file("*final_data.csv") into pacbio_final_data_ch
            file("*final_data.csv") into final_data_ch_pb

            file "all_assignments.csv" into all_assignments_ch1
            file "compare_pacbio_df.csv" into compare_pacbio_ch
        
        script:
        """
        gunzip -d --force *.gz
        Rscript ${COMPARE_DF} -s ${SYPH_R} -m ${METADATA_FILE} -d ./ --pacbio -c ${task.cpus}
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
    // Also grabs frequency tables from PacBio samples and merges the two,
    // creating allreads.csv file.
    process createFrequencyTables_Illumina {
    container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        input: 
            file(ILLUMINA_FILE) from illumina_ch.collect()
            file "all_assignments.csv" from all_assignments_ch1
            file("compare_pacbio_df.csv") from compare_pacbio_ch
            file(METADATA_FILE)
            file(COMPARE_DF)
            file(SYPH_R)
        output:
            file("*final_data.csv") into illumina_final_data_ch
            file("*final_data.csv") into final_data_ch_ill
            file("all_assignments.csv") into all_assignments_ch2
            file("allreads.csv") into allreads_ch
        
        script:
        """
        gunzip -d --force *.gz
        Rscript ${COMPARE_DF} -s ${SYPH_R} -m ${METADATA_FILE} -d ./ --illumina -c ${task.cpus}
        """
    }
}

//
// All reads section
// 

// Filters allreads.csv based on set parameters and 
// recalculates relative frequencies after filter.
process filterReads {
    container "quay.io/greninger-lab/tprk:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    publishDir params.OUTDIR, mode: 'copy'

    input:
      file "allreads.csv" from allreads_ch
      file METADATA_FILE
      file FILTER_ALL_READS
      file RECALCULATE_FREQUENCY
    output:
      file "allreads_filtered.csv" into allreads_filt_ch
      file "allreads_filtered_heatmap.csv" into allreads_filt_heatmap_ch
      file "allreads.csv" into allreads_ch2

    script:
    """
    # Creates allreads_filtered.csv and recalculates the relative frequencies.
    python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a "allreads.csv"
    Rscript ${RECALCULATE_FREQUENCY} -f "allreads_filtered.csv" -m ${METADATA_FILE}
    
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

if (INPUT_TYPE != "pacbio") {

    // Filters allreads.csv based on set parameters and 
    // recalculates relative frequencies after filter.
    // Create relative frequency plots for Illumina samples. Default html.
    process createFrequencyPlots_Illumina {
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}/Figures/Relative_Frequency_Plots", mode: 'copy', pattern: '*_RelativeFreqPlot*'

        input:
        file(FINAL_DATA) from final_data_ch_ill
        file(SYPH_VISUALIZER)
        file(FILTER_ALL_READS)
        val(sample_name) from sample_name_ch
        file(RECALCULATE_FREQUENCY)
        file(METADATA_FILE)

        output:
        tuple val(sample_name), file("Ill_${sample_name}_final_data_filtered.csv") into final_data_filtered_ch_ill
        file("*_RelativeFreqPlot*") into relative_freq_plot_ch_ill

        script:
        """
        python3 ${SYPH_VISUALIZER} Ill_${sample_name}_final_data.csv -t Ill_${sample_name} -o ./
        
        python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a Ill_${sample_name}_final_data.csv
        Rscript ${RECALCULATE_FREQUENCY} -f Ill_${sample_name}_final_data_filtered.csv -m ${METADATA_FILE}
        python3 ${SYPH_VISUALIZER} Ill_${sample_name}_final_data_filtered.csv -t Ill_${sample_name}_filtered -o ./

        """
        }
}

if (INPUT_TYPE != "illumina") {
    // Filters allreads.csv based on set parameters and 
    // recalculates relative frequencies after filter.
    // Create relative frequency plots for PacBio samples. Default html.
    process createFrequencyPlots_PacBio {
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}/Figures/Relative_Frequency_Plots", mode: 'copy', pattern: '*_RelativeFreqPlot*'

        input:
        file(FINAL_DATA) from final_data_ch_pb
        file(SYPH_VISUALIZER)
        file(FILTER_ALL_READS)
        file(METADATA_FILE)
        file(RECALCULATE_FREQUENCY)
        val(sample_name) from pacbio_sample_name_ch

        output:
        tuple val(sample_name), file("PB_${sample_name}.noprimers.filtered.RAD.nolines.fix_final_data_filtered.csv") into final_data_filtered_ch_pb
        file("*_RelativeFreqPlot*") into relative_freq_plot_pb


        script:
        """
        python3 ${SYPH_VISUALIZER} PB_${sample_name}.noprimers.filtered.RAD.nolines.fix_final_data.csv -t PB_${sample_name} -o ./
        
        python3 ${FILTER_ALL_READS} -f ${params.RF_FILTER} -c ${params.COUNT_FILTER} -a PB_${sample_name}.noprimers.filtered.RAD.nolines.fix_final_data.csv
        Rscript ${RECALCULATE_FREQUENCY} -f PB_${sample_name}.noprimers.filtered.RAD.nolines.fix_final_data_filtered.csv -m ${METADATA_FILE}
        python3 ${SYPH_VISUALIZER} PB_${sample_name}.noprimers.filtered.RAD.nolines.fix_final_data_filtered.csv -t PB_${sample_name}_filtered -o ./

        """
        }
}


// Generates PacBio vs. Illumina scatterplots for each sample. Compares filtered and non-filtered side by side,
// as well as automatically generates a zoomed in version from 0-10% relative frequency.
// Does not occur if running only PacBio or only Illumina files.
if (INPUT_TYPE == "both") {
    process createPacbioVsIlluminaPlots {
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}Figures/PacBio_vs_Illumina_Plots", mode: 'copy'

        input:
        file("allreads.csv") from allreads_ch
        file("allreads_filtered.csv") from allreads_filt_ch
        tuple val(sample_name), file(ILLUMINA_FILE), file(PACBIO_FILE), val(INPUT_TYPE) from metadata_ch
        file(PACBIO_VS_ILLUMINA)

        output:
        file("*.pdf") into pacbio_v_illumina_plots_ch
        file("*.RData") into pacbio_v_illumina_plots_ch2

        script:
        """
        Rscript ${PACBIO_VS_ILLUMINA} -p ./ -s ${sample_name}
        """
    }
}

// Generates dot-line plots for comparing variable regions between two samples.
// By default, does only filtered plots.
if (INPUT_TYPE != "pacbio") {
    process createVariableRegionComparisons {
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}Figures/Variable_Region_Comparisons", mode: 'copy'

        input:
        file("allreads.csv") from allreads_ch
        file("allreads_filtered.csv") from allreads_filt_ch
        file(VARIABLE_REGION_COMPARE)
        file(METADATA_FILE)

        output:
        file("*.pdf") into variable_region_ch
        file("*.RData") into variable_region_ch2

        script:
        if (params.REFERENCE != false) {
            """
            Rscript ${VARIABLE_REGION_COMPARE} -d ./ -m ${METADATA_FILE} -c ${task.cpus} -r ${params.REFERENCE}
            """
        }
        else if (params.LARGE == false) {
            """
            Rscript ${VARIABLE_REGION_COMPARE} -d ./ -m ${METADATA_FILE} -c ${task.cpus}
            """
        } else {
            """
            echo "--LARGE specified and no --REFERENCE given. Skipping making variable region comparisons between Illumina samples..."
            """
        }
        
    }
}


if (INPUT_TYPE != "illumina") {
    // Creates a ggtree of all the PacBio samples.
    // Currently automatically roots by midpoint, but will have to manually 
    // reorder by certain branch if needed.
    // Also saves ggtree in RDS which can be accessed with readRDS() in R for manual edits.
    process createPacBioTree{
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}Figures/Tree", mode: 'copy'

        input:
        file(METADATA_FILE)
        file(PACBIOTREE)
        file(PACBIO_FILE) from pacbio_ch2.collect()

        output:
        file("PacBio_Tree_Filtered.pdf") into tree_ch
        file("*_fullORFs.fasta") into tree_ch1
        file("*.tsv") into tree_ch2
        file("*.nwk") into tree_ch3
        file("*.RData") into tree_ch4

        script:
        """
        Rscript ${PACBIOTREE} -d . -m ${METADATA_FILE} -r ${params.RF_FILTER}
        """
    }
}

// Generates visualizations (heatmap and variable regions) for filtered
// allreads.csv. Default is html.
process visualizeAllData {
        container "quay.io/greninger-lab/tprk:latest"

        // Retry on fail at most three times 
        errorStrategy 'retry'
        maxRetries 3

        publishDir "${params.OUTDIR}Figures", mode: 'copy'

        input:
        file("allreads_filtered.csv") from allreads_filt_ch
        file(METADATA_FILE)
        file(ALLDATA_VISUALIZER)

        output:
        file("all_*") into alldata_visual_ch

        script:
        // Splits up htmls if large datasets, otherwise Bokeh can't open
        if (params.LARGE == false) {
            """
            python3 ${ALLDATA_VISUALIZER} allreads_filtered.csv ${METADATA_FILE}
            """
        } else {
            print("--LARGE flag specified. Splitting up alldata visualizations...")
            """
            python3 ${ALLDATA_VISUALIZER} allreads_filtered.csv ${METADATA_FILE} -large
            """
        }
}