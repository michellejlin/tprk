## Load packages
list.of.packages <- c("reticulate", "optparse", "ggplot2",
                      "reshape2", "optparse", "dada2", "ShortRead")
lapply(list.of.packages,library,character.only = TRUE)

##Specifying Illumina vs. PacBio files, and what the sample name is.
option_list <- list(make_option(c("-s", "--script_path"), type="character", default=NULL, help="Directory where scripts are located.", 
                                metavar="character"),
                    make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"),
                    make_option(c("-p", "--pacbio"), type="character", default=FALSE, help="Specify if these files are only PacBio.", 
                                metavar="character", action="store_true"),
                    make_option(c("-i", "--illumina"), type="character", default=FALSE, help="Specify if these files are only Illumina.", 
                                metavar="character", action="store_true"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- "./"
script.dir <- opt$script_path

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are. (This is the -directory option).
## script.dir refers to the folder where all the script files are located. This should point to where you saved the cloned GitHub.

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files_redo3"
#script.dir <- "/Users/uwvirongs/Documents/tprK-master/"

## This script can also be run from the command line.
## Usage: rscript \path\to\og_files_to_all_reads.R -s [script_path] -d [directory]

#####

metadata <- read.table(opt$metadata, sep=',', header=TRUE)
PacBio_fns <- c(as.character(metadata$PacBio))
Illumina_fns <- c(as.character(metadata$Illumina))
sample_names <- c(as.character(metadata$SampleName))
syph_path <- paste(script.dir,"/syph_r.py",sep='')

noprimer_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.fastq",sep ='')
nop <- file.path(noprimer_filenames)

if(opt$illumina == FALSE) {
  RAD_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.RAD.fasta",sep ='')
  RAD_files <- file.path(RAD_filenames)
  RAD_files_nolines <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fasta",sep ='')
  RAD_files_fix <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fix.fasta",sep ='')
  
  print("Making PacBio frequency files...")
  # Calls syph_r to make the final_data.csv for each PacBio sample.
  # TODO: right now this doesn't work if it's a repeat run because of previous .fastq files. Fix so you can rerun safely
  # without having to delete stuff when rerunning.
  for (num in c(1:length(RAD_files_fix))) {
    syphrPacBio_command <- paste("python3 ",syph_path," -i fasta -pacbio -d ./",sep='')
    system(syphrPacBio_command)
  }
  
  ## Make PacBio frequency comparison file
  print("Making PacBio comparison dataframe...")
  PacBio_freq_files = list.files("./",pattern="*_final_data.csv")
  PacBio_freq_files_fullpath = list.files("./",pattern="*_final_data.csv",full.names=T)
  compare_PacBio_df <- data.frame(Region=character(),Read=character())
  
  sample_names <- sort(sample_names)
  
  # Renames the columns attaching PB and sample name to relative frequency and count so we can 
  # mash everything together into a big dataframe.
  for (i in 1:length(PacBio_freq_files)) {
    PacBioFreqtitle <- paste("PB_",sample_names[i],"_RelativeFreq",sep='')
    PacBioCounttitle <- paste("PB_",sample_names[i],"_Count",sep='')
    pacbiodf <- read.csv(PacBio_freq_files_fullpath[i],col.names = c("Region","Read",PacBioFreqtitle,PacBioCounttitle),check.names = FALSE)
    pacbiodf <- pacbiodf[order(pacbiodf$Region,-pacbiodf[[PacBioFreqtitle]]),][]    
    compare_PacBio_df <- merge(compare_PacBio_df,pacbiodf,all=TRUE)
    compare_PacBio_df_out <- paste("compare_pacbio_df.csv", sep='')
    write.csv(compare_PacBio_df,file=compare_PacBio_df_out,row.names=FALSE,quote=FALSE)
    # assign(PacBio_freq_files[i], read.csv(PacBio_freq_files_fullpath[i],col.names = c("Region","Read",PacBioFreqtitle,PacBioCounttitle)))
    # compare_PacBio_df <- merge(compare_PacBio_df,get(PacBio_freq_files[i]),all=TRUE)
  } 
} else {
  print("Illumina option specified. Skipping making PacBio frequency files...")
}

if (opt$pacbio == FALSE) {
  print("Making Illumina frequency files...")
  Illumina_freq_path <- "./"

  # Creates frequency files for Illumina (final_data.csvs).
  syphrIllumina_command <- paste("python3 ",syph_path," -i fastq -illumina -d ",Illumina_freq_path,"/ ;",sep='')
  system(syphrIllumina_command)

  print("Making Illumina comparison dataframe...")

  # Makes Illumina frequency comparison file
  Illumina_freq_files = list.files(Illumina_freq_path,pattern="*_final_data.csv")
  Illumina_freq_files_fullpath = list.files(Illumina_freq_path,pattern="*_final_data.csv",full.names=T)
  compare_Illumina_df <- data.frame(Region=character(),Read=character())
  #svgIllumina_command <- paste("/usr/bin/python /Users/uwvirongs/Documents/Michelle/syphilis/syph_visualizer_single_svg.py ",Illumina_freq_files_fullpath,sep='')
  #for (num in length(svgIllumina_command)) system(svgIllumina_command[num])
  
  # Renames columns to attach Ill_ and sample name onto relative frequency and count columns so we can mash everything
  # together into a big dataframe.
  for (i in 1:length(Illumina_freq_files)) {
    IlluminaFreqtitle <- paste("Ill_",sample_names[i],"_RelativeFreq",sep='')
    IlluminaCounttitle <- paste("Ill_",sample_names[i],"_Count",sep='')
    Illuminadf <- read.csv(Illumina_freq_files_fullpath[i],col.names = c("Region","Read",IlluminaFreqtitle,IlluminaCounttitle),check.names = FALSE)
    Illuminadf <- Illuminadf[order(Illuminadf$Region,-Illuminadf[[IlluminaFreqtitle]]),][]
    compare_Illumina_df <- merge(compare_Illumina_df,Illuminadf,all=TRUE)
    compare_Illumina_df_out <- paste("compare_illumina_df.csv", sep='')
    write.csv(compare_Illumina_df,file=compare_Illumina_df_out,row.names=FALSE,quote=FALSE)
    
    if(file.size("compare_pacbio_df.csv")>0) {
        print("compare_pacbio_df.csv is empty, using compare_illumina_df as allreads.csv...")
        allreads <- merge(compare_Illumina_df,read.csv("compare_pacbio_df.csv"),all=T)
    } else {
        allreads <- compare_Illumina_df
    }
    print("Making allreads.csv...")
    allreads_out <- paste("allreads.csv",sep='')
    write.csv(allreads,file=allreads_out,row.names=FALSE,quote=FALSE)
    } 
  } else {
  print("Pacbio option specified. Skipping making Illumina frequency comparison files.")
}