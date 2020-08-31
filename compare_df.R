## Load packages
list.of.packages <- c("reticulate", "optparse", "ggplot2",
                      "reshape2", "optparse", "dada2", "ShortRead", "foreach", "iterators", "doParallel")
lapply(list.of.packages,library,character.only = TRUE)

##Specifying Illumina vs. PacBio files, and what the sample name is.
option_list <- list(make_option(c("-s", "--script_path"), type="character", default=NULL, help="Directory where scripts are located.", 
                                metavar="character"),
                    make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"),
                    make_option(c("-p", "--pacbio"), type="character", default=FALSE, help="Specify if these files are only PacBio.", 
                                metavar="character", action="store_true"),
                    make_option(c("-i", "--illumina"), type="character", default=FALSE, help="Specify if these files are only Illumina.", 
                                metavar="character", action="store_true"),
                    make_option(c("-c", "--cpus"), type="character", default=FALSE, help="task.cpus", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

cpus = opt$cpus
num_cores <- strsplit(cpus, "[ gb]")[[1]][1]
print(num_cores)
# Specify number of cores for parallelizing frequency table creation
registerDoParallel(cores=num_cores)

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
syph_path <- opt$script_path

noprimer_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-8),"noprimers.fastq",sep ='')
nop <- file.path(noprimer_filenames)

if(opt$illumina == FALSE) {
  RAD_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-8),"noprimers.filtered.RAD.fasta",sep ='')
  RAD_files <- file.path(RAD_filenames)
  RAD_files_nolines <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fasta",sep ='')
  RAD_files_fix <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fix.fasta",sep ='')
  
  print("Making PacBio frequency files...")
  # Calls syph_r to make the final_data.csv for each PacBio sample. Parallelized.
  foreach(i=1:length(RAD_files_fix)) %dopar% {
    print(RAD_files_fix)
    PB_file_name <- substr(RAD_files_fix[i],1,nchar(RAD_files_fix[i]))
    print(PB_file_name)
    syphrPacBio_command <- paste("python3 ",syph_path," -i fasta -pacbio -d . -s ",PB_file_name,sep='')
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
  Illumina_freq_path <- "."
  
  # Creates frequency files for Illumina (final_data.csvs). Parallelized.
  foreach(i=1:length(Illumina_fns)) %dopar% {
    Illumina_file_name <- substr(Illumina_fns[i],1,nchar(Illumina_fns[i])-3)
    print(Illumina_file_name)
    syphrIllumina_command <- paste("python3 ",syph_path," -i fastq -illumina -d . -s ",Illumina_file_name,sep='')
    system(syphrIllumina_command)
  }

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
    
    # If there are any PacBio files, then compare_pacbio_df.csv should have been made, and we merge with compare_illumina_df to obtain allreads.csv.
    # Otherwise, we copy compare_illumina_df for allreads.csv.
    if(file.size("compare_pacbio_df.csv")>0) {
        allreads <- merge(compare_Illumina_df,read.csv("compare_pacbio_df.csv", check.names=FALSE),all=T)
        print("Making allreads.csv...")
        allreads_out <- paste("allreads.csv",sep='')
        write.csv(allreads,file=allreads_out,row.names=FALSE,quote=FALSE)
    } else {
        print("compare_pacbio_df.csv is empty, using compare_illumina_df as allreads.csv...")
        allreads <- compare_Illumina_df
        system("cp compare_illumina_df.csv allreads.csv")
    }
  }
}