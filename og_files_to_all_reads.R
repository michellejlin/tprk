##Install packages if haven't already
list.of.packages <- c("JuliaCall", "reticulate", "devtools", "optparse", "devtools", "dada2", "ggplot2", "ShortRead",
                      "reshape2", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install(c("dada2"))
lapply(list.of.packages,library,character.only=T)


##Specifying Illumina vs. PacBio files, and what the sample name is.
option_list <- list(make_option(c("-s", "--script_path"), type="character", default=NULL, help="Directory where scripts are located.", 
                                metavar="character"),
                    make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-p", "--pacbio"), type="character", default=FALSE, help="Specify if these files are only PacBio.", 
                                metavar="character", action="store_true"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$directory
script.dir <- opt$script_path

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.

## path refers to the folder your metadata.csv and sequencing files (.fastq) are.
#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files"
## script.dir refers to the folder where all the script files are located. This should point to where you saved the cloned GitHub.
#script.dir <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline"

metadata <- read.table(paste(path,"/metadata.csv", sep=''), sep=',', header=TRUE)
setwd(path)
PacBio_fns <- c(as.character(metadata$PacBio))
Illumina_fns <- c(as.character(metadata$Illumina))
sample_names <- c(as.character(metadata$SampleName))

## Identify primers to go from ATG to stop in tprK
tprKF <- "GGAAAGAAAAGAACCATACATCC"
tprKR <- "CGCAGTTCCGGATTCTGA"
rc <- dada2:::rc
noprimer_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.fastq",sep ='')
nop <- file.path(path, noprimer_filenames)

## Remove primers
for (count in c(1:length(nop))) {
  if(file.exists(nop[count])) {
    print(paste(noprimer_filenames[count], " already exists. Skipping removing primers step...", sep=""))
  } else {
    print("Removing primers from PacBio...")
    prim <- removePrimers(PacBio_fns, nop, primer.fwd=tprKF, primer.rev=rc(tprKR), orient=TRUE, verbose=TRUE)
  }
}
  

print("Filtering PacBio reads...")
## Setting up file names to filter.
filter_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.fastq",sep ='')
filterEE1_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.EE1.fastq",sep ='')
filt <- file.path(path,filter_filenames)
filtEE1 <- file.path(path,filterEE1_filenames)

##Filter reads for tprK length and do not worry about expected errors.
for (count in c(1:length(filt))) {
  if (file.exists(filt[count])) {
    print(paste(filter_filenames[count]," already exists. Skipping filtering step..."), sep="")
    } else {
      print(paste("Filtering ",nop[count],"...",sep=""))
      track <- fastqFilter(nop[count], filt[count], minLen=1400,maxLen=1800,
                       maxN=0,
                       compress=FALSE, multithread=TRUE)
    }
}

##Consider: Filter reads for tprK length and allow only 1 expected error for the entire read.
# for (count in c(1:length(filtEE1))) {
#   track <- fastqFilter(nop[count], filtEE1[count], minLen=1400,maxLen=1800,
#                        maxN=0, maxEE=1,
#                        compress=FALSE, multithread=TRUE)
# }

RAD_filenames <- paste(substr(basename(PacBio_fns),1,nchar(basename(PacBio_fns))-5),"noprimers.filtered.RAD.fasta",sep ='')
RAD_files <- file.path(path, RAD_filenames)

print("Constructing RAD files...")
## Build RAD files for each PacBio sample. This step takes forever!!!
for (count in c(1:length(filt))) {
  to_rad_name <- paste(RAD_filenames[count])
  # Skips RAD step if files already exist, because it takes forever.
  if(file.exists(to_rad_name)) {
    print(paste(to_rad_name, " already exists. Skipping RAD step...", sep=""))
  } else{
    ## mind the path. In a typical Mac installation, it points to the Julia application in the Application folder
    print("Setting up Julia...")
    julia <- julia_setup(JULIA_HOME = "/Applications/Julia-1.2.app/Contents/Resources/julia/bin/")
    julia_command("using Pkg")
    julia_command("Pkg.add(PackageSpec(name=\"NextGenSeqUtils\", rev= \"1.0\", url = \"https://github.com/MurrellGroup/NextGenSeqUtils.jl.git\"))")
    julia_command("Pkg.add(PackageSpec(name=\"DPMeansClustering\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/DPMeansClustering.jl.git\"))")
    julia_command("Pkg.add(PackageSpec(name=\"RobustAmpliconDenoising\", rev=\"1.0\", url = \"https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git\"))")
    julia_command("using RobustAmpliconDenoising")
    
    julia_readfastq <- paste("seqs, QVs, seq_names = read_fastq(\"",filt[count],'")',sep="")
    julia_command(julia_readfastq)
    julia_command("templates,template_sizes,template_indices = denoise(seqs)")
    julia_writefasta <- paste("write_fasta(\"",RAD_files[count],'",templates,names = ["seqs$(j)_$(template_sizes[j])" for j in 1:length(template_sizes)])',sep="")
    julia_command(julia_writefasta)
  }
}

## RAD denoised files are written.  Let's get some frequencies of different variable regions
RAD_files_nolines <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fasta",sep ='')
RAD_files_fix <- paste(substr(RAD_files,1,nchar(RAD_files)-5),"nolines.fix.fasta",sep ='')
PacBio_freq_path <- paste(path,"/PacBio_frequencies",sep='')
RAD_files_fix_newdir <- file.path(PacBio_freq_path,basename(RAD_files_fix))
mkdir_freq <- paste("mkdir ",PacBio_freq_path,sep='')
system(mkdir_freq)
copyPacBio <- paste("cp ",PacBio_fns,PacBio_freq_path,";")
for (num in 1:length(copyPacBio)) {
  system(copyPacBio[num])
}

# Fixes up the fastas so they wrap and don't have awkward new lines.
awk_command <- paste("awk '/^>/ {printf(\"\\n%s\\n\",$0);next; } { printf(\"%s\",$0);}  END {printf(\"\\n\");}' < ",RAD_files," > ",RAD_files_nolines," ;")
fix_firstline <- paste("tail -n +2 ",RAD_files_nolines," > ",RAD_files_fix_newdir)
for (count in c(1:length(awk_command))) {
  system(awk_command[count])
  system(fix_firstline[count])
}

print("Making PacBio frequency files...")
# Calls syph_r to make the final_data.csv for each PacBio sample.
syph_path <- paste(script.dir,"/syph_r.py",sep='')
for (num in c(1:length(PacBio_freq_path))) {
  syphrPacBio_command <- paste("/Library/Frameworks/Python.framework/Versions/3.7/bin/python3 ",syph_path,
                               " -i fasta -pacbio -d ",PacBio_freq_path[num],"/",sep='')
  system(syphrPacBio_command)
}

## Make PacBio frequency comparison file
PacBio_freq_files = list.files(PacBio_freq_path,pattern="*_final_data.csv")
PacBio_freq_files_fullpath = list.files(PacBio_freq_path,pattern="*_final_data.csv",full.names=T)
compare_PacBio_df <- data.frame(Region=character(),Read=character())

for (i in 1:length(PacBio_freq_files)) {
  PacBioFreqtitle <- paste("PB_",sample_names[i],"_RelativeFreq",sep='')
  PacBioCounttitle <- paste("PB_",sample_names[i],"_Count",sep='')
  assign(PacBio_freq_files[i], read.csv(PacBio_freq_files_fullpath[i],col.names = c("Region","Read",PacBioFreqtitle,PacBioCounttitle)))
  compare_PacBio_df <- merge(compare_PacBio_df,get(PacBio_freq_files[i]),all=TRUE)
} 

if (opt$pacbio) {
  print("Pacbio option specified. Skipping making Illumina frequency files...")
} else {
  print("Making Illumina frequency files...")
  Illumina_freq_path <- paste(path,"/Illumina_frequencies",sep='')
  mkdir_freq <- paste("mkdir ",Illumina_freq_path,sep='')
  system(mkdir_freq)
  copyIllumina <- paste("cp ",Illumina_fns,Illumina_freq_path,";")
  for (num in 1:length(copyIllumina)) {
    system(copyIllumina[num])
  }
  syphrIllumina_command <- paste("/Library/Frameworks/Python.framework/Versions/3.7/bin/python3 ",syph_path," -i fastq -illumina -d ",Illumina_freq_path,"/ ;",sep='')
  system(syphrIllumina_command)
}

if (opt$pacbio) {
  print("Pacbio option specified. Skipping making Illumina frequency comparison files...")
  allreads <- compare_PacBio_df
  allreads_out <- paste(path,"/allreads.csv",sep='')
  write.csv(allreads,file=allreads_out,row.names=FALSE,quote=FALSE)
} else {
  ##Make Illumina frequency comparison file
  Illumina_freq_files = list.files(Illumina_freq_path,pattern="*_final_data.csv")
  Illumina_freq_files_fullpath = list.files(Illumina_freq_path,pattern="*_final_data.csv",full.names=T)
  
  compare_Illumina_df <- data.frame(Region=character(),Read=character())
  #svgIllumina_command <- paste("/usr/bin/python /Users/uwvirongs/Documents/Michelle/syphilis/syph_visualizer_single_svg.py ",Illumina_freq_files_fullpath,sep='')
  #for (num in length(svgIllumina_command)) system(svgIllumina_command[num])
  
  for (i in 1:length(Illumina_freq_files)) {
    IlluminaFreqtitle <- paste("Ill_",sample_names[i],"_RelativeFreq",sep='')
    IlluminaCounttitle <- paste("Ill_",sample_names[i],"_Count",sep='')
    Illuminadf <- read.csv(Illumina_freq_files_fullpath[i],col.names = c("Region","Read",IlluminaFreqtitle,IlluminaCounttitle),check.names = FALSE)
    Illuminadf <- Illuminadf[order(Illuminadf$Region,-Illuminadf[[IlluminaFreqtitle]]),][]
    compare_Illumina_df <- merge(compare_Illumina_df,Illuminadf,all=TRUE)
  }
  
  # Merges the PacBio and Illumina comparison tables into the all-important allreads.csv!
  print("Merging to allreads.csv...")
  allreads <- merge(compare_PacBio_df,compare_Illumina_df,all=T)
  allreads_out <- paste(path,"/allreads.csv",sep='')
  write.csv(allreads,file=allreads_out,row.names=FALSE,quote=FALSE)
}


# allreads$`Illumina
# allreads$`IlluminaRelativeFreq_148B-tprK` <- as.numeric(allreads$`IlluminaRelativeFreq_148B-tprK`)
# allreads$`IlluminaRelativeFreq_148B2-tprK` <- as.numeric(allreads$`IlluminaRelativeFreq_148B2-tprK`)
# ggplot(allreads,aes(x=X148B_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + xlim(0,10) + ylim(0,10)
# 
# ggplot(allreads,aes(x=X148B_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + scale_y_log10(limits = c(0.05,100)) + scale_x_log10(limits = c(0.05,100)) #+ xlim(0.01,100) + ylim(0.01,100)
# 
# ggplot(allreads,aes(x=allreads$X148B2.tprK_PacBioRelativeFreq,y=`IlluminaRelativeFreq_148B2-tprK`,color=Region)) + geom_point(size=2) + geom_abline(linetype="dashed",color="grey") +
#   theme_classic() + xlim(0,10) + ylim(0,10)