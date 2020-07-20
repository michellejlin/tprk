# Generates dot-line plots for comparing the variable regions between two samples.
# Currently takes the path and goes through the allreads.csv.

list.of.packages <- c("ggplot2", "grid", "nplr", "plyr", "dplyr", "scales", "gridExtra", "RColorBrewer", "optparse","randomcoloR", "cowplot",
                      "tidyr", "tibble")
lapply(list.of.packages,library,character.only=T)

option_list <- list(make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Specify metadata", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

#path <- opt$directory
path <- "."

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding line of path
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are.

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files"

## This script can also be run from the command line.
## Usage: rscript \path\to\Variable_region_compare.R -p [path]

#####

allreads <- paste(path,"/allreads_filtered.csv", sep="")
allreads_filtered <- paste(path,"/allreads_filtered.csv", sep="")
# Grabs sample names from metadata.csv.
metadata <- read.table(opt$metadata, sep=',', header=TRUE)
sample_names = as.character(metadata$SampleName)
  
alldata <- read.csv(allreads,header=TRUE,sep=",",stringsAsFactors = FALSE)

# Loops through and generates variable region comparisons for all the combinations of files.
# This means a lot of files if list is long.
for (i in 1:(length(sample_names) - 1)) {
  for (j in (i+1):(length(sample_names))) {
    print(paste("Generating figure for ",sample_names[i]," and ",sample_names[j],"...",sep=""))
    rfcol <- paste("Ill_",sample_names[i],"_RelativeFreq",sep = "")
    rfcol2 <- paste("Ill_",sample_names[j],"_RelativeFreq",sep = "")

    commondfIllumina <- select(alldata,Region,Read,rfcol,rfcol2)
    commondfIllumina <- filter(commondfIllumina,!((commondfIllumina[[rfcol]] == 0) & (commondfIllumina[[rfcol2]] == 0)))
    sortedIllumina <- commondfIllumina[order(commondfIllumina$Region,-commondfIllumina[[rfcol]],-commondfIllumina[[rfcol2]]),]
    sortedIllumina <- gather(sortedIllumina,rfcol,rfcol2,key="Sample",value="Frequency")
    sortedIllumina$Sample[sortedIllumina$Sample == rfcol] <- sample_names[i]
    sortedIllumina$Sample[sortedIllumina$Sample == rfcol2] <- sample_names[j]
    
    myColors <- distinctColorPalette(1068)
    names(myColors) <- levels(sortedIllumina$Read)
    colScale <- scale_colour_manual(name = NULL, guide = FALSE, values = myColors)
    h <- ggplot(sortedIllumina[which(sortedIllumina$Frequency>0),]) + geom_point(aes(y = Frequency, x = Sample, color=Read)) + geom_line(aes(y = Frequency, x = Sample, group=Read, color=Read)) +  
      facet_wrap(~Region,nrow=1) + theme(legend.position = "none") + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
      theme(axis.text.x=element_text(angle=45,hjust=1))
    h1 <- h + colScale
    suppressMessages(ggsave(paste(sample_names[i],"_vs_",sample_names[j],"_VariableRegions_DotLine_Filtered.pdf",sep=""),
           path="./",plot=h1,width=5,height=4,units="in"))
  }
}
